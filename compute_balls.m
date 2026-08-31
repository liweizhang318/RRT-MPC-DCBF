function compute_balls
% COMPUTE_BALLS  Compute the dynamic-obstacle parameters for visualize.m so that
%   ONE ball travels from afar, crosses the robot's path exactly at the point
%   (0.5,-0.5), and the robot must detour around it. Prints the computed dyn_obs row.
%   (Runs the same simulation as visualize.m but WITHOUT dynamic obstacles.)

    % ---- parameters (same as visualize.m) ----
    DT = 0.1; N = 12; GAMMA = 0.3;
    V_MAX = 1.0; V_MIN = 0.0; OMEGA_MAX = 1.0;
    ROBOT_R = 0.15; SAFE_D = 0.25; GOAL_R = 0.25; V_REF = 0.8;
    Q_XY = 10; Q_TH = 0.1; R_V = 0.01; R_W = 0.02; P_XY = 50;
    static_obs = [-2.0  1.5 0.6;
                   1.0  2.0 0.7;
                   2.5 -1.0 0.5;
                  -1.5 -2.5 0.6];
    start = [-4 -4 0]; goal = [4 4];
    r_eff_all = static_obs(:,3)' + ROBOT_R + SAFE_D;

    % ---- RRT + shortcut ----
    path = rrt(start(1:2), goal, static_obs);
    path = shortcut(path, static_obs);

    % ---- simulation WITHOUT dynamic obstacles; record trajectory + time ----
    curr = start; t = 0; last_z = [];
    traj = start(1:2); times = 0;
    for step = 1:600
        if norm(curr(1:2) - goal) < GOAL_R, break; end
        ref = make_ref(path, curr, goal, N, DT, V_REF);
        obs_hor = obstacle_horizon(static_obs, t, N, DT);
        [u, last_z] = solve_mpc(curr, ref, obs_hor, last_z, ...
            N, DT, V_MAX, V_MIN, OMEGA_MAX, V_REF, GAMMA, r_eff_all, ...
            Q_XY, Q_TH, R_V, R_W, P_XY);
        curr = curr + [u(1)*cos(curr(3))*DT, u(1)*sin(curr(3))*DT, u(2)*DT];
        t = t + DT;
        traj(end+1,:) = curr(1:2); times(end+1) = t;
    end
    fprintf('No-dyn-ball run: reached goal @ t=%.1fs, %d steps\n', t, step);

    % ---- fixed encounter point: where the removed static obstacle used to be ----
    Pc = [0.5, -0.5];
    [~, ic] = min(sum((traj - Pc).^2, 2));
    t_cross = times(ic);
    fprintf('Robot passes closest to (0.5,-0.5) at t=%.1fs, pos=(%.2f,%.2f)\n', ...
            t_cross, traj(ic,1), traj(ic,2));

    % ---- ball: speed chosen freely, direction = normal to local path ----
    v_ball = 0.50;
    n = normal_at(traj, ic);
    % choose the normal orientation whose start point lies farther outside the map edge
    s1 = Pc - n * v_ball * t_cross;
    s2 = Pc + n * v_ball * t_cross;
    if min(abs(s1)) < min(abs(s2))   % prefer the start that is closer to the map border
        s = s2; v = -n * v_ball;
    else
        s = s1; v = n * v_ball;
    end
    fprintf('\ndyn_obs row (radius 0.45): [%.2f %.2f 0.45 %.2f %.2f]\n', s(1), s(2), v(1), v(2));
end

function n = normal_at(traj, i)
    i0 = max(1, i-5); i1 = min(size(traj,1), i+5);
    d = traj(i1,:) - traj(i0,:);
    if norm(d) < 1e-9, d = [1 0]; end
    n = [-d(2), d(1)] / norm(d);
end

% ==================== helpers (same as visualize.m) ====================
function p = rrt(s, g, obs)
    rng(0);
    nodes = s; parent = 0;
    for it = 1:4000
        if rand < 0.1, samp = g; else, samp = [rand*10-5, rand*10-5]; end
        [~, idx] = min(vecnorm(nodes - samp, 2, 2));
        near = nodes(idx,:); vec = samp - near; d = norm(vec);
        if d < 1e-6, continue; end
        new = near + vec/d * min(0.5, d);
        if ~seg_free(near, new, obs), continue; end
        nodes = [nodes; new]; parent = [parent; idx];
        if norm(new - g) < 0.25   % GOAL_R: 0.25, same value as visualize.m
            p = [g]; cur = size(nodes,1);
            while cur > 0, p = [nodes(cur,:); p]; cur = parent(cur); end
            return;
        end
    end
    p = [];
end

function ok = seg_free(a, b, obs)
    d = norm(b-a); n = max(1, round(d/0.1));
    for i = 0:n
        p = a + (b-a)*i/n;
        for j = 1:size(obs,1)
            if norm(p - obs(j,1:2)) < obs(j,3) + 0.15 + 0.25
                ok = false; return;
            end
        end
    end
    ok = true;
end

function p = shortcut(p, obs)
    out = p(1,:); i = 1; n = size(p,1);
    while i < n
        next_i = i + 1;
        for j = n:-1:i+1
            if seg_free(p(i,:), p(j,:), obs), next_i = j; break; end
        end
        out = [out; p(next_i,:)]; i = next_i;
    end
    p = out;
end

function ref = make_ref(p, c, g, N, DT, V_REF)
    full = [p; g];
    best_i = 1; best_proj = full(1,:); best_d = inf;
    for i = 1:size(full,1)-1
        a = full(i,:); b = full(i+1,:); ab = b-a;
        L2 = dot(ab,ab);
        if L2 < 1e-9, continue; end
        tt = min(max(dot(c(1:2)-a, ab)/L2, 0), 1);
        proj = a + tt*ab; dd = norm(c(1:2)-proj);
        if dd < best_d, best_d = dd; best_i = i; best_proj = proj; end
    end
    pts = [best_proj; full(best_i+1:end,:)];
    cum = zeros(size(pts,1),1);
    for i = 2:size(pts,1), cum(i) = cum(i-1) + norm(pts(i,:)-pts(i-1,:)); end
    ref = zeros(N+1,3); ref(1,:) = c;
    for k = 1:N
        target = V_REF*DT*k;
        i = find(cum >= target, 1) - 1;
        if isempty(i) || i < 1, i = 1; end
        if i > size(pts,1)-1, i = size(pts,1)-1; end
        seg = pts(i+1,:) - pts(i,:); L = norm(seg);
        if L < 1e-9, pt = pts(i,:); else
            frac = min(max((target-cum(i))/(cum(i+1)-cum(i)+1e-9),0),1);
            pt = pts(i,:) + seg*frac;
        end
        ref(k+1,1:2) = pt;
        ref(k+1,3) = atan2(ref(k+1,2)-ref(k,2), ref(k+1,1)-ref(k,1));
    end
end

function obs_hor = obstacle_horizon(sobs, tt, N, DT)
    n_s = size(sobs,1);
    obs_hor = zeros(n_s, N, 2);
    for j = 1:n_s
        for k = 1:N, obs_hor(j,k,:) = sobs(j,1:2); end
    end
end

function [u, z] = solve_mpc(c, ref, obs_hor, last_z, N, DT, V_MAX, V_MIN, OMEGA_MAX, V_REF, GAMMA, r_eff_all, Q_XY, Q_TH, R_V, R_W, P_XY) %#ok<INUSD>
    dth = atan2(sin(ref(2,3)-c(3)), cos(ref(2,3)-c(3)));
    w0 = min(max(dth/DT, -OMEGA_MAX), OMEGA_MAX);
    lb = [repmat(V_MIN,N,1); repmat(-OMEGA_MAX,N,1)];
    ub = [repmat(V_MAX,N,1); repmat( OMEGA_MAX,N,1)];
    z0 = [repmat(V_REF,N,1); repmat(w0,N,1)];
    opts = optimoptions('fmincon', 'Display','off', 'Algorithm','sqp', ...
                        'MaxIterations', 150, 'OptimalityTolerance', 1e-3, ...
                        'ConstraintTolerance', 1e-3, 'StepTolerance', 1e-6);
    [z, ~, ef] = fmincon(@(z) mpc_cost(z, c, ref, N, DT, Q_XY, Q_TH, R_V, R_W, P_XY), ...
                         z0, [], [], [], [], lb, ub, ...
                         @(z) mpc_nonlcon(z, c, obs_hor, N, DT, GAMMA, r_eff_all), opts);
    if ef == -2, z = [zeros(N,1); zeros(N,1)]; end
    u = [z(1); z(N+1)];
end

function J = mpc_cost(z, c, ref, N, DT, Q_XY, Q_TH, R_V, R_W, P_XY)
    v = z(1:N); w = z(N+1:2*N);
    x = c; J = 0;
    for k = 1:N
        x = [x(1)+v(k)*cos(x(3))*DT, x(2)+v(k)*sin(x(3))*DT, x(3)+w(k)*DT];
        ex = x(1)-ref(k+1,1); ey = x(2)-ref(k+1,2);
        eth = atan2(sin(x(3)-ref(k+1,3)), cos(x(3)-ref(k+1,3)));
        J = J + Q_XY*(ex^2+ey^2) + Q_TH*eth^2 + R_V*v(k)^2 + R_W*w(k)^2;
    end
    ex = x(1)-ref(N+1,1); ey = x(2)-ref(N+1,2);
    J = J + P_XY*(ex^2+ey^2);
end

function [cc, ceq] = mpc_nonlcon(z, c, obs_hor, N, DT, GAMMA, r_eff_all)
    v = z(1:N); w = z(N+1:2*N);
    n_o = size(obs_hor,1);
    h_prev = zeros(n_o,1);
    for j = 1:n_o
        h_prev(j) = norm([c(1)-obs_hor(j,1,1), c(2)-obs_hor(j,1,2)]) - r_eff_all(j);
    end
    x = c; cc = [];
    for k = 1:N
        x = [x(1)+v(k)*cos(x(3))*DT, x(2)+v(k)*sin(x(3))*DT, x(3)+w(k)*DT];
        ob_idx = min(k+1, N);
        h_cur = zeros(n_o,1);
        for j = 1:n_o
            h_cur(j) = norm([x(1)-obs_hor(j,ob_idx,1), x(2)-obs_hor(j,ob_idx,2)]) - r_eff_all(j);
        end
        cc = [cc; (1-GAMMA)*h_prev - h_cur]; %#ok<AGROW>
        cc = [cc; -h_cur]; %#ok<AGROW>
        h_prev = h_cur;
    end
    ceq = [];
end
