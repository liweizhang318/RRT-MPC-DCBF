function test_suite(n_runs, include_headon)
% TEST_SUITE  Multi-scenario stress test of the pure RRT + MPC-D-CBF framework.
%   Random obstacle layouts (static + dynamic; dynamic obstacles in random
%   directions: head-on / crossing / same-direction), one full simulation per
%   scenario, statistics: arrival rate, collision rate, stuck rate, min surface
%   distance.
%
%   Usage: test_suite(50)            run 50 random scenarios (with head-on obstacles)
%          test_suite(50, false)     without head-on obstacles (control experiment)

    if nargin < 1, n_runs = 20; end
    if nargin < 2, include_headon = true; end

    % ---- fixed parameters (same as simulate.m) ----
    DT = 0.1; N = 12; GAMMA = 0.3;
    V_MAX = 1.0; V_MIN = 0.0; OMEGA_MAX = 1.0;
    ROBOT_R = 0.15; SAFE_D = 0.15; GOAL_R = 0.25; V_REF = 0.8;
    Q_XY = 10; Q_TH = 0.1; R_V = 0.01; R_W = 0.02; P_XY = 50;

    n_arrive = 0; n_collision = 0; n_stuck = 0;
    min_dists = []; steps_arr = [];

    for run = 1:n_runs
        rng(run);                       % fixed seed per run for reproducibility
        [sobs, dobs, start, goal] = gen_scene(include_headon);

        r_eff_all = [sobs(:,3)' + ROBOT_R + SAFE_D, dobs(:,3)' + ROBOT_R + SAFE_D];

        % RRT
        path = rrt(start(1:2), goal, sobs, ROBOT_R, SAFE_D, GOAL_R);
        if isempty(path)
            fprintf('[run %2d] RRT failed (no path found)\n', run);
            continue;
        end
        path = shortcut(path, sobs, ROBOT_R, SAFE_D);

        % simulation
        curr = start; t = 0; last_z = []; min_dist = inf;
        collided = false; arrived = false;
        for step = 1:600
            if norm(curr(1:2) - goal) < GOAL_R
                arrived = true; break;
            end
            ref = make_ref(path, curr, goal, N, DT, V_REF);
            obs_hor = obstacle_horizon(sobs, dobs, t, N, DT);
            [u, last_z] = solve_mpc(curr, ref, obs_hor, last_z, N, DT, ...
                V_MAX, V_MIN, OMEGA_MAX, V_REF, GAMMA, r_eff_all, ...
                Q_XY, Q_TH, R_V, R_W, P_XY);
            curr = curr + [u(1)*cos(curr(3))*DT, u(1)*sin(curr(3))*DT, u(2)*DT];
            t = t + DT;
            % collision detection
            for i = 1:size(sobs,1)
                d = norm(curr(1:2)-sobs(i,1:2)) - sobs(i,3) - ROBOT_R;
                if d < 0, collided = true; end
                min_dist = min(min_dist, d);
            end
            for i = 1:size(dobs,1)
                x = dobs(i,1)+dobs(i,4)*t; y = dobs(i,2)+dobs(i,5)*t;
                d = norm(curr(1:2)-[x y]) - dobs(i,3) - ROBOT_R;
                if d < 0, collided = true; end
                min_dist = min(min_dist, d);
            end
        end

        if collided, n_collision = n_collision + 1; end
        if arrived
            n_arrive = n_arrive + 1;
        else
            n_stuck = n_stuck + 1;
        end
        min_dists = [min_dists, min_dist]; %#ok<AGROW>
        if arrived, steps_arr = [steps_arr, step]; end %#ok<AGROW>

        fprintf('[run %2d] %s  %s  min surface dist=%.3fm\n', run, ...
                ternary(arrived, 'arrived', 'stuck'), ...
                ternary(collided, 'collision!', 'safe'), min_dist);
    end

    fprintf('\n========== Summary (%d scenarios) ==========\n', n_runs);
    fprintf('Arrived: %d / %d (%.0f%%)\n', n_arrive, n_runs, 100*n_arrive/n_runs);
    fprintf('Collision: %d / %d (%.0f%%)\n', n_collision, n_runs, 100*n_collision/n_runs);
    fprintf('Stuck: %d / %d (%.0f%%)\n', n_stuck, n_runs, 100*n_stuck/n_runs);
    if ~isempty(min_dists)
        fprintf('Min surface dist: worst=%.3fm  median=%.3fm\n', min(min_dists), median(min_dists));
    end
    if ~isempty(steps_arr)
        fprintf('Avg steps to arrive: %.0f (%.1fs)\n', mean(steps_arr), mean(steps_arr)*DT);
    end

    function s = ternary(cond, a, b)
        if cond, s = a; else, s = b; end
    end
end

% ==================== Scenario generation ====================
function [sobs, dobs, start, goal] = gen_scene(include_headon)
    start = [-4 -4 0];
    goal  = [4 4];
    % static obstacles: 5 random circles, clear of start/goal and the central corridor
    sobs = zeros(5,3);
    for i = 1:5
        for try_k = 1:200
            x = rand*8 - 4; y = rand*8 - 4; r = 0.4 + rand*0.4;
            if norm([x y]-start(1:2)) < 1.2 || norm([x y]-goal) < 1.2, continue; end
            % do not completely block the diagonal corridor
            if abs(x - y) < 0.5 && abs(x) < 3 && abs(y) < 3, continue; end
            sobs(i,:) = [x y r]; break;
        end
        if sobs(i,3) == 0, sobs(i,:) = [10 10 0.3]; end  % if generation failed, move far away
    end
    % dynamic obstacles: 2, random start/direction, covering head-on / crossing / same-direction
    dobs = zeros(2,5);
    for i = 1:2
        r = 0.4 + rand*0.15;
        % random motion mode (head-on replaced by same-direction drift when include_headon=false)
        mode = randi(3);
        if mode == 1        % crossing (perpendicular through the main corridor)
            x = (rand*2-1)*3; y = (rand*2-1)*3;
            vx = 0; vy = sign(0 - y)* (0.3 + rand*0.3);   % toward the center
            if abs(vy) < 0.05, vy = 0.35; end
        elseif mode == 2 && include_headon    % head-on (from goal direction toward start)
            x = 3 + rand; y = 3 + rand*2 - 1;
            vx = -0.25 - rand*0.2; vy = -0.1 - rand*0.15;
        else                % same-direction / oblique drift
            x = (rand*2-1)*3.5; y = (rand*2-1)*3.5;
            vx = 0.2*rand; vy = 0.2*rand;
        end
        dobs(i,:) = [x y r vx vy];
    end
end

% ==================== RRT ====================
function p = rrt(s, g, obs, ROBOT_R, SAFE_D, GOAL_R)
    rng(0);
    nodes = s; parent = 0;
    for it = 1:4000
        if rand < 0.1, samp = g; else, samp = [rand*10-5, rand*10-5]; end
        [~, idx] = min(vecnorm(nodes - samp, 2, 2));
        near = nodes(idx,:);
        vec = samp - near; d = norm(vec);
        if d < 1e-6, continue; end
        new = near + vec/d * min(0.5, d);
        if ~seg_free(near, new, obs, ROBOT_R, SAFE_D), continue; end
        nodes = [nodes; new]; parent = [parent; idx];
        if norm(new - g) < GOAL_R
            p = [g]; cur = size(nodes,1);
            while cur > 0, p = [nodes(cur,:); p]; cur = parent(cur); end
            return;
        end
    end
    p = [];
end

function ok = seg_free(a, b, obs, ROBOT_R, SAFE_D)
    d = norm(b-a); n = max(1, round(d/0.1));
    for i = 0:n
        p = a + (b-a)*i/n;
        for j = 1:size(obs,1)
            if norm(p - obs(j,1:2)) < obs(j,3) + ROBOT_R + SAFE_D
                ok = false; return;
            end
        end
    end
    ok = true;
end

function p = shortcut(p, obs, ROBOT_R, SAFE_D)
    out = p(1,:); i = 1; n = size(p,1);
    while i < n
        next_i = i + 1;
        for j = n:-1:i+1
            if seg_free(p(i,:), p(j,:), obs, ROBOT_R, SAFE_D)
                next_i = j; break;
            end
        end
        out = [out; p(next_i,:)]; i = next_i;
    end
    p = out;
end

% ==================== Reference trajectory ====================
function ref = make_ref(p, c, g, N, DT, V_REF)
    full = [p; g];
    best_i = 1; best_proj = full(1,:); best_d = inf;
    for i = 1:size(full,1)-1
        a = full(i,:); b = full(i+1,:); ab = b-a;
        L2 = dot(ab,ab);
        if L2 < 1e-9, continue; end
        tt = min(max(dot(c(1:2)-a, ab)/L2, 0), 1);
        proj = a + tt*ab;
        dd = norm(c(1:2)-proj);
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

% ==================== Obstacle prediction horizon ====================
function obs_hor = obstacle_horizon(sobs, dobs, tt, N, DT)
    n_s = size(sobs,1); n_d = size(dobs,1); n_o = n_s + n_d;
    obs_hor = zeros(n_o, N, 2);
    for j = 1:n_s
        for k = 1:N, obs_hor(j,k,:) = sobs(j,1:2); end
    end
    for j = 1:n_d
        for k = 1:N
            obs_hor(n_s+j,k,:) = [dobs(j,1)+dobs(j,4)*(tt+(k-1)*DT), ...
                                  dobs(j,2)+dobs(j,5)*(tt+(k-1)*DT)];
        end
    end
end

% ==================== MPC-D-CBF solver ====================
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
    x = c;
    cc = [];
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
