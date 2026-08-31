function simulate(animate)
% SIMULATE  RRT + MPC + D-CBF 2D real-time simulation animation (MATLAB + fmincon)
%   Robot: unicycle [x, y, theta], control [v, omega]
%   Obstacles: circles; static ones stay fixed; dynamic balls (same size as the
%     robot) are placed at random positions along the path and cross it back and
%     forth perpendicularly (random spacing, staggered phases)
%   Safety: D-CBF decreasing chain  h(x_{k+1}, ob_{k+1}) >= (1-gamma) h(x_k, ob_k)
%     plus h >= 0 hard constraints
%   Goal: a sequence of random goals; upon reaching one, sample the next
%     (pure RRT + MPC + D-CBF)
%
%   Two-stage pipeline:
%     Stage 1 (offline): run the full simulation, store per-step states only,
%       no plotting
%     Stage 2 (playback): replay from stored data quickly; frame rate is
%       independent of solver speed
%
%   Usage: simulate               (compute first, then playback animation)
%          simulate(false)        (compute only, no plotting)
%   (The playback window has an "Export video" button that saves the current
%     run to MP4)
%   (Requires fmincon from the Optimization Toolbox, bundled with MATLAB)

    if nargin < 1, animate = true; end

    % ============ Parameters ============
    DT = 0.1;                 % control period (s)
    N  = 12;                  % MPC prediction horizon
    GAMMA = 0.3;              % D-CBF decay rate
    V_MAX = 1.0; V_MIN = 0.0; % linear speed [0, 1] (no reversing)
    OMEGA_MAX = 1.0;          % angular speed limit
    ROBOT_R = 0.15;           % robot radius
    SAFE_D  = 0.10;           % extra safety margin
    GOAL_R  = 0.25;           % arrival detection radius
    V_REF   = 0.8;            % reference linear speed
    N_GOALS = 5;              % number of consecutive random goals
    % cost weights
    Q_XY = 10; Q_TH = 0.1; R_V = 0.01; R_W = 0.02; P_XY = 50;

    % ============ Static obstacles ============
    static_obs = [-2.0  1.5 0.6;
                   1.0  2.0 0.7;
                   2.5 -1.0 0.5;
                  -1.5 -2.5 0.6;
                   0.5 -0.5 0.4;
                  -3.5  0.5 0.5;
                   3.0  1.5 0.55;
                   0.0  3.0 0.55;
                  -2.5 -0.5 0.45;
                   1.5 -3.0 0.5];

    % ============ Dynamic obstacles: balls (same size as the robot) constantly
    % crossing the path perpendicularly at random positions ============
    n_dyn  = 7;                                     % number of dynamic balls (smaller count keeps the MPC constraint set lighter)
    dyn_r  = ROBOT_R;                               % radius = 0.15, same as the robot
    dyn_A  = 1.2;                                   % crossing amplitude (+-1.2 m: balls turn around closer to the path, crossing it more often)
    dyn_v  = 0.30;                                  % crossing speed m/s (quick crossings reduce time spent blocking the path)

    start = [-4 -4 0];

    % CBF effective radius (obstacle radius + robot radius + margin)
    r_eff_all = [static_obs(:,3)' + ROBOT_R + SAFE_D, ...
                 (dyn_r + ROBOT_R + SAFE_D)*ones(1,n_dyn)];

    % ============ Stage 1: offline computation (store data only, no plotting) ============
    curr = start;
    t = 0;
    last_z = [];
    min_dist = inf;
    collision_printed = false;

    rng('shuffle');           % different goal sequence + ball positions each run

    % random ball positions (random spacing) + random phases (staggered crossings)
    dyn_s   = sort(0.05 + 0.90*rand(n_dyn,1));
    dyn_phi = rand(n_dyn,1) * (4*dyn_A/dyn_v);

    dyn_n = zeros(n_dyn,2); dyn_c = zeros(n_dyn,2); t_path = 0;

    % recording arrays (for playback)
    rec_traj = zeros(0,2);        % robot position at each step
    rec_dyn  = zeros(0,2*n_dyn);  % dynamic obstacle positions at each step (flattened)
    rec_g    = zeros(0,1);        % goal index at each step
    rec_t    = zeros(0,1);        % time at each step
    rec_paths  = {};              % path for each goal
    rec_goals  = zeros(0,2);      % goal position for each goal
    rec_dyn_c  = {};              % crossing center per goal
    rec_dyn_n  = {};              % crossing direction per goal (per-ball, n_dyn x 2)
    rec_t_path = zeros(0,1);      % t_path per goal

    n_reached = 0;
    goal_idx = 0;
    for g = 1:N_GOALS
        % ---- sample the current goal ----
        goal = sample_goal(curr, static_obs, ROBOT_R, SAFE_D, GOAL_R);
        path = rrt(curr(1:2), goal, static_obs, ROBOT_R, SAFE_D, GOAL_R);
        attempt = 0;
        while isempty(path) && attempt < 10
            goal = sample_goal(curr, static_obs, ROBOT_R, SAFE_D, GOAL_R);
            path = rrt(curr(1:2), goal, static_obs, ROBOT_R, SAFE_D, GOAL_R);
            attempt = attempt + 1;
        end
        if isempty(path)
            fprintf('goal %d/%d: RRT found no path, skipping\n', g, N_GOALS);
            continue;
        end
        path = shortcut(path, static_obs, ROBOT_R, SAFE_D);
        goal_idx = goal_idx + 1;

        % ---- align dynamic-obstacle crossings (updated whenever the path is rebuilt, covering every polyline segment) ----
        % balls spread by cumulative arc length along the polyline; crossing
        % direction = perpendicular to the containing polyline segment
        n_seg = size(path,1) - 1;
        seg_len = zeros(n_seg,1);
        for ii = 1:n_seg
            seg_len(ii) = norm(path(ii+1,:) - path(ii,:));
        end
        total_len = sum(seg_len);
        cum_len = [0; cumsum(seg_len)];
        dyn_c = zeros(n_dyn,2);
        dyn_n = zeros(n_dyn,2);
        for j = 1:n_dyn
            s = dyn_s(j) * total_len;                % arc-length position of this ball along the polyline
            idx = find(cum_len >= s, 1) - 1;
            if idx < 1, idx = 1; end
            if idx > n_seg, idx = n_seg; end
            seg = path(idx+1,:) - path(idx,:);
            L = norm(seg);
            if L < 1e-9
                frac = 0;
            else
                frac = min(max((s - cum_len(idx)) / L, 0), 1);
            end
            dyn_c(j,:) = path(idx,:) + seg * frac;
            th = atan2(seg(2), seg(1));
            dyn_n(j,:) = [-sin(th), cos(th)];
        end
        t_path = t;

        % record goal data
        rec_paths{goal_idx} = path;         %#ok<AGROW>
        rec_goals(goal_idx,:) = goal;
        rec_dyn_c{goal_idx} = dyn_c;        %#ok<AGROW>
        rec_dyn_n{goal_idx} = dyn_n;        %#ok<AGROW>
        rec_t_path(goal_idx) = t_path;

        % ---- inner loop: drive to the current goal (or give up if stuck / timeout) ----
        best_dist = inf; stall_cnt = 0;
        for step = 1:600
            if norm(curr(1:2) - goal) < GOAL_R
                n_reached = n_reached + 1;
                fprintf('goal %d/%d reached @ t=%.1fs pos=(%.2f,%.2f)\n', g, N_GOALS, t, curr(1), curr(2));
                break;
            end

            ref = make_ref(path, curr, goal, N, DT, V_REF);
            obs_hor = obstacle_horizon(static_obs, dyn_c, dyn_n, dyn_A, dyn_v, dyn_phi, t, t_path, N, DT);
            [u, last_z] = solve_mpc(curr, ref, obs_hor, last_z, N, DT, V_MAX, V_MIN, OMEGA_MAX, V_REF, GAMMA, r_eff_all, Q_XY, Q_TH, R_V, R_W, P_XY); %#ok<ASGLU>

            curr = curr + [u(1)*cos(curr(3))*DT, u(1)*sin(curr(3))*DT, u(2)*DT];
            t = t + DT;

            % stall detection
            d_goal = norm(curr(1:2) - goal);
            if d_goal < best_dist - 0.15
                best_dist = d_goal; stall_cnt = 0;
            else
                stall_cnt = stall_cnt + 1;
            end
            if stall_cnt >= 300
                fprintf('goal %d/%d stuck (no progress for %.0fs), giving up and resampling\n', g, N_GOALS, stall_cnt*DT);
                break;
            end

            % collision detection + record this step's dynamic obstacle positions
            dyn_row = zeros(1, 2*n_dyn);
            for i = 1:size(static_obs,1)
                d = norm(curr(1:2)-static_obs(i,1:2)) - static_obs(i,3) - ROBOT_R;
                if d < 0 && ~collision_printed
                    fprintf('  [collision] t=%.1fs pos=(%.2f,%.2f) static obstacle %d surface dist=%.3f\n', t, curr(1), curr(2), i, d);
                    collision_printed = true;
                end
                min_dist = min(min_dist, d);
            end
            for i = 1:n_dyn
                pj = dyn_pos(dyn_c(i,:), dyn_n(i,:), dyn_A, dyn_v, dyn_phi(i), t - t_path);
                dyn_row(2*i-1:2*i) = pj;
                d = norm(curr(1:2)-pj) - dyn_r - ROBOT_R;
                if d < 0 && ~collision_printed
                    fprintf('  [collision] t=%.1fs pos=(%.2f,%.2f) dynamic obstacle %d@(%.2f,%.2f) surface dist=%.3f\n', ...
                            t, curr(1), curr(2), i, pj(1), pj(2), d);
                    collision_printed = true;
                end
                min_dist = min(min_dist, d);
            end

            % record this step
            rec_traj(end+1,:) = curr(1:2); %#ok<AGROW>
            rec_dyn(end+1,:) = dyn_row;    %#ok<AGROW>
            rec_g(end+1) = goal_idx;       %#ok<AGROW>
            rec_t(end+1) = t;              %#ok<AGROW>

            if mod(step, 50) == 0
                fprintf('  goal %d/%d step %d: t=%.1fs pos=(%.2f,%.2f)\n', g, N_GOALS, step, t, curr(1), curr(2));
            end
        end
    end

    if min_dist < 0
        fprintf('Done @ t=%.1fs, reached %d/%d goals, min surface dist=%.3fm (collision!)\n', t, n_reached, N_GOALS, min_dist);
    else
        fprintf('Done @ t=%.1fs, reached %d/%d goals, min surface dist=%.3fm (safe)\n', t, n_reached, N_GOALS, min_dist);
    end

    % ============ Stage 2: playback animation (from stored data, no re-solving) ============
    if animate && ~isempty(rec_traj)
        playback(rec_traj, rec_dyn, rec_g, rec_t, rec_paths, rec_goals, ...
                 static_obs, dyn_r, n_dyn, N_GOALS, start);
    end
end

% ==================== Top-level helper functions ====================

function playback(rec_traj, rec_dyn, rec_g, rec_t, rec_paths, rec_goals, ...
                  static_obs, dyn_r, n_dyn, N_GOALS, start)
 
    % with "Replay" and "Export video" buttons
    figure('Name', 'RRT + MPC + D-CBF Playback', 'Position', [80 80 720 720], 'Color', 'w');
    cla; hold on; axis equal; grid on;
    xlim([-5 5]); ylim([-5 5]);
    ax = gca;

    % static obstacles (drawn once)
    for i = 1:size(static_obs,1)
        rectangle('Position', [static_obs(i,1)-static_obs(i,3), static_obs(i,2)-static_obs(i,3), ...
                  2*static_obs(i,3), 2*static_obs(i,3)], 'Curvature',[1 1], ...
                  'FaceColor',[0.72 0.72 0.72], 'EdgeColor','k');
    end
    % dynamic obstacle handles (Position updated)
    dyn_h = gobjects(n_dyn,1);
    for i = 1:n_dyn
        dyn_h(i) = rectangle('Position', [0 0 2*dyn_r 2*dyn_r], 'Curvature',[1 1], ...
                   'FaceColor',[1 0.65 0.25], 'EdgeColor','k');
    end
    path_h = plot(0,0,'g--','LineWidth',1.5);
    traj_h = plot(0,0,'r-','LineWidth',2);
    robot_h = plot(0,0,'ro','MarkerSize',9,'MarkerFaceColor','r');
    goal_h = plot(0,0,'kp','MarkerSize',18,'MarkerFaceColor','k');
    hist_h = plot(0,0,'o','MarkerSize',10,'MarkerEdgeColor',[0.55 0.55 0.55],'LineWidth',1.5);

    % top-left legend (vehicle / obstacle / planned path)
    h_leg = [robot_h, traj_h, path_h, ...
             plot(nan, nan, 'o', 'MarkerSize', 9, ...         % static obstacle (gray circle)
                  'MarkerFaceColor',[0.72 0.72 0.72], 'MarkerEdgeColor','k'), ...
             plot(nan, nan, 'o', 'MarkerSize', 9, ...         % dynamic obstacle (orange circle)
                  'MarkerFaceColor',[1 0.65 0.25], 'MarkerEdgeColor','k'), ...
             goal_h];
    legend(h_leg, {'Vehicle', 'Actual trajectory', 'Planned path (RRT)', ...
                   'Static obstacle', 'Dynamic obstacle', 'Goal'}, ...
           'Location', 'northwest', 'FontSize', 8.5, 'Box', 'on');

    % Replay and Export-video buttons (side by side, top left)
    uicontrol('Style','pushbutton','String','Replay','FontSize',11, ...
              'Position',[10 690 100 32], 'Callback', @replay_cb);
    uicontrol('Style','pushbutton','String','Export video','FontSize',11, ...
              'Position',[120 690 100 32], 'Callback', @export_cb);

    n_steps = size(rec_traj,1);

    replay(false);   % first playback (no recording)

    % ---- nested function: replay once (records MP4 when do_record=true) ----
    function replay(do_record)
        if nargin < 1, do_record = false; end
        vid = [];
        if do_record
            vid = VideoWriter('mpc_cbf_animation', 'MPEG-4');
            vid.FrameRate = 20;
            open(vid);
        end
        set(traj_h,'XData',[],'YData',[]);       % clear trajectory
        set(hist_h,'XData',[],'YData',[]);       % clear past goals
        prev_g = 0;
        for k = 1:n_steps
            g = rec_g(k);
            % update path / goal / history when the goal switches
            if g ~= prev_g
                set(path_h, 'XData', rec_paths{g}(:,1), 'YData', rec_paths{g}(:,2));
                set(goal_h, 'XData', rec_goals(g,1), 'YData', rec_goals(g,2));
                if g > 1
                    set(hist_h, 'XData', rec_goals(1:g-1,1), 'YData', rec_goals(1:g-1,2));
                end
                prev_g = g;
            end
            % dynamic obstacles
            for i = 1:n_dyn
                pj = rec_dyn(k, 2*i-1:2*i);
                dyn_h(i).Position = [pj(1)-dyn_r, pj(2)-dyn_r, 2*dyn_r, 2*dyn_r];
            end
            % trajectory (up to the current step) + robot
            set(traj_h, 'XData', rec_traj(1:k,1), 'YData', rec_traj(1:k,2));
            set(robot_h, 'XData', rec_traj(k,1), 'YData', rec_traj(k,2));
            title(sprintf('goal %d/%d    t = %.1f s', g, N_GOALS, rec_t(k)), 'FontSize', 13, 'Interpreter','none');
            drawnow;
            if do_record
                writeVideo(vid, getframe(ax));
            end
            pause(0.02);   % ~5x playback speed; increase to slow down, decrease to speed up
        end
        title(sprintf('Playback finished @ t = %.1f s', rec_t(end)), 'FontSize', 13, 'Interpreter','none');
        if do_record
            close(vid);
            fprintf('Video saved: %s\n', fullfile(pwd, 'mpc_cbf_animation.mp4'));
        end
    end

    % ---- nested function: replay (no recording) ----
    function replay_cb(~, ~)
        replay(false);
    end

    % ---- nested function: export video (replay once and record MP4) ----
    function export_cb(~, ~)
        replay(true);
    end
end

function pos = dyn_pos(c, n, A, v, phi, dt)
    % ball crossing back and forth along perpendicular direction n
    % (dt = time relative to when the path was built); triangle-wave constant-speed motion
    pos = c + n * tri_wave(v*dt + phi, A);
end

function y = tri_wave(x, A)
    % triangle wave: maps x into [-A, A] with constant-speed back-and-forth motion (period 4A)
    x = mod(x + A, 4*A);
    y = A - abs(x - 2*A);
end

function g = sample_goal(c, static_obs, ROBOT_R, SAFE_D, GOAL_R)
    % sample a random goal in [-3.8,3.8]: far enough from the current position,
    % clear of static obstacles, and with 0.6m clearance from obstacle surfaces
    % (avoids edge crowding and getting stuck in static-obstacle pockets)
    for try_k = 1:500
        x = rand*7.6 - 3.8; y = rand*7.6 - 3.8;
        if norm([x y]-c(1:2)) < 2.5, continue; end
        ok = true;
        for j = 1:size(static_obs,1)
            if norm([x y]-static_obs(j,1:2)) < static_obs(j,3) + ROBOT_R + SAFE_D + 0.6
                ok = false; break;
            end
        end
        if ~ok, continue; end
        g = [x y]; return;
    end
    g = [4 4];
end

function p = rrt(s, g, obs, ROBOT_R, SAFE_D, GOAL_R)
    nodes = s; parent = 0;
    max_iter = 4000; step_len = 0.5;
    for it = 1:max_iter
        if rand < 0.1
            samp = g;
        else
            samp = [rand*10-5, rand*10-5];
        end
        [~, idx] = min(vecnorm(nodes - samp, 2, 2));
        near = nodes(idx, :);
        vec = samp - near;
        d = norm(vec);
        if d < 1e-6, continue; end
        new = near + vec/d * min(step_len, d);
        if ~seg_free(near, new, obs, ROBOT_R, SAFE_D), continue; end
        nodes = [nodes; new]; %#ok<AGROW>
        parent = [parent; idx]; %#ok<AGROW>
        if norm(new - g) < GOAL_R
            p = [g];
            cur = size(nodes,1);
            while cur > 0
                p = [nodes(cur,:); p]; %#ok<AGROW>
                cur = parent(cur);
            end
            return;
        end
    end
    p = [];
end

function ok = seg_free(a, b, obs, ROBOT_R, SAFE_D)
    d = norm(b-a); nseg = max(1, round(d/0.1));
    for ii = 0:nseg
        q = a + (b-a)*ii/nseg;
        for j = 1:size(obs,1)
            if norm(q - obs(j,1:2)) < obs(j,3) + ROBOT_R + SAFE_D
                ok = false; return;
            end
        end
    end
    ok = true;
end

function p = shortcut(p, obs, ROBOT_R, SAFE_D)
    out = p(1,:);
    ic = 1;
    np = size(p,1);
    while ic < np
        next_i = ic + 1;
        for j = np:-1:ic+1
            if seg_free(p(ic,:), p(j,:), obs, ROBOT_R, SAFE_D)
                next_i = j;
                break;
            end
        end
        out = [out; p(next_i,:)]; %#ok<AGROW>
        ic = next_i;
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
        tt = dot(c(1:2)-a, ab)/L2; tt = min(max(tt,0),1);
        proj = a + tt*ab;
        dd = norm(c(1:2)-proj);
        if dd < best_d, best_d = dd; best_i = i; best_proj = proj; end
    end
    pts = [best_proj; full(best_i+1:end,:)];
    cum = zeros(size(pts,1),1);
    for i = 2:size(pts,1)
        cum(i) = cum(i-1) + norm(pts(i,:)-pts(i-1,:));
    end
    ref = zeros(N+1, 3);
    ref(1,:) = c;
    for k = 1:N
        target = V_REF*DT*k;
        i = find(cum >= target, 1) - 1;
        if isempty(i) || i < 1, i = 1; end
        if i > size(pts,1)-1, i = size(pts,1)-1; end
        seg = pts(i+1,:) - pts(i,:);
        L = norm(seg);
        if L < 1e-9
            pt = pts(i,:);
        else
            frac = min(max((target-cum(i))/(cum(i+1)-cum(i)+1e-9),0),1);
            pt = pts(i,:) + seg*frac;
        end
        ref(k+1,1:2) = pt;
        ref(k+1,3) = atan2(ref(k+1,2)-ref(k,2), ref(k+1,1)-ref(k,1));
    end
    % reference-heading rate limiting: forward-difference heading jumps sharply
    % at corners, making MPC stop to turn; after rate limiting, heading changes
    % gradually and the robot turns while moving, removing corner stalls
    for k = 2:N+1
        d = atan2(sin(ref(k,3)-ref(k-1,3)), cos(ref(k,3)-ref(k-1,3)));
        max_d = 0.1;    % at most 0.1 rad per step (= OMEGA_MAX*DT, matches the robot's turning capability)
        if abs(d) > max_d
            ref(k,3) = ref(k-1,3) + sign(d)*max_d;
        end
    end
end

function obs_hor = obstacle_horizon(sobs, dyn_c, dyn_n, dyn_A, dyn_v, dyn_phi, tt, t_path, N, DT)
    n_s = size(sobs,1); n_d = size(dyn_c,1); n_o = n_s + n_d;
    obs_hor = zeros(n_o, N, 2);
    for j = 1:n_s
        for k = 1:N
            obs_hor(j,k,:) = sobs(j,1:2);
        end
    end
    for j = 1:n_d
        for k = 1:N
            obs_hor(n_s+j,k,:) = dyn_pos(dyn_c(j,:), dyn_n(j,:), dyn_A, dyn_v, dyn_phi(j), (tt+(k-1)*DT) - t_path);
        end
    end
end

function [u, z, ef] = solve_mpc(c, ref, obs_hor, last_z, N, DT, V_MAX, V_MIN, OMEGA_MAX, V_REF, GAMMA, r_eff_all, Q_XY, Q_TH, R_V, R_W, P_XY) %#ok<INUSD>
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
    if ef == -2
        z = [zeros(N,1); zeros(N,1)];
    end
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
