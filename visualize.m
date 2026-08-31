function visualize(play)
% VISUALIZE  Run one full simulation and export three key-frame figures (PNG)
%   showing the RRT + MPC + D-CBF architecture.
%   Usage: visualize           export the three PNGs
%          visualize(true)     export the three PNGs, then loop-play them
%                              in the same figure (1 s per frame, 3 loops)
%   Output (current folder), one run produces all three:
%     mpc_cbf_before.png  : moment BEFORE the robot changes course (trajectory still straight)
%     mpc_cbf_result.png  : moment the robot is closest to the dynamic ball (just passing it)
%     mpc_cbf_final.png   : goal reached, full trajectory
%   The dynamic ball's CBF safety ring follows the ball in every frame.
%   (Same architecture as simulate.m; the live animation is replaced by three
%    static key frames. The dynamic-ball parameters can be recomputed with
%    compute_balls.m)

    if nargin < 1, play = false; end

    % ============ Parameters (same as simulate.m) ============
    DT = 0.1; N = 12; GAMMA = 0.3;
    V_MAX = 1.0; V_MIN = 0.0; OMEGA_MAX = 1.0;
    ROBOT_R = 0.15; SAFE_D = 0.25; GOAL_R = 0.25; V_REF = 0.8;
    % SAFE_D = 0.25: larger safety margin than simulate.m's 0.10, making the
    %   D-CBF detours visually pronounced in the demo figures
    % V_REF = 0.8: simulate.m's default; fast enough to clear narrow passages
    %   between static obstacles
    Q_XY = 10; Q_TH = 0.1; R_V = 0.01; R_W = 0.02; P_XY = 50;

    % ============ Obstacles (default scenario) ============
    static_obs = [-2.0  1.5 0.6;
                   1.0  2.0 0.7;
                   2.5 -1.0 0.5;
                  -1.5 -2.5 0.6];
    dyn_obs = [ 2.64 -4.00  0.45 -0.26  0.43];   % single ball at 0.5 m/s, timed to cross the robot's path at (0.5,-0.5) at t=8.2s (parameters from compute_balls.m)
    start = [-4 -4 0];
    goal  = [4 4];
    r_eff_all = [static_obs(:,3)' + ROBOT_R + SAFE_D, dyn_obs(:,3)' + ROBOT_R + SAFE_D];

    % ============ RRT reference path ============
    path = rrt(start(1:2), goal, static_obs);
    if isempty(path), error('RRT found no path'); end
    path = shortcut(path, static_obs);

    % ============ Simulation loop (headless, collect trajectory) ============
    curr = start; traj = start(1:2); t = 0; last_z = [];
    times = 0;                        % time at each trajectory point
    dyn_traj = cell(1, size(dyn_obs,1));   % record dynamic-obstacle history positions
    for i = 1:size(dyn_obs,1), dyn_traj{i} = dyn_obs(i,1:2); end
    min_dist = inf; collision_printed = false;
    capture_i = 1; min_ball_dist = inf;   % key-frame: step where the robot is closest to the dynamic ball (just passing it)

    for step = 1:1500
        if norm(curr(1:2) - goal) < GOAL_R, break; end
        ref = make_ref(path, curr, goal);
        obs_hor = obstacle_horizon(static_obs, dyn_obs, t);
        [u, last_z] = solve_mpc(curr, ref, obs_hor, last_z); %#ok<ASGLU>
        curr = curr + [u(1)*cos(curr(3))*DT, u(1)*sin(curr(3))*DT, u(2)*DT];
        traj = [traj; curr(1:2)]; %#ok<AGROW>
        t = t + DT;
        times = [times; t]; %#ok<AGROW>
        for i = 1:size(dyn_obs,1)
            dyn_traj{i} = [dyn_traj{i}; dyn_obs(i,1) + dyn_obs(i,4)*t, dyn_obs(i,2) + dyn_obs(i,5)*t];
        end
        % key-frame tracking: robot closest to the dynamic ball (moment of passing it)
        db = norm(curr(1:2) - dyn_traj{1}(end,:));
        if db < min_ball_dist
            min_ball_dist = db; capture_i = size(traj,1);
        end
        for i = 1:size(static_obs,1)
            d = norm(curr(1:2)-static_obs(i,1:2)) - static_obs(i,3) - ROBOT_R;
            if d < 0 && ~collision_printed
                fprintf('  [collision] t=%.1fs static obstacle %d surface dist=%.3f\n', t, i, d); collision_printed = true;
            end
            min_dist = min(min_dist, d);
        end
        for i = 1:size(dyn_obs,1)
            xy = dyn_traj{i}(end,:);
            d = norm(curr(1:2)-xy) - dyn_obs(i,3) - ROBOT_R;
            if d < 0 && ~collision_printed
                fprintf('  [collision] t=%.1fs dynamic obstacle %d surface dist=%.3f\n', t, i, d); collision_printed = true;
            end
            min_dist = min(min_dist, d);
        end
    end

    % ============ Render three key frames ============
    %   1) mpc_cbf_before.png : a moment BEFORE the robot changes course (trajectory still straight)
    %   2) mpc_cbf_result.png : the moment the robot is closest to the dynamic ball (just passing it)
    %   3) mpc_cbf_final.png  : goal reached, full trajectory
    fig = figure('Color','w','Position',[100 100 860 820]);
    end_i = size(traj, 1);
    pre_i = max(1, capture_i - 25);   % ~2.5 s before the encounter
    if play
        % draw the three frames IN THIS figure (exporting PNGs), then loop-play
        % them in the same window: 1 s per frame, 3 loops total
        draw_frame(pre_i,     'mpc_cbf_before.png', 'before the detour', true);
        draw_frame(capture_i, 'mpc_cbf_result.png', 'just passing the dynamic obstacle', true);
        draw_frame(end_i,     'mpc_cbf_final.png',  'goal reached', true);
        for rep = 2:3
            draw_frame(pre_i,     '', 'before the detour', false);
            draw_frame(capture_i, '', 'just passing the dynamic obstacle', false);
            draw_frame(end_i,     '', 'goal reached', false);
        end
    else
        draw_frame(pre_i,     'mpc_cbf_before.png', 'before the detour', true);
        draw_frame(capture_i, 'mpc_cbf_result.png', 'just passing the dynamic obstacle', true);
        draw_frame(end_i,     'mpc_cbf_final.png',  'goal reached', true);
    end
    fprintf('Reached goal @ t=%.1fs, %d steps, min surface dist=%.3fm (%s)\n', t, step, min_dist, ...
            ternary(min_dist<0,'collision!','safe'));

    % ---- nested function: render one frame at trajectory index idx ----
    %   save_png=true  -> also export the frame to fname
    %   save_png=false -> draw in the figure and pause 1 s (loop playback)
    function draw_frame(idx, fname, tag, save_png)
        if nargin < 4, save_png = true; end
        clf(fig); hold on; axis equal; grid on; box on;
        xlim([-5.2 5.2]); ylim([-5.2 5.2]);
        % safety inflation zones: static obstacles fixed, dynamic ball at its CURRENT frame position (follows the ball)
        for i = 1:size(static_obs,1)
            draw_circle(static_obs(i,1), static_obs(i,2), static_obs(i,3)+ROBOT_R+SAFE_D, [1 0.82 0.82], 0.5);
        end
        for i = 1:size(dyn_obs,1)
            bp = dyn_traj{i}(idx,:);
            draw_circle(bp(1), bp(2), dyn_obs(i,3)+ROBOT_R+SAFE_D, [1 0.82 0.82], 0.5);
        end
        % static obstacles
        for i = 1:size(static_obs,1)
            draw_circle(static_obs(i,1), static_obs(i,2), static_obs(i,3), [0.62 0.62 0.62], 1);
        end
        % dynamic obstacle: history up to this frame + position at this frame
        for i = 1:size(dyn_obs,1)
            dtr = dyn_traj{i}(1:idx, :);
            plot(dtr(:,1), dtr(:,2), 'Color',[1 0.6 0.2 0.45], 'LineWidth',1, 'LineStyle','-');
            draw_circle(dtr(end,1), dtr(end,2), dyn_obs(i,3), [1 0.55 0.15], 1);
        end
        % RRT reference path (green dashed) + actual trajectory up to this frame (red solid)
        plot(path(:,1), path(:,2), 'g--', 'LineWidth', 1.8);
        plot(traj(1:idx,1), traj(1:idx,2), 'r-', 'LineWidth', 2.2);
        % start / goal / robot at this frame
        plot(start(1), start(2), 's', 'MarkerSize', 14, 'MarkerFaceColor',[0.2 0.8 0.2], 'MarkerEdgeColor','k');
        plot(goal(1), goal(2), 'p', 'MarkerSize', 22, 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
        plot(traj(idx,1), traj(idx,2), 'ro', 'MarkerSize', 12, 'MarkerFaceColor','r');
        draw_circle(traj(idx,1), traj(idx,2), ROBOT_R, 'none', 1);  % robot outline
        add_legend(gca);
        xlabel('x (m)'); ylabel('y (m)');
        title(sprintf('RRT + MPC + D-CBF   %s (t=%.1fs)   min surface dist=%.3fm (%s)', tag, times(idx), min_dist, ...
              ternary(min_dist<0, 'collision!', 'safe')), ...
              'FontSize',12, 'Interpreter','none');
        if save_png
            exportgraphics(fig, fname, 'Resolution', 150);
            fprintf('Saved %s (frame t=%.1fs)\n', fname, times(idx));
        end
        drawnow;
        pause(1.0);   % 1 s per frame: continuous playback whether or not we are exporting PNGs
    end

    % ==================== Helpers ====================
    function draw_circle(x, y, r, fc, alpha)
        if strcmp(fc,'none')
            rectangle('Position',[x-r y-r 2*r 2*r],'Curvature',[1 1],'EdgeColor','r','LineWidth',1.5);
        else
            rectangle('Position',[x-r y-r 2*r 2*r],'Curvature',[1 1],'FaceColor',fc,'EdgeColor','k',...
                      'FaceAlpha',alpha,'LineWidth',0.8);
        end
    end

    % ============ Top-left legend ============
    function add_legend(ax)
        % Place a legend box in the top-left corner explaining each graphic element:
        %   Vehicle (red dot) / Actual trajectory (red line) / Planned path (green
        %   dashed) / Safe distance (light-red CBF ring) / Static obstacle (gray
        %   circle) / Dynamic obstacle (orange circle) / Goal (black star)
        % Hidden handles created with plot(nan,nan,...) serve as legend entries
        % without interfering with the actual plot
        h = [plot(ax, nan, nan, 'o', 'MarkerSize', 9, ...          % Vehicle (red dot)
                  'MarkerFaceColor','r', 'MarkerEdgeColor','r');
             plot(ax, nan, nan, 'r-', 'LineWidth', 2.2);          % Actual trajectory (red solid)
             plot(ax, nan, nan, 'g--', 'LineWidth', 1.8);         % Planned path (green dashed)
             plot(ax, nan, nan, 'o', 'MarkerSize', 12, ...        % Safe distance (light-red CBF ring)
                  'MarkerFaceColor',[1 0.82 0.82], 'MarkerEdgeColor','k');
             plot(ax, nan, nan, 'o', 'MarkerSize', 9, ...         % static obstacle (gray circle)
                  'MarkerFaceColor',[0.62 0.62 0.62], 'MarkerEdgeColor','k');
             plot(ax, nan, nan, 'o', 'MarkerSize', 9, ...         % dynamic obstacle (orange circle)
                  'MarkerFaceColor',[1 0.55 0.15], 'MarkerEdgeColor','k');
             plot(ax, nan, nan, 'p', 'MarkerSize', 11, ...        % goal (black star)
                  'MarkerFaceColor','k', 'MarkerEdgeColor','k')];
        legend(h, {'Vehicle', 'Actual trajectory', 'Planned path (RRT)', ...
                   'Safe distance (CBF)', 'Static obstacle', 'Dynamic obstacle', 'Goal'}, ...
               'Location', 'northwest', 'FontSize', 8.5, 'Box', 'on');
    end

    function s = ternary(cond, a, b)
        if cond, s = a; else, s = b; end
    end

    % ==================== RRT ====================
    function p = rrt(s, g, obs)
        rng(0);   % fixed seed: demo figures are reproducible run-to-run
        nodes = s; parent = 0;
        for it = 1:4000
            if rand < 0.1, samp = g; else, samp = [rand*10-5, rand*10-5]; end
            [~, idx] = min(vecnorm(nodes - samp, 2, 2));
            near = nodes(idx,:); vec = samp - near; d = norm(vec);
            if d < 1e-6, continue; end
            new = near + vec/d * min(0.5, d);
            if ~seg_free(near, new, obs), continue; end
            nodes = [nodes; new]; %#ok<AGROW>
            parent = [parent; idx]; %#ok<AGROW>
            if norm(new - g) < GOAL_R
                p = [g]; cur = size(nodes,1);
                while cur > 0, p = [nodes(cur,:); p]; cur = parent(cur); end %#ok<AGROW>
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
                if norm(p - obs(j,1:2)) < obs(j,3) + ROBOT_R + SAFE_D
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
            out = [out; p(next_i,:)]; %#ok<AGROW>
            i = next_i;
        end
        p = out;
    end

    function ref = make_ref(p, c, g)
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

    function obs_hor = obstacle_horizon(sobs, dobs, tt)
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

    function [u, z] = solve_mpc(c, ref, obs_hor, last_z) %#ok<INUSD>
        dth = atan2(sin(ref(2,3)-c(3)), cos(ref(2,3)-c(3)));
        w0 = min(max(dth/DT, -OMEGA_MAX), OMEGA_MAX);
        lb = [repmat(V_MIN,N,1); repmat(-OMEGA_MAX,N,1)];
        ub = [repmat(V_MAX,N,1); repmat( OMEGA_MAX,N,1)];
        z0 = [repmat(V_REF,N,1); repmat(w0,N,1)];
        opts = optimoptions('fmincon', 'Display','off', 'Algorithm','sqp', ...
                            'MaxIterations', 150, 'OptimalityTolerance', 1e-3, ...
                            'ConstraintTolerance', 1e-3, 'StepTolerance', 1e-6);
        [z, ~, ef] = fmincon(@(z) mpc_cost(z, c, ref), z0, [], [], [], [], lb, ub, ...
                             @(z) mpc_nonlcon(z, c, obs_hor), opts);
        if ef == -2, z = [zeros(N,1); zeros(N,1)]; end
        u = [z(1); z(N+1)];
    end

    function J = mpc_cost(z, c, ref)
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

    function [cc, ceq] = mpc_nonlcon(z, c, obs_hor)
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
end
