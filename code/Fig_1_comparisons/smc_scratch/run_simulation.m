function [ErgA, ErgB] = run_simulation(poses, opt, teamA_strategy, teamB_strategy, X, Y)
    Nsteps = opt.sim.Nsteps;
    dt = opt.sim.dt;
    traj = zeros(Nsteps, opt.nagents, 3);
    erg_flags = false(opt.nagents,1);
    detections = cell(opt.nagents, 1);
    detection_map = zeros(size(X));
    % epsilon = 1e-10;
    % detection_map = detection_map + epsilon;
    bhatt_distance = zeros(Nsteps,1);
    ErgA = zeros(Nsteps, 1);
    ErgB = zeros(Nsteps, 1);

    % plot the team A infomap
    figure(1); set(gcf,'color','w'); hold on
    surface(X,Y,zeros(size(X)),reshape(opt.teamA.map,size(X)),'FaceColor','interp','EdgeColor','none');
    axis tight; axis equal; title('Multiagent Simulation (Team A map)');

    % plot start locations of every agent
    scatter(poses.x(opt.teamA.idx), poses.y(opt.teamA.idx), 120, 'go', 'filled', 'DisplayName', 'Team A');
    scatter(poses.x(opt.teamB.idx), poses.y(opt.teamB.idx), 120, 'rs', 'filled', 'DisplayName', 'Team B');

    erg_flags(opt.teamA.idx(1)) = true;

    fprintf('***Simulation Starts***\n')
    pause(2)
    for step = 1:Nsteps
        time = step * dt;
        fprintf('Time step (dt = %f) : %d \n',opt.sim.dt,step)
        for i = 1:opt.nagents
            % Check which agent is i?
            isTeamA = ismember(i, opt.teamA.idx);
            strategy = teamA_strategy * isTeamA + teamB_strategy * ~isTeamA;
            
            % Update Ck
            xrel = poses.x(i) - opt.DomainBounds.xmin;
            yrel = poses.y(i) - opt.DomainBounds.ymin;
            HK = isTeamA * opt.teamA.HK + ~isTeamA * opt.teamB.HK;
            Ck_update = cos(opt.erg.KX * pi * xrel / opt.L(1)) .* cos(opt.erg.KY * pi * yrel / opt.L(2)) * dt ./ HK';

            if isTeamA
                opt.Ck_teamA = opt.Ck_teamA + Ck_update;
            else
                opt.Ck_teamB = opt.Ck_teamB + Ck_update;
            end

            if isTeamA && strategy == 3
                if erg_flags(i)
                    for j = opt.teamB.idx
                        if is_agent_visible(poses.x(i), poses.y(i), poses.x(j), poses.y(j), opt.visibility_range)
                            % order of false and true matters matters.
                            new_erg = find_idle_ergodic(erg_flags(opt.teamA.idx));
                            disp(new_erg)
                            erg_flags(i) = false;
                            if ~isempty(new_erg)
                                erg_flags(opt.teamA.idx(new_erg)) = true;
                            end
                            break;
                        end
                    end
                end
                if erg_flags(i)
                    poses = apply_ergodic_step(poses, i, time, opt, true);
                else
                    poses = apply_random_step(poses, i, opt);
                end

            elseif ~isTeamA && (strategy == 1 || strategy == 2 || strategy == 3)
                visible_a = [];
                for j = opt.teamA.idx
                    if is_agent_visible(poses.x(i), poses.y(i), poses.x(j), poses.y(j), opt.visibility_range)
                        visible_a = [visible_a; poses.x(j), poses.y(j)];
                    end
                end
                detections{i} = [detections{i}; visible_a];

                if ~isempty(visible_a)
                    [xq, yq] = meshgrid(X(1,:), Y(:,1));
                    for pt = 1:size(visible_a, 1)
                        mu = visible_a(pt,:);
                        Sigma = eye(2) * 30; % adjust the spread as needed
                        gauss = mvnpdf([xq(:) yq(:)], mu, Sigma);
                        detection_map = detection_map + reshape(gauss, size(X));
                        if strategy == 3
                            opt.teamB.map = opt.teamB.map + reshape(gauss, size(opt.teamB.map));
                        end
                    end
                    % inds = knnsearch(opt.kdOBJ, visible_a);
                    % detection_map(inds) = detection_map(inds) + 1;
                    % 
                    % if strategy == 3
                    %     % Team B uses detected points to update its own info map
                    %     for idx = inds'
                    %         opt.teamB.map(idx) = opt.teamB.map(idx) + 1;
                    %     end
                    % end
                end

                if (strategy == 2 || strategy == 3) && ~isempty(visible_a)
                    tx = visible_a(1,1);
                    ty = visible_a(1,2);
                    [poses.x(i), poses.y(i), poses.theta(i)] = pursue_target(poses.x(i), poses.y(i), poses.theta(i), tx, ty, opt, i);
                else
                    poses = apply_ergodic_step(poses, i, time, opt, false);
                end

            else
                if strategy == 1
                    poses = apply_ergodic_step(poses, i, time, opt, isTeamA);
                else
                    if erg_flags(i)
                        poses = apply_ergodic_step(poses, i, time, opt, true);
                    else
                        poses = apply_random_step(poses, i, opt);
                    end
                end
            end

            traj(step,i,:) = [poses.x(i), poses.y(i), poses.theta(i)];
        end

        ckA = opt.Ck_teamA / (length(opt.teamA.idx) * time);
        ckB = opt.Ck_teamB / (length(opt.teamB.idx) * time);
        ErgA(step) = sum(sum(opt.erg.LK .* (ckA - opt.teamA.muk).^2));
        ErgB(step) = sum(sum(opt.erg.LK .* (ckB - opt.teamB.muk).^2));

        bhatt_distance(step) = evaluateBhattacharyyaDist(opt.teamA.map,detection_map);

        % Use if plotting trajectories every 5 time step
        % if mod(step, 5) == 0
        %     for i = 1:opt.nagents
        %         plot(traj(1:step,i,1), traj(1:step,i,2));
        %     end
        %     drawnow;
        % end
    end
    fprintf('***Simulation Ends***\n')

    % Use if plot trajectories only at the end
    for i = 1:opt.nagents
        plot(traj(1:step,i,1), traj(1:step,i,2));
    end
    drawnow;

    detection_map = detection_map / sum(detection_map(:));

    % Normalize updated teamB map if strategy 3 was used
    if teamB_strategy == 3
        opt.teamB.map = opt.teamB.map / sum(opt.teamB.map(:));
    end

    % Plot estimated map of A by B
    figure(2); set(gcf,'color','w');
    surface(X,Y,zeros(size(X)),detection_map,'FaceColor','interp','EdgeColor','none');
    axis tight; axis equal; title('Team A estimated map by B'); colorbar;

    % Plot Info map of B
    figure(3); set(gcf,'color','w');
    surface(X,Y,zeros(size(X)),reshape(opt.teamB.original_map,size(X)),'FaceColor','interp','EdgeColor','none');
    title('Original Team B Map'); axis equal; colorbar;

    assignin('base','teamB_detections',detections);
    assignin('base','teamB_detection_map',detection_map);
    assignin('base','bhatt_distance',bhatt_distance);

    % Print difference between estimated map of A and team A map.
    map_diff = detection_map(:) - opt.teamA.map;
    detection_error = sum((map_diff(:)).^2);
    fprintf('Detection error (Team A estimated map vs TeamA map): %e \n', detection_error);
    assignin('base','detection_error',detection_error);


    % Compare original and updated Team B map if strategy 3 used
    if teamB_strategy == 3
        figure(4); set(gcf,'color','w');
        subplot(1,2,1);
        surface(X,Y,zeros(size(X)),reshape(opt.teamB.original_map,size(X)),...
                'FaceColor','interp','EdgeColor','none');
        title('Original Team B Map'); axis equal; colorbar;

        subplot(1,2,2);
        surface(X,Y,zeros(size(X)),reshape(opt.teamB.map,size(X)),...
                'FaceColor','interp','EdgeColor','none');
        title('Updated Team B Map'); axis equal; colorbar;
    end
end
