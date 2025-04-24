function poses = apply_random_step(poses, i, opt)
    dt = opt.sim.dt;
    v = rand * (opt.vub(i) - opt.vlb(i)) + opt.vlb(i);
    w = rand * (opt.wub(i) - opt.wlb(i)) + opt.wlb(i);
    poses.theta(i) = poses.theta(i) + w * dt;

    % Compute new position
    new_x = poses.x(i) + v * cos(poses.theta(i)) * dt;
    new_y = poses.y(i) + v * sin(poses.theta(i)) * dt;

    % Reflect the agent if it hits domain boundary
    if new_x < opt.DomainBounds.xmin || new_x > opt.DomainBounds.xmax
        poses.theta(i) = pi - poses.theta(i); % reflect horizontally
        new_x = poses.x(i); % don't move out
    end
    if new_y < opt.DomainBounds.ymin || new_y > opt.DomainBounds.ymax
        poses.theta(i) = -poses.theta(i); % reflect vertically
        new_y = poses.y(i); % don't move out
    end

    % Update position
    poses.x(i) = new_x;
    poses.y(i) = new_y;
end
