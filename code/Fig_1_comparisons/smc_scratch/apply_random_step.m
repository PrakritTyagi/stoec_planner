function poses = apply_random_step(poses, i, opt)
    dt = opt.sim.dt;
    v = rand * (opt.vub(i) - opt.vlb(i)) + opt.vlb(i);
    w = rand * (opt.wub(i) - opt.wlb(i)) + opt.wlb(i);
    poses.theta(i) = poses.theta(i) + w * dt;
    poses.x(i) = poses.x(i) + v * cos(poses.theta(i)) * dt;
    poses.y(i) = poses.y(i) + v * sin(poses.theta(i)) * dt;
end