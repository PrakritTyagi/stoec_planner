function [x_new, y_new, theta_new] = pursue_target(x, y, theta, tx, ty, opt, i)
    dt = opt.sim.dt;
    dx = tx - x;
    dy = ty - y;
    target_theta = atan2(dy, dx);
    w = wrapToPi(target_theta - theta);
    v = opt.vub(i);
    x_new = x + v * cos(theta) * dt;
    y_new = y + v * sin(theta) * dt;
    theta_new = theta + w * dt;
end