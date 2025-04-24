function poses = apply_ergodic_step(poses, i, time, opt, isTeamA)
    KX = opt.erg.KX;
    KY = opt.erg.KY;
    LK = opt.erg.LK;
    Lx = opt.L(1);
    Ly = opt.L(2);
    xmin = opt.DomainBounds.xmin;
    ymin = opt.DomainBounds.ymin;
    dt = opt.sim.dt;

    if isTeamA
        HK = opt.teamA.HK;
        muk = length(opt.teamA.idx) * time * opt.teamA.muk;
        % Ck = opt.Ck_teamA / (length(opt.teamA.idx) * time);
        Ck = opt.Ck_teamA;
    else
        HK = opt.teamB.HK;
        muk = length(opt.teamB.idx) * time * opt.teamB.muk;
        % Ck = opt.Ck_teamB / (length(opt.teamB.idx) * time);
        Ck = opt.Ck_teamB;
    end

    xrel = poses.x(i) - xmin;
    yrel = poses.y(i) - ymin;

    Bjx = sum(sum(LK./HK' .* (Ck - muk) .* (-KX * pi / Lx .* sin(KX * pi * xrel/Lx) .* cos(KY * pi * yrel/Ly))));
    Bjy = sum(sum(LK./HK' .* (Ck - muk) .* (-KY * pi / Ly .* cos(KX * pi * xrel/Lx) .* sin(KY * pi * yrel/Ly))));

    GammaV = Bjx * cos(poses.theta(i)) + Bjy * sin(poses.theta(i));
    GammaW = -Bjx * sin(poses.theta(i)) + Bjy * cos(poses.theta(i));

    if GammaV >= 0
        v = opt.vlb(i);
    else
        v = opt.vub(i);
    end
    if GammaW >= 0
        w = opt.wlb(i);
    else
        w = opt.wub(i);
    end



     %velocity motion model
    if(abs(w) < 1e-10 )
        poses.x(i) = poses.x(i) + v*dt*cos(poses.theta(i));   
        poses.y(i) = poses.y(i) + v*dt*sin(poses.theta(i));
    else
        poses.x(i) = poses.x(i) + v/w*(sin(poses.theta(i) + w*dt) - sin(poses.theta(i)));   
        poses.y(i) = poses.y(i) + v/w*(cos(poses.theta(i)) - cos(poses.theta(i)+ w*dt));    
    end
    poses.theta(i) = poses.theta(i) + w * dt;
    % poses.x(i) = poses.x(i) + v * cos(poses.theta(i)) * dt;
    % poses.y(i) = poses.y(i) + v * sin(poses.theta(i)) * dt;
end