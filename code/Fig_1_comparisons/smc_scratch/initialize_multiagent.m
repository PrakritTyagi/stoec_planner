function [poses, opt] = initialize_multiagent()
    opt = struct();

    DomainBounds.xmin = 0.0;
    DomainBounds.xmax = 150.0;
    DomainBounds.ymin = 0.0;
    DomainBounds.ymax = 150.0;
    opt.DomainBounds = DomainBounds;
    opt.L = [DomainBounds.xmax - DomainBounds.xmin; DomainBounds.ymax - DomainBounds.ymin];

    N_teamB = 1; % Example
    N_teamA = N_teamB + 1;
    opt.nagents = N_teamA + N_teamB;
    opt.teamA.idx = 1:N_teamA;
    opt.teamB.idx = N_teamA+1:opt.nagents;

    opt.vlb = 0.1 * ones(opt.nagents,1);
    opt.vub = 5.0 * ones(opt.nagents,1);
    opt.wlb = -0.2 * ones(opt.nagents,1);
    opt.wub = 0.2 * ones(opt.nagents,1);

    % poses.x = opt.DomainBounds.xmax * rand(opt.nagents,1);
    % poses.y = opt.DomainBounds.ymax * rand(opt.nagents,1);
    % poses.theta = 2*pi*rand(opt.nagents,1);
    poses.x = [140,100,25];
    poses.y = [20,40,125];
    poses.theta = 120*pi/180*ones(opt.nagents,1);

    opt.sim.Nsteps =1000;
    opt.sim.dt = 0.1;

    opt.erg.s = 1.5;
    opt.erg.Nkx = 50;
    opt.erg.Nky = 50;
    opt.erg.KX = (0:opt.erg.Nkx-1)' * ones(1,opt.erg.Nky);
    opt.erg.KY = ones(opt.erg.Nkx,1) * (0:opt.erg.Nky-1);
    opt.erg.LK = 1.0 ./ ((1.0 + opt.erg.KX.^2 + opt.erg.KY.^2).^opt.erg.s);
end