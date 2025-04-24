clear; clc; close all;

%% Initialization
[init_poses, opt] = initialize_multiagent(); % change parameters inside this function

% Generate different utility maps for both teams
[X, Y, infoMapA] = GenerateUtilityMapA(opt, 0);
[~, ~, infoMapB] = GenerateUtilityMapB(opt, 0); % You can change parameters for different maps

% Assign utility maps to each team
opt.teamA.map = infoMapA / sum(infoMapA(:));
opt.teamB.map = infoMapB / sum(infoMapB(:));
opt.teamB.original_map = opt.teamB.map;

% KDTree for spatial queries
opt.kdOBJ = KDTreeSearcher([X(:),Y(:)]);

% Fourier coefficients for each team
[opt.teamA.muk, opt.teamA.HK] = GetFourierCoeff(opt, X, Y, opt.teamA.map);
[opt.teamB.muk, opt.teamB.HK] = GetFourierCoeff(opt, X, Y, opt.teamB.map);

% Assign team strategies
teamA_strategy = 3; % 1 / 2 / 3
teamB_strategy = 3; % 1 / 2 / 3

%% Initialize cooperative Fourier coefficients for teams
opt.Ck_teamA = zeros(opt.erg.Nkx, opt.erg.Nky);
opt.Ck_teamB = zeros(opt.erg.Nkx, opt.erg.Nky);

%% Run simulation
[ErgA, ErgB] = run_simulation(init_poses, opt, teamA_strategy, teamB_strategy, X, Y);

fprintf('Erodicity of team A : %e \n',ErgA(end))
%% Plot ergodicity metrics over time
timeVec = (1:length(ErgA)) * opt.sim.dt;
figure; plot(timeVec, ErgA, 'b', 'LineWidth', 1.5); hold on;
plot(timeVec, ErgB, 'r', 'LineWidth', 1.5);
legend('Team A', 'Team B');
xlabel('Time (s)'); ylabel('Ergodicity Metric');
title('Ergodicity Metric Over Time'); grid on;