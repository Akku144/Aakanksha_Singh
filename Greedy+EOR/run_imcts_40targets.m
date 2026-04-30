% run_imcts_40targets.m
% IMCTS + EOR reproduction scaffold for Ye et al. 40-target / 2-SSc case.
% Units: km, s, rad. Delta-v output: km/s.
clear; clc; rng(1);

P = constants_eor();
[targets, ssc] = data_40targets();        % oe = [a e i RAAN omega M]
cfg.Ns = 2;
cfg.Nt = 40;
cfg.m = 5;                            % Table 9
cfg.nchild = 5;                          % Table 9
cfg.Dvmax = 0.500;                        % 500 m/s = 0.5 km/s
cfg.cp = sqrt(2);
cfg.NmaxPhase = 5;                        % paper does not state Nmax for IMCTS case
cfg.verbose = true;

[result, root] = imcts_allocate(ssc, targets(1:40,:), cfg, P);

disp('Best allocation found:');
for k = 1:cfg.Ns
    fprintf('SSc%d: ', k); disp(result.seq{k});
end
fprintf('Total fuel = %.6f km/s\n', result.totalDV);

% Paper Table 11 reported:
% I-MCTS SSc1 [30 6 16 28 33 24 20 19 2 11 26 29 18 1 12 14 39 25 17 31 21 7 22 35 23 37 13 40 34 9 5 32]
% I-MCTS SSc2 [36 3 4 27 15 10 8 38]
% Fuel = 1.2352 km/s, Time = 126.535 s
