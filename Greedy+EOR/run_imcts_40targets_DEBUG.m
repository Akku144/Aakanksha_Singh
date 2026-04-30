% run_imcts_40targets_DEBUG.m
% Fast debug runner. Run this first.

clear; clc; rng(1);

P = constants_eor();
[targets, ssc] = data_40targets();

cfg.Ns = 2;
cfg.Nt = 40;              % first test only 5 targets
cfg.m = 500;               % not used in simplified debug allocator
cfg.nchild = 3;
cfg.Dvmax = 5.000;
cfg.cp = sqrt(2);
cfg.NmaxPhase = 1;
cfg.verbose = true;

[result, root] = imcts_allocate(ssc, targets(1:cfg.Nt,:), cfg, P);

disp('Best allocation found:');
for k = 1:cfg.Ns
    fprintf('SSc%d: ', k);
    disp(result.seq{k});
end
fprintf('Total fuel = %.6f km/s\n', result.totalDV);
