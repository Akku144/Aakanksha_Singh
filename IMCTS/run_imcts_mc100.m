function stats = run_imcts_mc100(nTargets, nRuns, m)
% Monte-Carlo robustness for I-MCTS.
if nargin < 1, nTargets = 50; end
if nargin < 2, nRuns = 100; end
if nargin < 3, m = 10000; end

fuels = zeros(nRuns, 1);
times = zeros(nRuns, 1);

for k = 1:nRuns
    seed = 2025 + k - 1;
    out = run_imcts_repro(nTargets, m, seed, false);
    fuels(k) = out.best.fuel;
    times(k) = out.best.time;
end

stats.nTargets = nTargets;
stats.nRuns = nRuns;
stats.m = m;
stats.fuel_min = min(fuels);
stats.fuel_max = max(fuels);
stats.fuel_mean = mean(fuels);
stats.fuel_std = std(fuels);
stats.time_min = min(times);
stats.time_max = max(times);
stats.time_mean = mean(times);
stats.time_std = std(times);

fprintf('\n=== I-MCTS Monte-Carlo (%d targets, %d runs) ===\n', nTargets, nRuns);
fprintf('Fuel  min/max/mean/std: %.4f / %.4f / %.4f / %.4f km/s\n', ...
    stats.fuel_min, stats.fuel_max, stats.fuel_mean, stats.fuel_std);
fprintf('Time  min/max/mean/std: %.3f / %.3f / %.3f / %.3f s\n', ...
    stats.time_min, stats.time_max, stats.time_mean, stats.time_std);

if nTargets == 40
    fprintf('\nPaper Table 12 (I-MCTS, 40 targets):\n');
    fprintf('Fuel min/max: 1.1970 / 1.2395 km/s\n');
    fprintf('Time min/max: 110.71 / 180.21 s\n');
end
if nTargets == 50
    fprintf('\nPaper Table 13 average time (I-MCTS, 50 targets): 230.3963 s\n');
end
end
