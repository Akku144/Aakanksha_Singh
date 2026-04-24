function bestAll = run_imcts_bestof(nTargets, nRuns, m, seed0)
% Run multiple I-MCTS seeds and keep best-fuel solution.
if nargin < 1, nTargets = 40; end
if nargin < 2, nRuns = 20; end
if nargin < 3, m = 20000; end
if nargin < 4, seed0 = 2025; end

bestFuel = inf;
bestOut = [];
fuels = zeros(nRuns,1);
ct = zeros(nRuns,1);

for k = 1:nRuns
    seed = seed0 + k - 1;
    out = run_imcts_repro(nTargets, m, seed, false);
    fuels(k) = out.best.fuel;
    ct(k) = out.compute_time_s;
    if out.best.fuel < bestFuel
        bestFuel = out.best.fuel;
        bestOut = out;
    end
    fprintf('run %d/%d seed=%d fuel=%.4f compute=%.2fs\n', k, nRuns, seed, out.best.fuel, out.compute_time_s);
end

fprintf('\n=== I-MCTS Best-Of Summary (%d targets) ===\n', nTargets);
fprintf('Best fuel: %.4f km/s\n', bestOut.best.fuel);
fprintf('Best compute time: %.3f s\n', bestOut.compute_time_s);
fprintf('Fuel min/mean/max over runs: %.4f / %.4f / %.4f\n', min(fuels), mean(fuels), max(fuels));
fprintf('Compute min/mean/max over runs: %.2f / %.2f / %.2f s\n', min(ct), mean(ct), max(ct));

if nTargets == 40
    fprintf('Paper I-MCTS (Table 11): fuel=1.2352, compute time=126.535 s\n');
end

bestAll = bestOut;
end
