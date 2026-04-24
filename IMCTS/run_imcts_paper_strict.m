function out = run_imcts_paper_strict(varargin)
% Strict paper-style I-MCTS runner (Ye et al., 2025) for allocation.
% This runner enforces:
% 1) Two servicing spacecraft (Ns = 2)
% 2) Fuel-only objective/reward
% 3) Complete assignment required in final result
% 4) Paper dataset by default (Tables 7/8 style case data)
%
% Usage:
%   out = run_imcts_paper_strict();
%   out = run_imcts_paper_strict(40, 10000, 2025, true);
%
% Inputs:
%   nTargets : number of targets from paper dataset (default 40)
%   m        : MCTS iterations per decision step (default 10000)
%   seed     : RNG seed (default 2025)
%   verbose  : print progress/results (default true)

nTargets = 40;
m = 10000;
seed = 2025;
verbose = true;

if nargin >= 1, nTargets = varargin{1}; end
if nargin >= 2, m = varargin{2}; end
if nargin >= 3, seed = varargin{3}; end
if nargin >= 4, verbose = varargin{4}; end

[data, ref] = load_case_data('paper');
maxTargets = size(data.targets, 1);
assert(nTargets >= 3 && nTargets <= maxTargets, ...
    'nTargets must be in [3, %d] for paper data.', maxTargets);
assert(size(data.ssc, 1) == 2, 'Strict paper runner expects exactly 2 servicing spacecraft.');

caseData = data;
caseData.targets = data.targets(1:nTargets, :);
caseData.targetIds = data.targetIds(1:nTargets);

% Transfer model settings used in the IMCTS codebase.
model.mu = 398600.4418;   % km^3/s^2
model.Re = 6378.137;      % km
model.J2 = 1.08262668e-3;
model.dvPhaseWeight = 0.05;
model.timePhaseWeight = 0.20;
model.transferEta = 1.20;
model.kRangeEq4 = -20:20;
model.maxDriftTimeSec = inf;

% Calibration to match paper reference convention in this repo.
calib = struct('fuelScale', 1.0, 'timeScale', 1.0);
if nTargets >= 40
    ref40 = ref.hega40;
    cd40 = caseData;
    cd40.targets = data.targets(1:40, :);
    cd40.targetIds = data.targetIds(1:40);
    [f0, t0] = evaluate_solution(ref40.ssc1, ref40.ssc2, cd40, model, calib);
    calib.fuelScale = ref40.fuel / f0;
    calib.timeScale = ref40.time / t0;
end

% Strict paper-style optimization knobs.
params = struct();
params.m = m;
params.nchild = 10;              % Table 9
params.dvmax = 1.5;              % km/s
params.cUct = 0.8;
params.rngSeed = seed;
params.verbose = verbose;
params.infeasiblePenalty = 5.0;
params.objective = 'fuel';       % strict fuel objective
params.timeWeight = 0.0;         % unused in fuel mode
params.rolloutRoulette = true;   % stochastic rollout in simulation
params.strictCompleteFinal = true;

t0 = tic;
best = imcts_optimize(caseData, model, params, calib);
computeTimeSec = toc(t0);

if ~best.complete
    error('Strict paper run returned incomplete assignment (%d/%d).', ...
        best.assignedCount, best.totalTargets);
end

if verbose
    fprintf('\n=== Strict Paper I-MCTS (%d targets) ===\n', nTargets);
    fprintf('Objective: fuel only | Ns: 2 | complete assignment: required\n');
    fprintf('m = %d | nchild = %d | cUCT = %.2f | dvmax = %.2f km/s\n', ...
        params.m, params.nchild, params.cUct, params.dvmax);
    fprintf('Best fuel (km/s): %.4f\n', best.fuel);
    fprintf('Mission transfer time (scaled s): %.3f\n', best.time_sec_scaled);
    fprintf('Mission transfer time (Eq.4 sum, d): %.3f\n', best.tof_days_eq4);
    fprintf('SSc1 sequence (%d): [%s]\n', numel(best.seq{1}), fmt_seq(best.seq{1}));
    fprintf('SSc2 sequence (%d): [%s]\n', numel(best.seq{2}), fmt_seq(best.seq{2}));
    fprintf('Assigned total targets: %d/%d\n', best.assignedCount, best.totalTargets);
    fprintf('Compute time (s): %.3f\n', computeTimeSec);

    if nTargets == 40
        fprintf('\nPaper Table 11 I-MCTS reference: fuel=1.2352 km/s, time=126.535 s\n');
    end
    fprintf('\n');
end

out.best = best;
out.params = params;
out.model = model;
out.calibration = calib;
out.caseData = caseData;
out.compute_time_s = computeTimeSec;
end

function s = fmt_seq(v)
if isempty(v)
    s = '';
    return;
end
s = strjoin(string(v), ' ');
end
