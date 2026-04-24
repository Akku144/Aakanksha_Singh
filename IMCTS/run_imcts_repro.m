function out = run_imcts_repro(varargin)
% Reproduce improved MCTS (I-MCTS) from Ye et al. (2025)
% Usage:
%   out = run_imcts_repro();                 % default: 40 targets (paper setup)
%   out = run_imcts_repro(40);               % 40 targets
%   out = run_imcts_repro(50, 10000, 2025);  % nTargets, m, seed
%   out = run_imcts_repro(40,10000,2025,true,'fuel',0.01,'paper',0.5,5);
%   out = run_imcts_repro(40,10000,2025,true,'time',0.01,'pmdt',0.8,100);

nTargets = 40;
m = 1500;        % faster default for interactive runs
seed = 2025;
verbose = true;
objective = 'weighted'; % fuel | time | weighted
timeWeight = 2e-6;     % balances fuel (km/s) vs time (s)
dataSource = 'paper'; % paper | pmdt
dvmax = 1.5;          % km/s per action leg (less over-pruning for PMDT data)
infeasiblePenalty = []; % defaults set after source selection
if nargin >= 1, nTargets = varargin{1}; end
if nargin >= 2, m = varargin{2}; end
if nargin >= 3, seed = varargin{3}; end
if nargin >= 4, verbose = varargin{4}; end
if nargin >= 5, objective = varargin{5}; end
if nargin >= 6, timeWeight = varargin{6}; end
if nargin >= 7, dataSource = varargin{7}; end
if nargin >= 8, dvmax = varargin{8}; end
if nargin >= 9, infeasiblePenalty = varargin{9}; end

[data, ref] = load_case_data(dataSource);
maxTargets = size(data.targets, 1);
assert(nTargets >= 3 && nTargets <= maxTargets, ...
    'nTargets must be in [3, %d] for source "%s".', maxTargets, dataSource);
caseData = data;
caseData.targets = data.targets(1:nTargets, :);
caseData.targetIds = data.targetIds(1:nTargets);

model.mu = 398600.4418;   % km^3/s^2
model.Re = 6378.137;      % km
model.J2 = 1.08262668e-3;
model.dvPhaseWeight = 0.05;
model.timePhaseWeight = 0.20;
model.transferEta = 1.20;
model.kRangeEq4 = -20:20;
if strcmpi(dataSource, 'pmdt')
    model.maxDriftTimeSec = 2000 * 86400; % align with PMDT-style wait horizon
else
    model.maxDriftTimeSec = inf;
end

% Keep same calibration convention as HEGA code for comparability.
calib = struct('fuelScale', 1.0, 'timeScale', 1.0);
if strcmpi(dataSource, 'paper') && nTargets >= 40
    ref40 = ref.hega40;
    cd40 = caseData;
    cd40.targets = data.targets(1:40, :);
    cd40.targetIds = data.targetIds(1:40);
    [f0, t0] = evaluate_solution(ref40.ssc1, ref40.ssc2, cd40, model, calib);
    calib.fuelScale = ref40.fuel / f0;
    calib.timeScale = ref40.time / t0;
end

params.m = m;
params.nchild = 10;      % Table 9
params.dvmax = dvmax;
params.cUct = 0.8;
params.rngSeed = seed;
params.verbose = verbose;
params.rolloutRoulette = false; % Greedy rollout is usually more stable for fuel minimization
if isempty(infeasiblePenalty)
    if strcmpi(dataSource, 'pmdt')
        params.infeasiblePenalty = 500.0;
    else
        params.infeasiblePenalty = 5.0;
    end
else
    params.infeasiblePenalty = infeasiblePenalty;
end
params.objective = objective;
params.timeWeight = timeWeight;
params.strictCompleteFinal = false; % report real values even if incomplete; completeness is flagged separately

t0 = tic;
best = imcts_optimize(caseData, model, params, calib);
computeTimeSec = toc(t0);

if verbose
    fprintf('\n=== I-MCTS Reproduction (%d targets) ===\n', nTargets);
    fprintf('Data source: %s | Objective: %s', dataSource, objective);
    if strcmpi(objective, 'weighted')
        fprintf(' (timeWeight=%.4f)', timeWeight);
    end
    fprintf('\n');
    fprintf('Best fuel (km/s): %.4f\n', best.fuel);
    fprintf('Compute time (s): %.3f\n', computeTimeSec);
    fprintf('Mission transfer time (scaled s): %.3f\n', best.time_sec_scaled);
    fprintf('Mission transfer time (scaled d): %.3f\n', best.time_days_scaled);
    fprintf('Mission transfer time (Eq.4 sum, d): %.3f\n', best.tof_days_eq4);
    fprintf('SSc1 sequence (%d targets): [%s]\n', numel(best.seq{1}), fmt_seq(best.seq{1}));
    fprintf('SSc2 sequence (%d targets): [%s]\n', numel(best.seq{2}), fmt_seq(best.seq{2}));
    allAssigned = [best.seq{1}(:); best.seq{2}(:)];
    fprintf('Assigned total targets: %d\n', numel(unique(allAssigned)));
    if isfield(best, 'complete') && ~best.complete
        fprintf('WARNING: Incomplete assignment (%d/%d). Missing IDs: [%s]\n', ...
            best.assignedCount, best.totalTargets, fmt_seq(best.missingIds));
    end

    if nTargets == 40 && strcmpi(dataSource, 'paper')
        fprintf('\nPaper Table 11 I-MCTS: fuel=1.2352 km/s, time=126.535 s\n');
        fprintf('Paper SSc1: [30 6 16 28 33 24 20 19 2 11 26 29 18 1 12 14 39 25 17 31 21 7 22 35 23 37 13 40 34 9 5 32]\n');
        fprintf('Paper SSc2: [36 3 4 27 15 10 8 38]\n');
    end
end

out.best = best;
out.params = params;
out.model = model;
out.calibration = calib;
out.caseData = caseData;
out.compute_time_s = computeTimeSec;
out.total_tof_days_eq4 = best.tof_days_eq4;
out.total_tof_sec_eq4 = best.tof_sec_eq4;
out.total_tof_days_scaled = best.time_days_scaled;
out.total_tof_sec_scaled = best.time_sec_scaled;
end

function s = fmt_seq(v)
if isempty(v)
    s = '';
    return;
end
s = strjoin(string(v), ' ');
end
