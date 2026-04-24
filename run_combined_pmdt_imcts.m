function out = run_combined_pmdt_imcts(varargin)

% Usage:
%   out = run_combined_pmdt_imcts();
%   out = run_combined_pmdt_imcts(40, 10000, 2025, true, 'weighted', 5e-9, 1.0, 100);
%   out = run_combined_pmdt_imcts(40, 10000, 2025, true, 'time');
%   out = run_combined_pmdt_imcts(40, 10000, 2025, true, 'fuel', 0.01, 2.0, 100, 3, 5);
%
% Inputs:
%   nTargets  - number of debris targets for IMCTS (default 40)
%   m         - iterations per decision step (default 10000)
%   seed      - RNG seed (default 2025)
%   verbose   - print progress (default true)
%   objective - fuel | time | weighted (default fuel)
%   timeWeight- weight for weighted objective (default 0.01)
%   dvmax     - max per-leg delta-v in km/s (default 0.5)
%   infeasiblePenalty - penalty per unassigned target in rollout (default 100 for PMDT)
%   refineTopK - number of top IMCTS candidates to PMDT-refine (default 3, 0 to disable)
%   nSeedRuns  - number of IMCTS seed runs used to build candidate pool (default 5)
%   maxLegTofDays - strict Eq.4 per-leg TOF cap for TOF scoring (default 365)

nTargets = 40;
m = 1500;
seed = 2025;
verbose = true;
objective = 'weighted';
timeWeight = 5e-9;
dvmax = 1.2;             % allow more direct transfers (often lowers TOF)
infeasiblePenalty = 100;
refineTopK = 0; % IMCTS-only by default (fast)
nSeedRuns = 8;
maxLegTofDays = 120;
if nargin >= 1, nTargets = varargin{1}; end
if nargin >= 2, m = varargin{2}; end
if nargin >= 3, seed = varargin{3}; end
if nargin >= 4, verbose = varargin{4}; end
if nargin >= 5, objective = varargin{5}; end
if nargin >= 6, timeWeight = varargin{6}; end
if nargin >= 7, dvmax = varargin{7}; end
if nargin >= 8, infeasiblePenalty = varargin{8}; end
if nargin >= 9, refineTopK = varargin{9}; end
if nargin >= 10, nSeedRuns = varargin{10}; end
if nargin >= 11, maxLegTofDays = varargin{11}; end

root = fileparts(mfilename('fullpath'));

% PMDT paths (copied project)
addpath(root);
addpath(fullfile(root, 'PMDT'));
addpath(fullfile(root, 'CoordinateConversions'));
addpath(fullfile(root, 'Pertubations'));
addpath(fullfile(root, 'Data'));
addpath(fullfile(root, 'SGP4routines_NAIF'));

% IMCTS paths (added into this new copy)
addpath(fullfile(root, 'IMCTS'));

if verbose
    fprintf('\n=== Combined PMDT + IMCTS Workspace ===\n');
    fprintf('Root: %s\n', root);
    fprintf('Running IMCTS on PMDT debris data (DebrisData.mat)...\n');
    fprintf('Objective=%s | timeWeight=%.3e | dvmax=%.2f | refineTopK=%d\n', ...
        objective, timeWeight, dvmax, refineTopK);
    fprintf('Strict Eq.4 leg TOF cap = %.1f d\n', maxLegTofDays);
end

% IMCTS run using PMDT debris dataset through load_case_data('pmdt').
out = run_imcts_repro(nTargets, m, seed, verbose, objective, timeWeight, 'pmdt', dvmax, infeasiblePenalty);
if refineTopK > 0
    out.best.strict_tof_days_eq4 = compute_strict_tof_days(out.best.seq, out.caseData, maxLegTofDays);
end

% Hybrid stage: PMDT-refine top IMCTS candidates (multi-seed pool).
if refineTopK > 0
    out.hybrid = hybrid_refine_with_pmdt(root, out, nTargets, m, seed, ...
        objective, timeWeight, dvmax, infeasiblePenalty, refineTopK, nSeedRuns, verbose, maxLegTofDays);
end
end

function hybrid = hybrid_refine_with_pmdt(root, baseOut, nTargets, m, seed, ...
    objective, timeWeight, dvmax, infeasiblePenalty, refineTopK, nSeedRuns, verbose, maxLegTofDays)

nSeedRuns = max(1, round(nSeedRuns));
refineTopK = max(1, round(refineTopK));

pool = cell(nSeedRuns, 1);
pool{1} = baseOut;
for k = 2:nSeedRuns
    sk = seed + (k - 1);
    pool{k} = run_imcts_repro(nTargets, m, sk, false, objective, timeWeight, ...
        'pmdt', dvmax, infeasiblePenalty);
end

fuelVals = cellfun(@(s) s.best.fuel, pool);
[~, ord] = sort(fuelVals, 'ascend');
pick = ord(1:min(refineTopK, numel(ord)));

if verbose
    fprintf('\n=== Hybrid Refinement (PMDT over top IMCTS candidates) ===\n');
    fprintf('Candidate pool: %d seed runs | PMDT refine top: %d\n', nSeedRuns, numel(pick));
end

rescored = repmat(struct(), numel(pick), 1);
for j = 1:numel(pick)
    idx = pick(j);
    cand = pool{idx};
    [pmFuel, pmTimeDays, strictTofDays, ok, msg] = pmdt_rescore_solution( ...
        root, cand.best, cand.caseData, objective, timeWeight, maxLegTofDays);

    rescored(j).poolIndex = idx;
    rescored(j).seed = cand.params.rngSeed;
    rescored(j).imctsFuelKmps = cand.best.fuel;
    rescored(j).imctsTimeNorm = cand.best.time;
    rescored(j).imctsComplete = cand.best.complete;
    rescored(j).pmdtFuelKmps = pmFuel;
    rescored(j).pmdtTimeDays = pmTimeDays;
    rescored(j).strictTofDaysEq4 = strictTofDays;
    rescored(j).pmdtOK = ok;
    rescored(j).message = msg;
    rescored(j).best = cand.best;
end

% Select by PMDT score if feasible, else fall back to IMCTS best.
score = inf(numel(rescored), 1);
for j = 1:numel(rescored)
    if ~rescored(j).pmdtOK
        continue;
    end
    tDays = rescored(j).strictTofDaysEq4;
    if ~isfinite(tDays)
        tDays = rescored(j).pmdtTimeDays;
    end
    switch lower(objective)
        case 'fuel'
            score(j) = rescored(j).pmdtFuelKmps;
        case 'time'
            score(j) = tDays * 86400.0;
        case 'weighted'
            score(j) = rescored(j).pmdtFuelKmps + timeWeight * tDays * 86400.0;
        otherwise
            score(j) = rescored(j).pmdtFuelKmps;
    end
end

[bestScore, ib] = min(score);
hybrid = struct();
hybrid.rescored = rescored;
hybrid.selectedBy = 'imcts';
hybrid.selected = baseOut.best;
hybrid.selectedPoolIndex = 1;
hybrid.selectedScore = baseOut.best.fuel;

if isfinite(bestScore)
    chosen = rescored(ib);
    hybrid.selectedBy = 'pmdt_refine';
    hybrid.selected = chosen.best;
    hybrid.selected.strict_tof_days_eq4 = chosen.strictTofDaysEq4;
    hybrid.selectedPoolIndex = chosen.poolIndex;
    hybrid.selectedSeed = chosen.seed;
    hybrid.selectedScore = bestScore;
end

if verbose
    for j = 1:numel(rescored)
        fprintf('cand %d seed=%d | IMCTS fuel=%.4f | PMDT fuel=%.4f km/s, time=%.2f d | strictTOF=%.2f d | ok=%d\n', ...
            j, rescored(j).seed, rescored(j).imctsFuelKmps, rescored(j).pmdtFuelKmps, ...
            rescored(j).pmdtTimeDays, rescored(j).strictTofDaysEq4, rescored(j).pmdtOK);
    end
    fprintf('Hybrid selected by: %s\n', hybrid.selectedBy);
    if isfield(hybrid.selected, 'strict_tof_days_eq4')
        fprintf('Hybrid selected strict Eq.4 TOF: %.2f d\n', hybrid.selected.strict_tof_days_eq4);
    end
end
end

function [fuelKmps, timeDays, strictTofDays, ok, msg] = pmdt_rescore_solution(root, best, caseData, objective, ~, maxLegTofDays)
fuelKmps = inf;
timeDays = inf;
strictTofDays = inf;
ok = false;
msg = '';

try
    strictTofDays = compute_strict_tof_days(best.seq, caseData, maxLegTofDays);

    % PMDT model defaults (match PMDT app conventions; conservative settings).
    pBase = struct();
    pBase.TU = 86400;                     % s
    pBase.mu = 3.986005e14;               % m^3/s^2
    pBase.Re = 6378137;                   % m
    pBase.J2 = 1.0826668e-3;
    pBase.N = 100;
    pBase.k = 3 * pBase.J2 * pBase.Re^2 / (2 * pBase.mu^3);
    pBase.g0 = 9.80665;

    pBase.TOFbefore = 0;
    pBase.m_wet_servicer = 400;           % kg
    pBase.T = 0.5;                        % N
    pBase.Isp = 3100;                     % s
    pBase.target_altitude = 360e3;        % m (requested shepherd orbit)
    pBase.dvLimit = 8000;                 % m/s
    pBase.waitTimeLimit = 365;            % days (lower TOF bias)
    pBase.a_pchangeLimit = 0.01;
    pBase.dutyRatio = 0.9;
    pBase.Cd = 2.2;
    pBase.Area = 2.0;
    pBase.RDV1 = 0 * pBase.TU;
    pBase.RDV2 = 0 * pBase.TU;
    pBase.krange = -20:1:20;
    pBase.dragfactor = 1;
    pBase.eclipses = false;
    pBase.drag = false;
    pBase.plots = false;
    pBase.guidance = false;
    pBase.t0 = juliandate(datetime(2022, 3, 25, 6, 37, 13));
    pBase.Topt = strcmpi(objective, 'time');

    names = map_target_ids_to_names(best.seq, caseData, root);

    dvTotMps = 0;
    tofTotDays = 0;
    for s = 1:numel(names)
        debrisID = names{s};
        if isempty(debrisID)
            continue;
        end

        param = pBase;
        param.debrisID = debrisID;
        [x0, ub, lb] = getInitGuess(param);
        if isempty(x0)
            ptxt = evalc('[~, ~, tofDays, dvMps] = PMDT([], debrisID, param);'); %#ok<NASGU>
        else
            opts = optimoptions('fmincon', 'Display', 'off', ...
                'MaxIterations', 40, 'MaxFunctionEvaluations', 4000);
            ptxt = evalc('[xOpt, ~] = fmincon(@(x) PMDT(x, debrisID, param), x0, [], [], [], [], lb, ub, @(x) PMDT_limits(x, debrisID, param), opts);'); %#ok<NASGU>
            ptxt = evalc('[~, ~, tofDays, dvMps] = PMDT(xOpt, debrisID, param);'); %#ok<NASGU>
        end
        dvTotMps = dvTotMps + dvMps;
        tofTotDays = tofTotDays + tofDays;
    end

    fuelKmps = dvTotMps / 1000.0;
    timeDays = tofTotDays;
    ok = isfinite(fuelKmps) && isfinite(timeDays);
    if ~ok
        msg = 'Non-finite PMDT score';
    end
catch ME
    msg = ME.message;
end
end

function names = map_target_ids_to_names(seqCell, caseData, root)
names = cell(size(seqCell));

if isfield(caseData, 'targetNames') && ~isempty(caseData.targetNames)
    tgtNames = caseData.targetNames;
else
    S = load(fullfile(root, 'Data', 'DebrisData.mat'), 'DebrisData');
    tgtNames = string({S.DebrisData.name})';
end

for s = 1:numel(seqCell)
    ids = seqCell{s};
    if isempty(ids)
        names{s} = {};
    else
        names{s} = cellstr(tgtNames(ids));
    end
end
end

function tofDays = compute_strict_tof_days(seqCell, caseData, maxLegTofDays)
model.mu = 398600.4418;
model.Re = 6378.137;
model.J2 = 1.08262668e-3;
model.kRangeEq4 = -100:100;
model.minDtSec = 1e-6 * 86400.0;
if isfinite(maxLegTofDays)
    model.maxLegTofSec = maxLegTofDays * 86400.0;
else
    model.maxLegTofSec = inf;
end

Ns = numel(seqCell);
tofSec = 0.0;
for s = 1:Ns
    curOE = caseData.ssc(s, :);
    curT = 0.0;
    ids = seqCell{s};
    for k = 1:numel(ids)
        tgt = caseData.targets(ids(k), :);
        [dt, ok] = strict_eq4_leg_tof(curOE, tgt, curT, model);
        if ~ok || ~isfinite(dt)
            tofDays = inf;
            return;
        end
        tofSec = tofSec + dt;
        curT = curT + dt;
        curOE = tgt;
    end
end
tofDays = tofSec / 86400.0;
end

function [dtBest, ok] = strict_eq4_leg_tof(oe1, oe2, tNow, model)
mu = model.mu; Re = model.Re; J2 = model.J2;
a1 = Re + oe1(1); e1 = oe1(2); i1 = deg2rad(oe1(3)); O1 = deg2rad(oe1(4));
a2 = Re + oe2(1); e2 = oe2(2); i2 = deg2rad(oe2(3)); O2 = deg2rad(oe2(4));

n1 = sqrt(mu / a1^3);
n2 = sqrt(mu / a2^3);
p1 = a1 * (1 - e1^2);
p2 = a2 * (1 - e2^2);
Odot1 = -1.5 * J2 * n1 * (Re / p1)^2 * cos(i1);
Odot2 = -1.5 * J2 * n2 * (Re / p2)^2 * cos(i2);

dO0 = wrap_to_pi_local(O2 - O1);
dOdot = Odot2 - Odot1;
dtBest = inf;
ok = false;
if abs(dOdot) < 1e-14
    return;
end

for kk = model.kRangeEq4
    dO = dO0 + 2*pi*kk;
    dt = dO / dOdot;
    if ~isfinite(dt) || dt <= model.minDtSec
        continue;
    end
    if isfinite(model.maxLegTofSec) && dt > model.maxLegTofSec
        continue;
    end
    if dt < dtBest
        dtBest = dt;
        ok = true;
    end
end

end

function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end
