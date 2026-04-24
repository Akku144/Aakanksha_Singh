function best = imcts_optimize(caseData, model, params, calib)
% Improved MCTS for multi-spacecraft multi-target allocation.

rng(params.rngSeed);
Ns = size(caseData.ssc, 1);
params = set_objective_defaults(params);

nodes = init_root(caseData, Ns);
root = 1;
step = 0;

while true
    if isempty(nodes(root).unassigned) || nodes(root).dead
        break;
    end

    step = step + 1;
    for it = 1:params.m
        leaf = improved_selection(nodes, root, params);

        simNode = leaf;
        if nodes(leaf).visits > 0
            [nodes, newNode] = improved_expansion(nodes, leaf, caseData, model, params, calib);
            if newNode > 0
                simNode = newNode;
            end
        end

        [qSim, ~] = rollout_value(nodes(simNode), caseData, model, params, calib);
        nodes(simNode).visits = nodes(simNode).visits + 1;
        if qSim > nodes(simNode).Q
            nodes(simNode).Q = qSim;
        end

        nodes = improved_backprop(nodes, simNode);
    end

    ch = nodes(root).children;
    if isempty(ch)
        nodes(root).dead = true;
        break;
    end

    qvals = [nodes(ch).Q];
    [~, ib] = max(qvals);
    root = ch(ib);

    if params.verbose
        fprintf('Step %d | remaining=%d | rootQ=%.6f\n', step, numel(nodes(root).unassigned), nodes(root).Q);
    end
end

chosen = root;
isComplete = arrayfun(@(n) isempty(n.unassigned), nodes);
if any(isComplete)
    idx = find(isComplete);
    [~, ib] = max([nodes(idx).Q]);
    chosen = idx(ib);
end

best.seq = nodes(chosen).seq;
best.seq_imcts = best.seq;
strictFinal = false;
if isfield(params, 'strictCompleteFinal') && ~isempty(params.strictCompleteFinal)
    strictFinal = logical(params.strictCompleteFinal);
end
[~, ~, evalMeta0] = evaluate_solution(best.seq{1}, best.seq{2}, caseData, model, calib, false);
best.repairedComplete = false;
if isfield(params, 'forceCompleteRepair') && params.forceCompleteRepair && ~evalMeta0.complete
    best.seq = repair_complete_greedy(best.seq, caseData, model, calib, params);
    best.repairedComplete = true;
end

[best.fuel, best.time, evalMeta] = evaluate_solution(best.seq{1}, best.seq{2}, caseData, model, calib, strictFinal);
best.tof_sec_eq4 = evalMeta.timeRawSec;
best.tof_days_eq4 = evalMeta.timeRawSec / 86400.0;
best.time_sec_scaled = evalMeta.timeScaledSec;
best.time_days_scaled = evalMeta.timeScaledSec / 86400.0;
best.fuel_kmps_raw = evalMeta.fuelRawKmps;
best.fuel_kmps_scaled = evalMeta.fuelScaledKmps;
best.node = nodes(chosen);
best.complete = evalMeta.complete;
best.assignedCount = evalMeta.assignedCount;
best.totalTargets = evalMeta.totalTargets;
best.missingCount = evalMeta.missingCount;
best.missingIds = evalMeta.missingIds;
end

function seqOut = repair_complete_greedy(seqIn, caseData, model, calib, params)
seqOut = {seqIn{1}(:)', seqIn{2}(:)'};
Ns = numel(seqOut);
assigned = unique([seqOut{1}, seqOut{2}]);
missing = setdiff(caseData.targetIds(:)', assigned, 'stable');
if isempty(missing)
    return;
end

% Cost-aware completion: place each missing target at the best insertion point.
while ~isempty(missing)
    bestTotalCost = inf;
    bestS = 1;
    bestIdx = 1;
    bestPos = 1;
    for mi = 1:numel(missing)
        tid = missing(mi);
        for s = 1:Ns
            L = numel(seqOut{s});
            for pos = 1:(L + 1)
                cand = seqOut;
                cand{s} = [cand{s}(1:pos-1), tid, cand{s}(pos:end)];
                cTot = full_objective_cost(cand, caseData, model, calib, params);
                if cTot < bestTotalCost
                    bestTotalCost = cTot;
                    bestS = s;
                    bestIdx = mi;
                    bestPos = pos;
                end
            end
        end
    end

    pick = missing(bestIdx);
    seqOut{bestS} = [seqOut{bestS}(1:bestPos-1), pick, seqOut{bestS}(bestPos:end)];
    missing(bestIdx) = [];
end

% Local improvement: relocate one target at a time if total objective drops.
seqOut = improve_complete_solution(seqOut, caseData, model, calib, params);
end

function seqOut = improve_complete_solution(seqIn, caseData, model, calib, params)
seqOut = seqIn;
Ns = numel(seqOut);
maxPasses = 2;
tol = 1e-10;

for pass = 1:maxPasses
    improved = false;
    baseCost = full_objective_cost(seqOut, caseData, model, calib, params);
    allIds = [seqOut{1}, seqOut{2}];

    for ii = 1:numel(allIds)
        tid = allIds(ii);
        [sFrom, posFrom] = find_target(seqOut, tid);
        if sFrom == 0
            continue;
        end

        removed = seqOut;
        removed{sFrom}(posFrom) = [];

        bestCost = baseCost;
        bestSeq = seqOut;
        for sTo = 1:Ns
            L = numel(removed{sTo});
            for posTo = 1:(L + 1)
                cand = removed;
                cand{sTo} = [cand{sTo}(1:posTo-1), tid, cand{sTo}(posTo:end)];
                cTot = full_objective_cost(cand, caseData, model, calib, params);
                if cTot + tol < bestCost
                    bestCost = cTot;
                    bestSeq = cand;
                end
            end
        end

        if bestCost + tol < baseCost
            seqOut = bestSeq;
            improved = true;
            break;
        end
    end

    if ~improved
        break;
    end
end
end

function cTot = full_objective_cost(seqs, caseData, model, calib, params)
fuelTot = 0;
timeTot = 0;
Ns = numel(seqs);
for s = 1:Ns
    curOE = caseData.ssc(s, :);
    curT = 0;
    ids = seqs{s};
    for k = 1:numel(ids)
        tgt = caseData.targets(ids(k), :);
        [dv, dt] = transfer_cost(curOE, tgt, curT, model);
        dv = dv * calib.fuelScale;
        dt = dt * calib.timeScale;
        fuelTot = fuelTot + dv;
        timeTot = timeTot + dt;
        curOE = tgt;
        curT = curT + dt;
    end
end
cTot = objective_cost(fuelTot, timeTot, params);
end

function [sIdx, pIdx] = find_target(seqs, tid)
sIdx = 0;
pIdx = 0;
for s = 1:numel(seqs)
    p = find(seqs{s} == tid, 1, 'first');
    if ~isempty(p)
        sIdx = s;
        pIdx = p;
        return;
    end
end
end

function nodes = init_root(caseData, Ns)
node.parent = 0;
node.children = [];
node.action = zeros(1, Ns);
node.edgeFuel = 0;
node.curOE = caseData.ssc;
node.curTime = zeros(Ns, 1);
node.unassigned = caseData.targetIds(:)';
node.seq = cell(1, Ns);
for i = 1:Ns, node.seq{i} = []; end
node.Q = 0;
node.visits = 0;
node.expanded = false;
node.dead = false;
nodes = node;
end

function leaf = improved_selection(nodes, root, params)
leaf = root;
while true
    ch = nodes(leaf).children;
    if isempty(ch)
        return;
    end

    if ~nodes(leaf).expanded
        leaf = uct_pick(nodes, leaf, params.cUct);
    else
        leaf = roulette_pick(nodes, ch);
    end
end
end

function idx = uct_pick(nodes, parent, c)
ch = nodes(parent).children;
np = max(nodes(parent).visits, 1);
score = -inf(size(ch));
for i = 1:numel(ch)
    ci = ch(i);
    if nodes(ci).visits == 0
        score(i) = inf;
    else
        score(i) = nodes(ci).Q + c * sqrt(log(np + 1) / nodes(ci).visits);
    end
end
[~, j] = max(score);
idx = ch(j);
end

function idx = roulette_pick(nodes, ch)
q = [nodes(ch).Q];
q(~isfinite(q) | q < 0) = 0;
s = sum(q);
if s <= 0
    idx = ch(randi(numel(ch)));
    return;
end
p = q / s;
r = rand;
c = cumsum(p);
idx = ch(find(r <= c, 1, 'first'));
end

function [nodes, newChild] = improved_expansion(nodes, leaf, caseData, model, params, calib)
newChild = 0;
if nodes(leaf).expanded
    ch = nodes(leaf).children;
    if ~isempty(ch), newChild = ch(randi(numel(ch))); end
    return;
end

state = nodes(leaf);
    [actions, actionCost] = generate_actions(state, caseData, model, params, calib);

if isempty(actions)
    nodes(leaf).expanded = true;
    nodes(leaf).dead = true;
    return;
end

    [~, ord] = sort(actionCost, 'ascend');
    ord = ord(1:min(params.nchild, numel(ord)));
    actions = actions(ord, :);
    actionCost = actionCost(ord);

childIdx = zeros(1, size(actions, 1));
for r = 1:size(actions, 1)
    [childState, ok] = apply_action(state, actions(r, :), caseData, model, params, calib);
    if ~ok
        continue;
    end
    childState.parent = leaf;
    childState.children = [];
    childState.action = actions(r, :);
    childState.edgeFuel = actionCost(r);
    childState.Q = 0;
    childState.visits = 0;
    childState.expanded = false;
    childState.dead = false;

    nodes(end + 1) = childState; %#ok<AGROW>
    childIdx(r) = numel(nodes);
end

childIdx = childIdx(childIdx > 0);
nodes(leaf).children = childIdx;
nodes(leaf).expanded = true;

if isempty(childIdx)
    nodes(leaf).dead = true;
else
    newChild = childIdx(randi(numel(childIdx)));
end
end

function [actions, actionCost] = generate_actions(state, caseData, model, params, calib)
Ns = size(state.curOE, 1);
ua = state.unassigned;
if isempty(ua)
    actions = [];
    actionCost = [];
    return;
end

subsets = cell(1, Ns);
for i = 1:Ns, subsets{i} = []; end

for t = ua
    phi = inf(1, Ns);
    dv = inf(1, Ns);
    tgt = caseData.targets(t, :);
    for k = 1:Ns
        dvk = transfer_cost(state.curOE(k, :), tgt, state.curTime(k), model);
        dv(k) = dvk * calib.fuelScale;

        ak = state.curOE(k, 1) + model.Re;
        aj = tgt(1) + model.Re;
        ik = deg2rad(state.curOE(k, 3));
        ij = deg2rad(tgt(3));
        phi(k) = ((ak - aj) / (2 * ak))^2 + (ik - ij)^2;
    end

    feasible = find(dv <= params.dvmax);
    if isempty(feasible)
        continue;
    end

    % Preserve the phi-based preference while not discarding targets that
    % are feasible for another spacecraft.
    [~, ordFeas] = sort(phi(feasible), 'ascend');
    feasOrdered = feasible(ordFeas);
    for kk = feasOrdered
        subsets{kk}(end + 1) = t; %#ok<AGROW>
    end
end

actions = build_action_product(subsets);
if isempty(actions)
    actionCost = [];
    return;
end

keep = true(size(actions, 1), 1);
for r = 1:size(actions, 1)
    nz = actions(r, actions(r, :) > 0);
    if isempty(nz)
        keep(r) = false;
    elseif numel(unique(nz)) < numel(nz)
        keep(r) = false;
    end
end
actions = actions(keep, :);

if isempty(actions)
    actionCost = [];
    return;
end

actionCost = inf(size(actions, 1), 1);
for r = 1:size(actions, 1)
    [~, ok, ~, ~, cost] = apply_action(state, actions(r, :), caseData, model, params, calib);
    if ok
        actionCost(r) = cost;
    end
end
valid = isfinite(actionCost);
actions = actions(valid, :);
actionCost = actionCost(valid);
end

function A = build_action_product(subsets)
Ns = numel(subsets);
A = 0;
for i = 1:Ns
    opts = [0, subsets{i}];
    if i == 1
        A = opts(:);
    else
        B = zeros(size(A, 1) * numel(opts), i);
        idx = 1;
        for r = 1:size(A, 1)
            for o = opts
                B(idx, :) = [A(r, :), o];
                idx = idx + 1;
            end
        end
        A = B;
    end
end
end

function [next, ok, fuel, dtime, cost] = apply_action(state, action, caseData, model, params, calib)
next = state;
ok = true;
fuel = 0;
dtime = 0;
cost = inf;
Ns = size(state.curOE, 1);

for i = 1:Ns
    tgtId = action(i);
    if tgtId <= 0
        continue;
    end
    if ~ismember(tgtId, next.unassigned)
        ok = false;
        return;
    end

    tgt = caseData.targets(tgtId, :);
    [dv, dt] = transfer_cost(next.curOE(i, :), tgt, next.curTime(i), model);
    dv = dv * calib.fuelScale;
    dt = dt * calib.timeScale;

    if dv > params.dvmax
        ok = false;
        return;
    end

    fuel = fuel + dv;
    dtime = dtime + dt;
    next.curOE(i, :) = tgt;
    next.curTime(i) = next.curTime(i) + dt;
    next.seq{i}(end + 1) = tgtId;
    next.unassigned(next.unassigned == tgtId) = [];
end

if fuel <= 0
    ok = false;
    return;
end
cost = objective_cost(fuel, dtime, params);
end

function [q, totalCost] = rollout_value(state, caseData, model, params, calib)
cur = state;
totalCost = 0;

while ~isempty(cur.unassigned)
    [actions, actionCost] = generate_actions(cur, caseData, model, params, calib);
    if isempty(actions)
        totalCost = totalCost + params.infeasiblePenalty * numel(cur.unassigned);
        break;
    end

    if params.rolloutRoulette
        rew = 1 ./ max(actionCost, 1e-9);
        p = rew / sum(rew);
        r = rand;
        c = cumsum(p);
        idx = find(r <= c, 1, 'first');
    else
        [~, idx] = min(actionCost);
    end

    [nxt, ok, ~, ~, cost] = apply_action(cur, actions(idx, :), caseData, model, params, calib);
    if ~ok
        totalCost = totalCost + params.infeasiblePenalty * numel(cur.unassigned);
        break;
    end

    totalCost = totalCost + cost;
    cur = nxt;
end

if totalCost <= 0
    q = 1e9;
else
    q = 1 / totalCost;
end
end

function nodes = improved_backprop(nodes, fromIdx)
child = fromIdx;
parent = nodes(child).parent;

while parent > 0
    nodes(parent).visits = nodes(parent).visits + 1;

    R = 1 / max(nodes(child).edgeFuel, 1e-9);
    Qc = max(nodes(child).Q, 1e-12);
    qCand = 1 / (1 / R + 1 / Qc);

    if qCand > nodes(parent).Q
        nodes(parent).Q = qCand;
    end

    child = parent;
    parent = nodes(child).parent;
end
end

function params = set_objective_defaults(params)
if ~isfield(params, 'objective') || isempty(params.objective)
    params.objective = 'fuel';
end
if ~isfield(params, 'timeWeight') || isempty(params.timeWeight)
    params.timeWeight = 0.01;
end
if ~isfield(params, 'forceCompleteRepair') || isempty(params.forceCompleteRepair)
    params.forceCompleteRepair = true;
end
end

function cost = objective_cost(fuel, dtime, params)
switch lower(params.objective)
    case 'fuel'
        cost = fuel;
    case 'time'
        cost = dtime;
    case 'weighted'
        cost = fuel + params.timeWeight * dtime;
    otherwise
        error('Unknown objective "%s". Use fuel, time, or weighted.', params.objective);
end
end
