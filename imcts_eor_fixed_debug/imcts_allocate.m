function [best, tree] = imcts_allocate(ssc0, targets0, cfg, P)
% imcts_allocate.m
% Full IMCTS-style implementation using an indexed node array.
%
% Implements the paper structure:
%   1. Selection with roulette wheel based on node value
%   2. Expansion with dynamic classification Eq. (16)-(18)
%   3. Simulation / rollout using fuel-reward criterion
%   4. Backpropagation using best value from child paths Eq. (21)-style
%   5. Final decision by repeatedly selecting most valuable child
%
% Inputs:
%   ssc0     : Ns x 6 initial servicing spacecraft orbital elements
%   targets0 : Nt x 6 target orbital elements
%   cfg      : cfg.Ns, cfg.Nt, cfg.m, cfg.nchild, cfg.Dvmax, cfg.NmaxPhase
%   P        : constants
%
% Output:
%   best.seq, best.totalDV, best.actions
%   tree: internal node array

Ns = cfg.Ns;
Nt = cfg.Nt;

if size(ssc0,1) ~= Ns
    error('cfg.Ns does not match number of servicing spacecraft rows.');
end

if size(targets0,1) ~= Nt
    targets0 = targets0(1:Nt,:);
end

% ---------------- Root node ----------------
root = make_node(ssc0, zeros(Ns,1), false(1,Nt), cell(Ns,1), [], 0, 0, []);
tree = root;

fprintf('Starting full IMCTS: Ns=%d, Nt=%d, m=%d, nchild=%d\n', Ns, Nt, cfg.m, cfg.nchild);
drawnow;

% ---------------- Build search tree ----------------
for iter = 1:cfg.m

    % 1) Selection
    leafIdx = select_node(tree, 1);

    % 2) Expansion
    if tree(leafIdx).n > 0 && ~all(tree(leafIdx).assigned)
        [tree, leafIdx] = expand_node(tree, leafIdx, targets0, cfg, P);
    elseif tree(leafIdx).n == 0
        % keep leaf as selected for first visit
    elseif all(tree(leafIdx).assigned)
        % terminal leaf
    end

    % 3) Simulation / rollout
    rolloutCost = simulate_rollout(tree(leafIdx), targets0, cfg, P);

    if ~isfinite(rolloutCost) || rolloutCost <= 0
        qLeaf = 0;
    else
        qLeaf = 1 / rolloutCost;
    end

    % 4) Backpropagation
    tree = backpropagate(tree, leafIdx, qLeaf);

    if isfield(cfg,'verbose') && cfg.verbose
        if mod(iter,50) == 0 || iter == 1
            fprintf('IMCTS iter %d/%d | nodes=%d | root Q=%.6g | root visits=%d\n', ...
                iter, cfg.m, numel(tree), tree(1).Q, tree(1).n);
            drawnow;
        end
    end
end

% ---------------- Extract best path ----------------
idx = 1;
actions = [];
while ~isempty(tree(idx).children)
    childIdxs = tree(idx).children;
    Qvals = arrayfun(@(ii) tree(ii).Q, childIdxs);

    % tie-breaker: lower step cost
    [maxQ,~] = max(Qvals);
    cand = childIdxs(abs(Qvals - maxQ) < 1e-15);
    if numel(cand) > 1
        stepCosts = arrayfun(@(ii) tree(ii).stepCost, cand);
        [~,ii] = min(stepCosts);
        nextIdx = cand(ii);
    else
        [~,ii] = max(Qvals);
        nextIdx = childIdxs(ii);
    end

    actions = [actions; tree(nextIdx).action]; %#ok<AGROW>
    idx = nextIdx;

    if all(tree(idx).assigned)
        break;
    end
end

best.seq = tree(idx).seq;
best.actions = actions;
best.totalDV = sequence_cost_direct(ssc0, targets0, best.seq, P, cfg);
best.finalNode = idx;

fprintf('\nFinished IMCTS.\n');
fprintf('Assigned targets = %d/%d\n', nnz(tree(idx).assigned), Nt);
fprintf('Total fuel = %.6f km/s\n', best.totalDV);
for k = 1:Ns
    fprintf('SSc%d: ', k);
    disp(best.seq{k});
end

end

% =========================================================================
function node = make_node(ssc, tnow, assigned, seq, parent, action, stepCost, depth)

if isempty(depth)
    depth = 0;
end

node = struct();
node.ssc = ssc;
node.tnow = tnow;
node.assigned = assigned;
node.seq = seq;
node.parent = parent;
node.action = action;
node.stepCost = stepCost;
node.depth = depth;

node.children = [];
node.untriedActions = [];
node.expanded = false;

node.Q = 0;          % value = reciprocal of best total future fuel
node.n = 0;          % visits
node.bestCost = inf; % best rollout cost from this node

end

% =========================================================================
function idx = select_node(tree, idx)

while tree(idx).expanded && ~isempty(tree(idx).children) && isempty(tree(idx).untriedActions)

    childIdxs = tree(idx).children;
    Q = arrayfun(@(ii) tree(ii).Q, childIdxs);

    % Roulette wheel selection. Shift if all zeros.
    Q = max(Q, 0);
    if sum(Q) <= 0
        probs = ones(size(Q))/numel(Q);
    else
        probs = Q/sum(Q);
    end

    r = rand;
    cs = cumsum(probs);
    j = find(r <= cs, 1, 'first');
    idx = childIdxs(j);
end

end

% =========================================================================
function [tree, childIdx] = expand_node(tree, idx, targets, cfg, P)

if all(tree(idx).assigned)
    childIdx = idx;
    return;
end

if ~tree(idx).expanded
    actions = generate_actions_dynamic(tree(idx), targets, cfg, P);
    tree(idx).untriedActions = actions;
    tree(idx).expanded = true;
end

if isempty(tree(idx).untriedActions)
    % no expansion possible
    childIdx = idx;
    return;
end

% Pick one untried action. Paper says randomly select child after expansion.
pick = randi(size(tree(idx).untriedActions,1));
act = tree(idx).untriedActions(pick,:);
tree(idx).untriedActions(pick,:) = [];

[childNode, ok] = transition(tree(idx), act, targets, cfg, P);

if ~ok
    childIdx = idx;
    return;
end

childNode.parent = idx;
childNode.depth = tree(idx).depth + 1;

tree(end+1) = childNode; %#ok<AGROW>
childIdx = numel(tree);
tree(idx).children(end+1) = childIdx;

end

% =========================================================================
function actions = generate_actions_dynamic(node, targets, cfg, P)

Ns = cfg.Ns;
unassigned = find(~node.assigned);

if isempty(unassigned)
    actions = [];
    return;
end

subsets = cell(Ns,1);

% Eq. (16)-(17): dynamic classification
for jj = unassigned

    phi = inf(Ns,1);
    dvv = inf(Ns,1);

    for k = 1:Ns
        phi(k) = sqrt(((node.ssc(k,1) - targets(jj,1))/(2*node.ssc(k,1)))^2 + ...
                      (node.ssc(k,3) - targets(jj,3))^2);

        dt = transfer_time(node.ssc(k,:), targets(jj,:), P);
        dvv(k) = eor_estimate(node.ssc(k,:), targets(jj,:), dt, P, cfg.NmaxPhase);
    end

    [~, ii] = min(phi);

    % Paper condition includes dv <= Dvmax. For numerical robustness, if
    % Dvmax is missing or <=0, keep classification only.
    if ~isfield(cfg,'Dvmax') || cfg.Dvmax <= 0 || dvv(ii) <= cfg.Dvmax
        subsets{ii}(end+1) = jj; %#ok<AGROW>
    end
end

% If strict Dvmax removed everything, fallback to nearest-orbit classification.
if all(cellfun(@isempty, subsets))
    for jj = unassigned
        phi = inf(Ns,1);
        for k = 1:Ns
            phi(k) = sqrt(((node.ssc(k,1) - targets(jj,1))/(2*node.ssc(k,1)))^2 + ...
                          (node.ssc(k,3) - targets(jj,3))^2);
        end
        [~,ii] = min(phi);
        subsets{ii}(end+1) = jj; %#ok<AGROW>
    end
end

% Eq. (18): ai in subset_i or ai = no target
C = cell(1,Ns);
for k = 1:Ns
    C{k} = [0 subsets{k}]; % 0 = no target
end

grids = cell(1,Ns);
[grids{:}] = ndgrid(C{:});

actions = zeros(numel(grids{1}), Ns);
for k = 1:Ns
    actions(:,k) = grids{k}(:);
end

% remove all-zero action
actions(all(actions == 0, 2), :) = [];

% remove duplicate targets in same action
keep = true(size(actions,1),1);
for r = 1:size(actions,1)
    nz = actions(r, actions(r,:) > 0);
    if numel(unique(nz)) < numel(nz)
        keep(r) = false;
    end
end
actions = actions(keep,:);

if isempty(actions)
    return;
end

% Sort actions by immediate reward, keep top nchild
costs = inf(size(actions,1),1);
for r = 1:size(actions,1)
    costs(r) = action_cost(node.ssc, actions(r,:), targets, cfg, P);
end

[~, ord] = sort(costs, 'ascend');
nkeep = min(cfg.nchild, numel(ord));
actions = actions(ord(1:nkeep), :);

end

% =========================================================================
function [child, ok] = transition(node, act, targets, cfg, P)

Ns = cfg.Ns;
ssc = node.ssc;
tnow = node.tnow;
assigned = node.assigned;
seq = node.seq;

stepCost = 0;
ok = true;

if all(act == 0)
    ok = false;
    child = node;
    return;
end

for k = 1:Ns
    j = act(k);

    if j > 0
        if assigned(j)
            ok = false;
            child = node;
            return;
        end

        dt = transfer_time(ssc(k,:), targets(j,:), P);
        dv = eor_estimate(ssc(k,:), targets(j,:), dt, P, cfg.NmaxPhase);

        if ~isfinite(dv) || dv <= 0
            ok = false;
            child = node;
            return;
        end

        stepCost = stepCost + dv;
        tnow(k) = tnow(k) + dt;

        % After rendezvous, SSc state follows target at arrival time.
        ssc(k,:) = propagate_oe_j2(targets(j,:), tnow(k), P);

        assigned(j) = true;
        seq{k}(end+1) = j;
    end
end

child = make_node(ssc, tnow, assigned, seq, [], act, stepCost, []);

end

% =========================================================================
function c = action_cost(ssc, act, targets, cfg, P)

c = 0;

for k = 1:cfg.Ns
    j = act(k);
    if j > 0
        dt = transfer_time(ssc(k,:), targets(j,:), P);
        dv = eor_estimate(ssc(k,:), targets(j,:), dt, P, cfg.NmaxPhase);
        if ~isfinite(dv) || dv <= 0
            c = inf;
            return;
        end
        c = c + dv;
    end
end

if c <= 0
    c = inf;
end

end

% =========================================================================
function rolloutCost = simulate_rollout(node, targets, cfg, P)

tmp = node;
rolloutCost = 0;

guard = 0;
while ~all(tmp.assigned)
    guard = guard + 1;
    if guard > cfg.Nt + 5
        rolloutCost = inf;
        return;
    end

    actions = generate_actions_dynamic(tmp, targets, cfg, P);

    if isempty(actions)
        rolloutCost = inf;
        return;
    end

    % Simulation in improved MCTS should use fuel consumption as value.
    % Choose probabilistically weighted by reciprocal immediate cost.
    costs = inf(size(actions,1),1);
    for r = 1:size(actions,1)
        costs(r) = action_cost(tmp.ssc, actions(r,:), targets, cfg, P);
    end

    weights = 1 ./ max(costs, eps);
    weights(~isfinite(weights)) = 0;

    if sum(weights) <= 0
        pick = randi(size(actions,1));
    else
        weights = weights / sum(weights);
        cs = cumsum(weights);
        pick = find(rand <= cs, 1, 'first');
    end

    act = actions(pick,:);
    [next, ok] = transition(tmp, act, targets, cfg, P);

    if ~ok
        rolloutCost = inf;
        return;
    end

    rolloutCost = rolloutCost + next.stepCost;
    tmp = next;
end

end

% =========================================================================
function tree = backpropagate(tree, idx, qLeaf)

% Convert leaf value to cost.
if qLeaf <= 0
    leafCost = inf;
else
    leafCost = 1/qLeaf;
end

while true

    % Eq. (21)-style: node value stores reciprocal of minimum found cost.
    if leafCost < tree(idx).bestCost
        tree(idx).bestCost = leafCost;
        tree(idx).Q = 1 / max(leafCost, eps);
    end

    tree(idx).n = tree(idx).n + 1;

    if isempty(tree(idx).parent)
        break;
    end

    parentIdx = tree(idx).parent;
    leafCost = tree(idx).stepCost + leafCost;
    idx = parentIdx;
end

end

% =========================================================================
function J = sequence_cost_direct(ssc0, targets, seq, P, cfg)

J = 0;
ssc = ssc0;
tnow = zeros(size(ssc0,1),1);

for k = 1:numel(seq)
    for j = seq{k}
        dt = transfer_time(ssc(k,:), targets(j,:), P);
        dv = eor_estimate(ssc(k,:), targets(j,:), dt, P, cfg.NmaxPhase);

        J = J + dv;
        tnow(k) = tnow(k) + dt;
        ssc(k,:) = propagate_oe_j2(targets(j,:), tnow(k), P);
    end
end

end
