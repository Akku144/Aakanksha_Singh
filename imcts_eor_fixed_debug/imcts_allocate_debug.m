function [best, root] = imcts_allocate(ssc0, targets0, cfg, P)
% Simplified robust IMCTS-style allocator for debugging.
% This version avoids MATLAB struct parent-pointer issues.
% It keeps the paper logic:
% - state = all SSc current orbital elements + assigned targets
% - action = one target or zero target per SSc
% - expansion = dynamic classification using Eq. (16)-(18)
% - reward = reciprocal of delta-v
%
% It is intentionally made stable first. After this runs, we can add the
% exact tree backpropagation version.

Ns = cfg.Ns;
Nt = cfg.Nt;

assigned = false(1,Nt);
ssc = ssc0;
tnow = zeros(Ns,1);
seq = cell(Ns,1);
totalDV = 0;

fprintf('Starting allocation: Ns = %d, Nt = %d\n', Ns, Nt);
drawnow;

step = 0;

while ~all(assigned)
    step = step + 1;

    actions = generate_actions_dynamic(ssc, assigned, targets0, cfg, P);

    if isempty(actions)
        error('No actions generated. Check cfg.Ns, cfg.Nt, targets size, or eor_estimate output.');
    end

    % Score each possible action using immediate reward + short random lookahead
    scores = zeros(size(actions,1),1);
    costs  = zeros(size(actions,1),1);

    for aidx = 1:size(actions,1)
        act = actions(aidx,:);
        c = action_cost(ssc, act, targets0, P, cfg);
        costs(aidx) = c;

        if ~isfinite(c) || c <= 0
            scores(aidx) = 0;
        else
            scores(aidx) = 1/c;
        end
    end

    % choose best immediate action for now
    [~, bestIdx] = max(scores);
    act = actions(bestIdx,:);
    stepCost = costs(bestIdx);

    fprintf('\nStep %d | Assigned %d/%d | Action = ', step, nnz(assigned), Nt);
    disp(act);
    fprintf('Step cost = %.6f km/s\n', stepCost);
    drawnow;

    % apply action
    for k = 1:Ns
        j = act(k);
        if j > 0 && ~assigned(j)
            dt = transfer_time(ssc(k,:), targets0(j,:), P);
            dv = eor_estimate(ssc(k,:), targets0(j,:), dt, P, cfg.NmaxPhase);

            totalDV = totalDV + dv;
            tnow(k) = tnow(k) + dt;

            % after rendezvous, servicing spacecraft follows target orbit
            ssc(k,:) = propagate_oe_j2(targets0(j,:), tnow(k), P);

            assigned(j) = true;
            seq{k}(end+1) = j;

            fprintf('  SSc%d -> Target %d | dv = %.6f km/s | t = %.2f days\n', ...
                k, j, dv, tnow(k)/P.day);
        end
    end

    fprintf('Total DV so far = %.6f km/s\n', totalDV);
    drawnow;
end

best.seq = seq;
best.totalDV = totalDV;

root = struct();
root.seq = seq;
root.totalDV = totalDV;
root.assigned = assigned;

fprintf('\nFinished allocation.\n');
fprintf('Total DV = %.6f km/s\n', totalDV);

end

% -------------------------------------------------------------------------
function actions = generate_actions_dynamic(ssc, assigned, targets, cfg, P)

unassigned = find(~assigned);
Ns = cfg.Ns;

if isempty(unassigned)
    actions = [];
    return;
end

subsets = cell(Ns,1);

% Dynamic classification Eq. (16)-(17)
for jj = unassigned
    phi = inf(Ns,1);

    for k = 1:Ns
        phi(k) = sqrt(((ssc(k,1)-targets(jj,1))/(2*ssc(k,1)))^2 + ...
                      (ssc(k,3)-targets(jj,3))^2);
    end

    [~, ii] = min(phi);

    % Do not reject by Dvmax in debug mode; otherwise tree can become empty.
    subsets{ii}(end+1) = jj;
end

% Eq. (18): ai belongs to subset_i or ai = -1/no target.
C = cell(1,Ns);
for k = 1:Ns
    C{k} = [0 subsets{k}];   % 0 means no target
end

if Ns == 1
    actions = C{1}(:);
elseif Ns == 2
    [A1,A2] = ndgrid(C{1}, C{2});
    actions = [A1(:), A2(:)];
else
    grids = cell(1,Ns);
    [grids{:}] = ndgrid(C{:});
    actions = zeros(numel(grids{1}), Ns);
    for k = 1:Ns
        actions(:,k) = grids{k}(:);
    end
end

% Remove all-zero action
actions(all(actions==0,2),:) = [];

% Remove duplicate target in same action
keep = true(size(actions,1),1);
for r = 1:size(actions,1)
    nonzero = actions(r, actions(r,:)>0);
    if numel(unique(nonzero)) < numel(nonzero)
        keep(r) = false;
    end
end
actions = actions(keep,:);

% Sort and keep top nchild by immediate cost
if isempty(actions)
    return;
end

costs = zeros(size(actions,1),1);
for r = 1:size(actions,1)
    costs(r) = action_cost(ssc, actions(r,:), targets, P, cfg);
end

[~,ord] = sort(costs,'ascend');
nkeep = min(cfg.nchild, numel(ord));
actions = actions(ord(1:nkeep),:);

end

% -------------------------------------------------------------------------
function c = action_cost(ssc, act, targets, P, cfg)

c = 0;

for k = 1:cfg.Ns
    j = act(k);

    if j > 0
        dt = transfer_time(ssc(k,:), targets(j,:), P);
        dv = eor_estimate(ssc(k,:), targets(j,:), dt, P, cfg.NmaxPhase);

        if ~isfinite(dv)
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
