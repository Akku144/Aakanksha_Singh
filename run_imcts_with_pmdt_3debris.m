function out = run_imcts_with_pmdt_3debris()
%==========================================================================
% run_imcts_with_pmdt_3debris

%
%==========================================================================

    clc;
    clearvars -except out;

    %% ---------------- USER SETTINGS -------------------------------------
    debrisList = {'H2AF15','ALOS2','GOSAT'};

    % PMDT settings 
    cfg.optimizeFor = 'fuel';             % 'fuel' or 'time'
    cfg.m_wet_servicer = 800;             % kg
    cfg.T = 60/1e3;                       % N (60 mN)
    cfg.Isp = 1300;                       % s
    cfg.target_altitude = 350e3;          % m
    cfg.dvLimit = 1500;                   % m/s
    cfg.a_pchangeLimit = 0.01;
    cfg.dutyRatio = 0.5;
    cfg.Cd = 2.2;
    cfg.Area = 5;                         % m^2
    cfg.waitTimeLimit = 50000;             % days
    cfg.RDV1_days = 45;
    cfg.RDV2_days = 30;
    cfg.launchDate = datetime(2022,03,25,6,37,13);

    cfg.krange = -20:1:20;
    cfg.dragfactor = 1;
    cfg.eclipses = false;
    cfg.drag = false;
    cfg.plots = false;
    cfg.refinement = true;
    cfg.guidance = false;
    cfg.gtype = 1;
    cfg.optimt0 = 0;

    % IMCTS settings
    cfg.mcts_iterations = 6;      % 3 debris -> at most 6 complete sequences
    cfg.c_uct = 1.4;
    cfg.rollouts_per_expansion = 1;
    cfg.verbose = false;

    % fmincon settings
    cfg.fmincon_options = optimoptions('fmincon', ...
        'Display', 'iter', ...
        'Algorithm', 'sqp', ...
        'MaxFunctionEvaluations', 3000, ...
        'MaxIterations', 80, ...
        'StepTolerance', 1e-6, ...
        'ConstraintTolerance', 1e-6, ...
        'OptimalityTolerance', 1e-6);

    %% ---------------- PATH SETUP ----------------------------------------
    root = fileparts(mfilename('fullpath'));
    addpath(root);
    addpath(fullfile(root, 'SGP4routines_NAIF'));
    addpath(fullfile(root, 'CoordinateConversions'));
    addpath(fullfile(root, 'Pertubations'));
    addpath(fullfile(root, 'Data'));
    addpath(fullfile(root, 'PMDT'));

    % Edit if needed
    miceRoot = '/Users/aakankshasingh/Downloads/MICE/mice';
    addpath(fullfile(miceRoot, 'src', 'mice'));
    addpath(fullfile(miceRoot, 'lib'));

    oldDir = pwd;
    cleaner = onCleanup(@() cd(oldDir)); %#ok<NASGU>
    cd(root);

    try
        cspice_furnsh('SGP4routines_NAIF/kernel.txt');
    catch ME
        warning('SPICE kernel load failed: %s', ME.message);
    end

    try
        constants;
    catch ME
        error('Could not load constants.m: %s', ME.message);
    end
    constParam = struct();
    if exist('param', 'var') == 1 && isstruct(param)
        constParam = param;
    end

 
%% ---------------- BUILD PMDT PARAM STRUCT ---------------------------
param = struct();
param.m_wet_servicer = cfg.m_wet_servicer;
param.T = cfg.T;
param.Isp = cfg.Isp;
param.target_altitude = cfg.target_altitude;
param.dvLimit = cfg.dvLimit;
param.t0 = juliandate(cfg.launchDate);
param.a_pchangeLimit = cfg.a_pchangeLimit;
param.dutyRatio = cfg.dutyRatio;
param.Cd = cfg.Cd;
param.Area = cfg.Area;
param.waitTimeLimit = cfg.waitTimeLimit;
param.krange = cfg.krange;
param.dragfactor = cfg.dragfactor;
param.eclipses = cfg.eclipses;
param.drag = cfg.drag;
param.plots = cfg.plots;
param.refinement = cfg.refinement;
param.guidance = cfg.guidance;
param.gtype = cfg.gtype;
param.optimt0 = cfg.optimt0;

if strcmpi(cfg.optimizeFor, 'time')
    param.Topt = true;
else
    param.Topt = false;
end

% ---------- Required PMDT physical constants ----------
% constants.m in this PMDT setup populates `param`, not `C`.
% Keep backward compatibility if a `C` struct exists in another setup.
if exist('C', 'var') == 1 && isstruct(C)
    csrc = C;
elseif exist('constParam', 'var') == 1 && isstruct(constParam)
    csrc = constParam;
else
    csrc = struct();
end

if isfield(csrc,'mu')
    param.mu = csrc.mu;
else
    param.mu = 3.986004418e14;   % fallback
end

if isfield(csrc,'Re')
    param.Re = csrc.Re;
else
    param.Re = 6378137.0;        % fallback
end

if isfield(csrc,'J2')
    param.J2 = csrc.J2;
else
    param.J2 = 1.08262668e-3;    % fallback
end

if isfield(csrc,'g0')
    param.g0 = csrc.g0;
else
    param.g0 = 9.80665;          % fallback
end

if isfield(csrc,'DU'); param.DU = csrc.DU; end
if isfield(csrc,'TU'); param.TU = csrc.TU; end
if isfield(csrc,'VU'); param.VU = csrc.VU; end
if isfield(csrc,'N');  param.N  = csrc.N;  end
if isfield(csrc,'k');  param.k  = csrc.k;  end
if isfield(csrc,'AU'); param.AU = csrc.AU; end
if isfield(csrc,'Rs'); param.Rs = csrc.Rs; end

% Safety fallback for Edelbaum discretization if constants did not provide N
if ~isfield(param, 'N') || isempty(param.N)
    param.N = 100;
end

% ---------- RDV durations ----------
% Use TU if available; otherwise keep in days for now
if isfield(param,'TU')
    param.RDV1 = cfg.RDV1_days * param.TU;
    param.RDV2 = cfg.RDV2_days * param.TU;
else
    warning('TU not found. Using RDV values directly in days. Check constants.m');
    param.RDV1 = cfg.RDV1_days;
    param.RDV2 = cfg.RDV2_days;
end
    %% ---------------- RUN IMCTS -----------------------------------------
    fprintf('\n============================================================\n');
    fprintf('Running IMCTS + PMDT for 3 debris\n');
    fprintf('Debris set: %s, %s, %s\n', debrisList{1}, debrisList{2}, debrisList{3});
    fprintf('Objective : %s\n', upper(cfg.optimizeFor));
    fprintf('Target alt: %.1f km\n', cfg.target_altitude/1e3);
    fprintf('============================================================\n');

    out = imcts_sequence_search_with_pmdt(debrisList, param, cfg);

    % Final detailed PMDT post-processing only for the best sequence
    bestSeq = out.best.sequence;
    bestX = out.best.xopt;
    paramEval = param;
    paramEval.debrisID = bestSeq;
    paramEval.eclipses = true;
    paramEval.drag = true;
    paramEval.plots = false;
    paramEval.refinement = true;
    paramEval.guidance = false;
    paramEval.gtype = 1;

    [resultRaw, RAANs, Total_TOF, Total_dV, fuelconsump] = ...
        PMDTAppPlot([], bestX, bestSeq, paramEval);
    out.best.resultRaw = resultRaw;
    out.best.RAANs = RAANs;
    out.best.Total_TOF = Total_TOF;
    out.best.Total_dV = Total_dV;
    out.best.fuelconsump = fuelconsump;

    %% ---------------- PRINT FINAL RESULT --------------------------------
    fprintf('\n============================================================\n');
    fprintf('FINAL BEST RESULT\n');
    fprintf('============================================================\n');
    fprintf('Best sequence : %s\n', strjoin(out.best.sequence, ' -> '));
    fprintf('Best DV       : %.6f\n', out.best.Total_dV);
    fprintf('Best TOF      : %.6f\n', out.best.Total_TOF);
    fprintf('Best fuel     : %.6f\n', out.best.fuelconsump);
    fprintf('Best xopt     : [%s]\n', num2str(out.best.xopt, ' %.6f'));
    fprintf('============================================================\n');

    save('IMCTS_PMDTOUT_3debris.mat', 'out');
end


%==========================================================================
% IMCTS SEARCH FUNCTION
%==========================================================================
function out = imcts_sequence_search_with_pmdt(debrisList, baseParam, cfg)

    nTargets = numel(debrisList);
    totalUniqueSequences = factorial(nTargets);

    % Node structure
    % node.sequence      -> current partial sequence
    % node.remaining     -> remaining debris
    % node.parent        -> parent index
    % node.children      -> child node indices
    % node.visits        -> visit count
    % node.totalReward   -> accumulated reward
    % node.isTerminal    -> true if sequence complete
    %
    tree = struct( ...
        'sequence', {{}}, ...
        'remaining', {{}}, ...
        'parent', 0, ...
        'children', [], ...
        'visits', 0, ...
        'totalReward', 0, ...
        'isTerminal', false, ...
        'cachedResult', []);

    % Root node
    tree(1).sequence = {};
    tree(1).remaining = debrisList;
    tree(1).parent = 0;
    tree(1).children = [];
    tree(1).visits = 0;
    tree(1).totalReward = 0;
    tree(1).isTerminal = false;
    tree(1).cachedResult = [];

    bestResult = [];
    bestMetric = inf;

    completeSeqMap = containers.Map();

    for iter = 1:cfg.mcts_iterations
        if completeSeqMap.Count >= totalUniqueSequences
            break;
        end

        %---------------- Selection ----------------
        nodeIdx = 1;
        while ~isempty(tree(nodeIdx).children) && ~tree(nodeIdx).isTerminal
            nodeIdx = select_child_uct(tree, nodeIdx, cfg.c_uct);
        end

        %---------------- Expansion ----------------
        if ~tree(nodeIdx).isTerminal
            if isempty(tree(nodeIdx).children)
                tree = expand_node(tree, nodeIdx);
            end

            if ~isempty(tree(nodeIdx).children)
                childList = tree(nodeIdx).children;
                unvisited = childList([tree(childList).visits] == 0);
                if ~isempty(unvisited)
                    nodeIdx = unvisited(1);
                else
                    nodeIdx = childList(randi(numel(childList)));
                end
            end
        end

        %---------------- Rollout / Simulation ----------------
        seq = rollout_to_complete_sequence(tree(nodeIdx).sequence, tree(nodeIdx).remaining);

        seqKey = strjoin(seq, '>');
        if isKey(completeSeqMap, seqKey)
            simResult = completeSeqMap(seqKey);
        else
            simResult = evaluate_sequence_with_pmdt(seq, baseParam, cfg);
            completeSeqMap(seqKey) = simResult;
            if completeSeqMap.Count >= totalUniqueSequences && cfg.verbose
                fprintf('\nAll unique sequences evaluated. Stopping early.\n');
            end
        end

        % Reward = negative cost (because MCTS usually maximizes reward)
        if strcmpi(cfg.optimizeFor, 'time')
            metric = simResult.Total_TOF;
        else
            metric = simResult.Total_dV;
        end
        reward = -metric;

        %---------------- Backpropagation ----------------
        backIdx = nodeIdx;
        while backIdx ~= 0
            tree(backIdx).visits = tree(backIdx).visits + 1;
            tree(backIdx).totalReward = tree(backIdx).totalReward + reward;
            backIdx = tree(backIdx).parent;
        end

        %---------------- Track best ----------------
        if metric < bestMetric
            bestMetric = metric;
            bestResult = simResult;
        end

        if cfg.verbose
            fprintf('\nIter %d / %d', iter, cfg.mcts_iterations);
            fprintf('\n  Rollout sequence : %s', strjoin(seq, ' -> '));
            fprintf('\n  Metric           : %.6f', metric);
            fprintf('\n');
        end

        if completeSeqMap.Count >= totalUniqueSequences
            break;
        end
    end

    % Collect all evaluated sequences
    seqKeys = keys(completeSeqMap);
    allResults = cell(numel(seqKeys),1);
    for i = 1:numel(seqKeys)
        allResults{i} = completeSeqMap(seqKeys{i});
    end

    out.tree = tree;
    out.allResults = allResults;
    out.best = bestResult;
end


%==========================================================================
% EXPAND NODE
%==========================================================================
function tree = expand_node(tree, nodeIdx)

    seq = tree(nodeIdx).sequence;
    rem = tree(nodeIdx).remaining;

    if isempty(rem)
        tree(nodeIdx).isTerminal = true;
        return;
    end

    for i = 1:numel(rem)
        childSeq = [seq, rem(i)];
        childRem = rem;
        childRem(i) = [];

        child.sequence = childSeq;
        child.remaining = childRem;
        child.parent = nodeIdx;
        child.children = [];
        child.visits = 0;
        child.totalReward = 0;
        child.isTerminal = isempty(childRem);
        child.cachedResult = [];

        tree(end+1) = child; %#ok<AGROW>
        childIdx = numel(tree);

        tree(nodeIdx).children(end+1) = childIdx;
    end
end


%==========================================================================
% SELECT CHILD BY UCT
%==========================================================================
function bestChildIdx = select_child_uct(tree, nodeIdx, c)

    children = tree(nodeIdx).children;
    parentVisits = max(tree(nodeIdx).visits, 1);

    uctVals = -inf(size(children));

    for k = 1:numel(children)
        cidx = children(k);
        if tree(cidx).visits == 0
            uctVals(k) = inf;
        else
            avgReward = tree(cidx).totalReward / tree(cidx).visits;
            explore = c * sqrt(log(parentVisits) / tree(cidx).visits);
            uctVals(k) = avgReward + explore;
        end
    end

    [~, ind] = max(uctVals);
    bestChildIdx = children(ind);
end


%==========================================================================
% RANDOM ROLLOUT TO COMPLETE SEQUENCE
%==========================================================================
function seq = rollout_to_complete_sequence(partialSeq, remaining)

    seq = partialSeq;
    rem = remaining;

    while ~isempty(rem)
        idx = randi(numel(rem));
        seq = [seq, rem(idx)]; %#ok<AGROW>
        rem(idx) = [];
    end
end


%==========================================================================
% PMDT EVALUATION FOR A GIVEN COMPLETE SEQUENCE
%==========================================================================
function simResult = evaluate_sequence_with_pmdt(seq, baseParam, cfg)

    param = baseParam;
    param.debrisID = seq;

    % Try PMDT's own initial guess generator first
    try
        [x0_guess, ub, lb] = getInitGuess(param);
    catch
        x0_guess = [];
        ub = [];
        lb = [];
    end

    % Fallback to your app's hardcoded x0 if needed
    if isempty(x0_guess) || ~isnumeric(x0_guess)
        x0 = [7668.55873489256, 1.73887318895310, ...
              7459.17412371716, 1.70295207272273];
    else
        x0 = x0_guess;
    end

    x0 = x0(:).';
    lb = lb(:).';
    ub = ub(:).';

    objfun = @(x) PMDT(x, seq, param);
    nonlcon = @(x) PMDT_limits(x, seq, param);

    [xopt, fval, exitflag, output] = fmincon( ...
        objfun, x0, [], [], [], [], lb, ub, nonlcon, cfg.fmincon_options);

    simResult = struct();
    simResult.sequence = seq;
    simResult.sequenceString = strjoin(seq, ' -> ');
    simResult.x0 = x0;
    simResult.xopt = xopt;
    simResult.fval = fval;
    simResult.exitflag = exitflag;
    simResult.output = output;

    % Use fmincon objective during IMCTS search; do detailed PMDTAppPlot only once for final best.
    if strcmpi(cfg.optimizeFor, 'time')
        simResult.Total_TOF = fval;
        simResult.Total_dV = NaN;
    else
        simResult.Total_dV = fval;
        simResult.Total_TOF = NaN;
    end
    simResult.fuelconsump = NaN;
    simResult.resultRaw = [];
    simResult.RAANs = [];
end
