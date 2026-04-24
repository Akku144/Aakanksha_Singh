function [fuelKmps, timeSec, meta] = evaluate_solution(seq1, seq2, caseData, model, calib, strictComplete)
if nargin < 6
    strictComplete = true;
end

seqs = {seq1(:)', seq2(:)'};
N = numel(caseData.targetIds);
allIds = [seqs{1}, seqs{2}];
uniqueIds = unique(allIds);
isValid = all(ismember(allIds, caseData.targetIds'));
isUnique = numel(uniqueIds) == numel(allIds);
isComplete = numel(allIds) == N && isUnique && isValid;

meta = struct();
meta.totalTargets = N;
meta.assignedCount = numel(uniqueIds);
meta.complete = isComplete;
meta.missingIds = setdiff(caseData.targetIds(:)', uniqueIds, 'stable');
meta.missingCount = numel(meta.missingIds);
meta.fuelRawKmps = nan;
meta.timeRawSec = nan;
meta.fuelScaledKmps = nan;
meta.timeScaledSec = nan;

if (strictComplete && ~isComplete) || (~isValid) || (~isUnique)
    fuelKmps = 1e9;
    timeSec = 1e9;
    meta.fuelRawKmps = fuelKmps;
    meta.timeRawSec = timeSec;
    meta.fuelScaledKmps = fuelKmps;
    meta.timeScaledSec = timeSec;
    return;
end

fuelRaw = 0;
timeRaw = 0;

for s = 1:2
    curOE = caseData.ssc(s, :);
    curT = 0;
    ids = seqs{s};
    for k = 1:numel(ids)
        tgt = caseData.targets(ids(k), :);
        [dv, dt] = transfer_cost(curOE, tgt, curT, model);
        fuelRaw = fuelRaw + dv;
        timeRaw = timeRaw + dt;
        curOE = tgt;
        curT = curT + dt;
    end
end

fuelKmps = fuelRaw * calib.fuelScale;
timeSec = timeRaw * calib.timeScale;
meta.fuelRawKmps = fuelRaw;
meta.timeRawSec = timeRaw;
meta.fuelScaledKmps = fuelKmps;
meta.timeScaledSec = timeSec;
end
