function out = run_combined_pmdt_imcts_pmdt_on(varargin)
% PMDT-ON launcher for combined IMCTS + PMDT workflow.
% Keeps run_combined_pmdt_imcts.m free for fast IMCTS-only defaults.
%
% Usage:
%   out = run_combined_pmdt_imcts_pmdt_on();
%   out = run_combined_pmdt_imcts_pmdt_on(40,1500,2025,true,'weighted',5e-9,1.2,100,3,5,120);
%
% Inputs are the same as run_combined_pmdt_imcts, but PMDT refinement is
% forced ON by default here (refineTopK >= 1).

nTargets = 40;
m = 1500;
seed = 2025;
verbose = true;
objective = 'weighted';
timeWeight = 5e-9;
dvmax = 1.2;
infeasiblePenalty = 100;
refineTopK = 3;   % PMDT ON
nSeedRuns = 3;    % smaller default than 5 to limit runtime
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

% Force PMDT refinement ON in this launcher.
refineTopK = max(1, round(refineTopK));

out = run_combined_pmdt_imcts(nTargets, m, seed, verbose, objective, timeWeight, ...
    dvmax, infeasiblePenalty, refineTopK, nSeedRuns, maxLegTofDays);
end
