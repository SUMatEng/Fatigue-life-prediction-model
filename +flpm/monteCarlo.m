function [T, sets] = monteCarlo(pops, mat, varargin)
%FLPM.MONTECARLO  Fatigue-strength distribution of a defect population.
%
%   T = flpm.monteCarlo(pops, mat, 'SurfaceLength', L)
%   [T, sets] = flpm.monteCarlo(pops, mat, 'SurfaceLength', L, ...
%                   'NumSimulations', 400, 'Seed', 1, ...
%                   'Options', flpm.options(), 'UseParallel', true)
%
%   For each simulation a defect set is drawn from the population
%   descriptors (flpm.samplePopulation) and its fatigue strength at the
%   target life is solved (flpm.fatigueStrength).
%
%   Name-value pairs
%     'SurfaceLength'   free-surface length per simulation, mm (required
%                       unless 'DefectSets' is given); pi*D for a round bar
%     'AxialLength'     gauge length over which defects are spread, mm [0]
%     'NumSimulations'  number of draws                              [400]
%     'Seed'            seed of the per-simulation substreams          [0]
%     'Options'         flpm.options(...)
%     'UseParallel'     parfor if a pool is open                    [false]
%     'DefectSets'      cell array of defect sets to evaluate instead of
%                       sampling; pops may then be []
%
%   Output table T, one row per simulation
%     stress (MPa), converged, stepInSN, N, status, nDefects, killer,
%     killerPopulation, killerSize (sqrt(area), mm), killerDepth (mm),
%     killerSpacing (edge-to-edge distance to the nearest defect, mm),
%     nEval, runtime (s)
%   An empty defect set is assigned mat.dSigma_e, status 'noDefects'.

ip = inputParser;
ip.addParameter('SurfaceLength', []);
ip.addParameter('AxialLength', 0);
ip.addParameter('NumSimulations', 400, @(x) isscalar(x) && x >= 1);
ip.addParameter('Seed', 0);
ip.addParameter('Options', []);
ip.addParameter('UseParallel', false);
ip.addParameter('DefectSets', []);
ip.parse(varargin{:});
r = ip.Results;
opts = r.Options;
if isempty(opts), opts = flpm.options(); end
sampleInside = isempty(r.DefectSets);
if sampleInside
    if isempty(r.SurfaceLength)
        error('flpm:monteCarlo:length', '''SurfaceLength'' (mm) is required for sampling.');
    end
    nSim = r.NumSimulations;
    sets = cell(nSim, 1);
else
    sets = r.DefectSets(:);
    nSim = numel(sets);
end
L = r.SurfaceLength;
Lz = r.AxialLength;
seed = r.Seed;

rows = cell(nSim, 1);
if r.UseParallel && canUseParallelPool()
    parfor k = 1:nSim
        [sets{k}, rows{k}] = runOne(k, sets{k}, sampleInside, pops, mat, opts, L, Lz, seed);
    end
else
    for k = 1:nSim
        [sets{k}, rows{k}] = runOne(k, sets{k}, sampleInside, pops, mat, opts, L, Lz, seed);
    end
end
T = struct2table(vertcat(rows{:}));
T = addvars(T, (1:nSim)', 'Before', 1, 'NewVariableNames', 'simulation');
T.Properties.VariableUnits = {'', 'MPa', '', '', 'cycles', '', '', '', '', ...
    'mm', 'mm', 'mm', '', 's'};
end

% -------------------------------------------------------------------------
function [set, row] = runOne(k, set, sampleInside, pops, mat, opts, L, Lz, seed)
t0 = tic;
if sampleInside
    set = flpm.samplePopulation(pops, 'simulation', k, 'SurfaceLength', L, ...
        'AxialLength', Lz, 'Seed', seed);
end
n = numel(set.size);
row = struct('stress', mat.dSigma_e, 'converged', true, 'stepInSN', false, ...
    'N', Inf, 'status', "noDefects", 'nDefects', n, 'killer', NaN, ...
    'killerPopulation', NaN, 'killerSize', NaN, 'killerDepth', NaN, ...
    'killerSpacing', NaN, 'nEval', 0, 'runtime', 0);
if n > 0
    res = flpm.fatigueStrength(set, mat, opts);
    kk = res.killer;
    row.stress = res.stress; row.converged = res.converged;
    row.stepInSN = res.stepInSN; row.N = res.N; row.status = string(res.status);
    row.killer = kk;
    if isfield(set, 'population'), row.killerPopulation = set.population(kk); end
    row.killerSize = set.size(kk); row.killerDepth = set.depth(kk);
    row.killerSpacing = nnSpacing(set, kk); row.nEval = res.nEval;
end
row.runtime = toc(t0);
end

function s = nnSpacing(set, k)
% Edge-to-edge distance between equivalent discs (radius sqrt(area)/sqrt(pi))
% from defect k to its nearest neighbour; 0 if touching.
if numel(set.size) < 2
    s = NaN; return
end
r = set.size / sqrt(pi);
y = zeros(size(r));
if isfield(set, 'axial') && ~isempty(set.axial), y = set.axial; end
d = sqrt((set.position - set.position(k)).^2 + (set.depth - set.depth(k)).^2 ...
    + (y - y(k)).^2) - r - r(k);
d(k) = Inf;
s = max(min(d), 0);
end

function tf = canUseParallelPool()
tf = false;
try
    tf = ~isempty(ver('parallel')) && ~isempty(gcp('nocreate'));
catch
end
end
