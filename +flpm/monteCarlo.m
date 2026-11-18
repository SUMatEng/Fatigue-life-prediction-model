function [T, sets] = monteCarlo(pops, mat, geom, varargin)
%FLPM.MONTECARLO  Fatigue-strength distribution of a defect population.
%
%   T = flpm.monteCarlo(pops, mat, geom)
%   [T, sets] = flpm.monteCarlo(pops, mat, geom, 'NumSimulations', 400, ...
%                   'Seed', 1, 'Options', flpm.modelOptions(), ...
%                   'UseParallel', true)
%
%   For each simulation a defect set is drawn from the population
%   descriptors (flpm.samplePopulation) and its fatigue strength at the
%   target life is solved (flpm.fatigueStrength). This is the stochastic
%   framework of Macallister, Deillon & Becker, Mater. Des. 266 (2026).
%
%   Inputs
%     pops   flpm.defectPopulation struct or cell array of them
%     mat    flpm.material(...)       geom  flpm.geometry(...)
%   Name-value pairs
%     'NumSimulations'  number of Monte Carlo draws               [400]
%     'Seed'            seed of the per-simulation substreams      [0]
%     'Options'         flpm.modelOptions(...)       ['paper' preset]
%     'UseParallel'     use parfor if a pool is available       [false]
%     'DefectSets'      cell array of defect sets to evaluate instead of
%                       sampling (e.g. legacy sim_XCT_* rows); pops may
%                       then be []
%
%   Output table T, one row per simulation
%     stress        fatigue strength dSigma_w (MPa)
%     converged     solver flag
%     N             cycles to failure at 'stress'
%     status        failure mode at 'stress'
%     nDefects      number of defects simulated
%     killer        index of the killer defect in the simulated set
%     killerPopulation, killerSize (mm), killerRadius (mm), killerDepth
%                   (mm, centre to surface), killerAspect, killerAngle,
%     killerSpacing nearest-neighbour edge spacing of the killer (mm)
%     nEval, runtime(s)
%   sets: the simulated defect sets (cell array).
%
%   An empty defect set (possible with an EVT threshold) is assigned the
%   defect-free fatigue strength mat.dSigma_e and flagged status 'noDefects'.

ip = inputParser;
ip.addParameter('NumSimulations', 400, @(x) isscalar(x) && x >= 1);
ip.addParameter('Seed', 0);
ip.addParameter('Options', []);
ip.addParameter('UseParallel', false);
ip.addParameter('DefectSets', []);
ip.parse(varargin{:});
r = ip.Results;
opts = r.Options;
if isempty(opts), opts = flpm.modelOptions(); end

if ~isempty(r.DefectSets)
    sets = r.DefectSets(:);
    nSim = numel(sets);
else
    nSim = r.NumSimulations;
    sets = cell(nSim, 1);
end
sampleInside = isempty(r.DefectSets);
seed = r.Seed;

stress = nan(nSim, 1); conv = false(nSim, 1); N = nan(nSim, 1);
status = strings(nSim, 1); nDef = zeros(nSim, 1); killer = nan(nSim, 1);
kPop = nan(nSim, 1); kSize = nan(nSim, 1); kRad = nan(nSim, 1);
kAsp = nan(nSim, 1); kAng = nan(nSim, 1); kSpc = nan(nSim, 1);
nEval = zeros(nSim, 1); rt = zeros(nSim, 1);
Rb = geom.D / 2;

useParallel = r.UseParallel && canUseParallelPool();
if useParallel
    parfor k = 1:nSim
        [sets{k}, row] = runOne(k, sets{k}, sampleInside, pops, mat, geom, opts, seed);
        [stress(k), conv(k), N(k), status(k), nDef(k), killer(k), kPop(k), ...
            kSize(k), kRad(k), kAsp(k), kAng(k), kSpc(k), nEval(k), rt(k)] = row{:};
    end
else
    for k = 1:nSim
        [sets{k}, row] = runOne(k, sets{k}, sampleInside, pops, mat, geom, opts, seed);
        [stress(k), conv(k), N(k), status(k), nDef(k), killer(k), kPop(k), ...
            kSize(k), kRad(k), kAsp(k), kAng(k), kSpc(k), nEval(k), rt(k)] = row{:};
    end
end

T = table((1:nSim)', stress, conv, N, status, nDef, killer, kPop, kSize, ...
    kRad, Rb - kRad, kAsp, kAng, kSpc, nEval, rt, ...
    'VariableNames', {'simulation', 'stress', 'converged', 'N', 'status', ...
    'nDefects', 'killer', 'killerPopulation', 'killerSize', 'killerRadius', ...
    'killerDepth', 'killerAspect', 'killerAngle', 'killerSpacing', ...
    'nEval', 'runtime'});
T.Properties.VariableUnits = {'', 'MPa', '', 'cycles', '', '', '', '', ...
    'mm', 'mm', 'mm', '', 'rad', 'mm', '', 's'};
end

% -------------------------------------------------------------------------
function [set, row] = runOne(k, set, sampleInside, pops, mat, geom, opts, seed)
t0 = tic;
if sampleInside
    set = flpm.samplePopulation(pops, geom, 'simulation', k, 'Seed', seed);
end
n = numel(set.size);
if n == 0
    row = {mat.dSigma_e, true, Inf, "noDefects", 0, NaN, NaN, NaN, NaN, ...
        NaN, NaN, NaN, 0, toc(t0)};
    return
end
res = flpm.fatigueStrength(set, mat, geom, opts);
kk = res.killer;
pop = NaN;
if isfield(set, 'population'), pop = set.population(kk); end
row = {res.stress, res.converged, res.N, string(res.status), n, kk, pop, ...
    set.size(kk), set.radius(kk), set.aspect(kk), set.angle(kk), ...
    nnSpacing(set, kk), res.nEval, toc(t0)};
end

function s = nnSpacing(set, k)
% Edge-to-edge in-plane distance from defect k to its nearest neighbour,
% treating each defect as a circle of radius sqrt(area) (as in
% KNN_spacing_cartesian.m); 0 if touching.
if numel(set.size) < 2
    s = NaN; return
end
[x, y] = pol2cart(set.angle, set.radius);
d = hypot(x - x(k), y - y(k)) - set.size - set.size(k);
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
