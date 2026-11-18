function sets = samplePopulation(pops, nSim, varargin)
%FLPM.SAMPLEPOPULATION  Draw Monte Carlo defect sets from population descriptors.
%
%   sets = flpm.samplePopulation(pops, nSim, 'SurfaceLength', L, 'Seed', 42)
%   set  = flpm.samplePopulation(pops, 'simulation', k, 'SurfaceLength', L)
%
%   pops  one flpm.defectPopulation struct or a cell array of them; the
%         sets of all populations are appended.
%   'SurfaceLength'  length L (mm) of free surface represented by one
%         simulation (for a round bar of diameter D, L = pi*D). It sets the
%         number of defects for 'CountPerLength' populations and the range
%         of positions along the surface, i.e. the size of the statistical
%         sample. It does not enter the stress-intensity solution. Required.
%   'AxialLength'    length (mm) along the load axis over which defects are
%         spread (3-D sampling; the gauge length). 0 puts every defect in
%         one plane, as in the published model.                     [0]
%   'Seed'  seed of the random streams                               [0]
%
%   Each set is a struct with column vectors size (sqrt(area), mm), depth
%   (mm), position (mm, uniform on [0, L)), axial (mm, uniform on
%   [0, AxialLength)) and population (index).
%   Surface valleys of depth Rv and shielding factor Y_R become surface-
%   breaking discs centred on the surface whose part below the surface has
%   sqrt(area) = (Y_R/0.65)^2 Rv, so that K = Y_R s sqrt(pi Rv).
%
%   Populations with 'Clustering' place defects by a Thomas process
%   (flpm.clustering): cluster centres uniform along the surface and the
%   axis, at depths from 'ParentDepth' (default the population's 'Depth');
%   Gaussian offsets; positions wrap around the surface length and the
%   axial length (periodic sample), depths are reflected at the surface.
%
%   Simulation k always uses substream k of an mrg32k3a stream seeded with
%   'Seed', so results do not depend on the order of evaluation.

ip = inputParser;
ip.addParameter('Seed', 0);
ip.addParameter('simulation', []);
ip.addParameter('SurfaceLength', [], @(x) isscalar(x) && x > 0);
ip.addParameter('AxialLength', 0, @(x) isscalar(x) && x >= 0);
if ischar(nSim) || isstring(nSim)
    varargin = [{nSim}, varargin];
    nSim = [];
end
ip.parse(varargin{:});
L = ip.Results.SurfaceLength;
if isempty(L)
    error('flpm:sample:length', '''SurfaceLength'' (mm) is required.');
end
if isstruct(pops), pops = num2cell(pops); end
seed = ip.Results.Seed;
Lz = ip.Results.AxialLength;
if ~isempty(ip.Results.simulation)
    sets = drawOne(pops, L, Lz, seed, ip.Results.simulation);
    return
end
sets = cell(nSim, 1);
for k = 1:nSim
    sets{k} = drawOne(pops, L, Lz, seed, k);
end
end

function set = drawOne(pops, L, Lz, seed, k)
rs = RandStream('mrg32k3a', 'Seed', seed);
rs.Substream = k;
Ys = flpm.internal.constants().Ysurface;
parts = cell(numel(pops), 1);
for j = 1:numel(pops)
    p = pops{j};
    if ~isempty(p.countPerLength)
        mu = p.countPerLength * L;
    elseif isfield(p, 'countPerArea') && ~isempty(p.countPerArea)
        if Lz <= 0
            error('flpm:sample:axial', '''CountPerArea'' needs ''AxialLength'' > 0.');
        end
        mu = p.countPerArea * L * Lz;
    else
        mu = p.count;
    end
    if isfield(p, 'clustering') && ~isempty(p.clustering)
        [dep, pos, ax] = clusteredLocations(p, mu, L, Lz, rs);
        n = numel(dep);
        sz = flpm.internal.sampleDescriptor(p.size, n, p.sizeBounds, rs) + p.sizeOffset;
    else
        if strcmp(p.countModel, 'poisson'), n = poissinv(rand(rs), mu); else, n = round(mu); end
        sz = flpm.internal.sampleDescriptor(p.size, n, p.sizeBounds, rs) + p.sizeOffset;
        if strcmp(p.type, 'volumetric')
            dep = flpm.internal.sampleDescriptor(p.depth, n, p.depthBounds, rs);
        else
            dep = zeros(n, 1);
        end
        pos = L * rand(rs, n, 1);
        ax = Lz * rand(rs, n, 1);
    end
    if ~isempty(p.threshold)
        keep = sz > p.threshold;
        sz = sz(keep); dep = dep(keep); pos = pos(keep); ax = ax(keep);
        n = numel(sz);
    end
    if strcmp(p.type, 'surface')
        sz = sqrt(2) * (p.shielding / Ys)^2 * sz;   % half-disc below the surface
    end
    parts{j} = [sz, dep, pos, ax, repmat(j, n, 1)];
end
M = vertcat(parts{:});
if isempty(M), M = zeros(0, 5); end
set = struct('size', M(:,1), 'depth', M(:,2), 'position', M(:,3), ...
    'axial', M(:,4), 'population', M(:,5));
end

function [dep, pos, ax] = clusteredLocations(p, mu, L, Lz, rs)
% Thomas process with mean total count mu and parent intensity mu/m.
c = p.clustering;
m = c.meanClusterSize;
if strcmp(p.countModel, 'poisson')
    np = poissinv(rand(rs), mu / m);
    k = poissinv(rand(rs, np, 1), m);
else                                   % fixed total, shared among the parents
    n = round(mu);
    np = min(n, max(1, round(mu / m)));
    k = accumarray(randi(rs, max(np, 1), n, 1), 1, [np 1]);
end
n = sum(k);
idx = reshape(repelem((1:np)', k), [], 1);   % column also for np = 1
pPos = L * rand(rs, np, 1);
pAx = Lz * rand(rs, np, 1);
pos = mod(pPos(idx) + c.spread(1) * randn(rs, n, 1), L);
if Lz > 0
    ax = mod(pAx(idx) + c.spread(2) * randn(rs, n, 1), Lz);
else
    ax = zeros(n, 1);
end
if strcmp(p.type, 'volumetric')
    d = c.parentDepth;
    if isempty(d), d = p.depth; end
    pDep = flpm.internal.sampleDescriptor(d, np, p.depthBounds, rs);
    lo = p.depthBounds(1); hi = p.depthBounds(2);
    dep = pDep(idx) + c.spread(1) * randn(rs, n, 1);
    dep = lo + abs(dep - lo);                     % reflect at the surface
    dep = min(dep, hi);
else
    dep = zeros(n, 1);
end
end
