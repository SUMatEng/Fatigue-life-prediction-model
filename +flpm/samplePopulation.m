function sets = samplePopulation(pops, geom, nSim, varargin)
%FLPM.SAMPLEPOPULATION  Draw Monte Carlo defect sets from population descriptors.
%
%   sets = flpm.samplePopulation(pops, geom, nSim)
%   sets = flpm.samplePopulation(pops, geom, nSim, 'Seed', 42)
%   set  = flpm.samplePopulation(pops, geom, 'simulation', k, 'Seed', 42)
%
%   pops  one flpm.defectPopulation struct, or a cell array/struct array of
%         them; the sampled sets of all populations are appended (Mater. Des.
%         Sec. 6.3: "sampled populations appended together").
%   Returns a cell array of nSim defect sets (or one set in the second form),
%   each a struct with column vectors
%       size (mm), radius (mm), angle (rad), aspect, yFloor, population
%   ready for flpm.fatigueLife / flpm.fatigueStrength.
%
%   Reproducibility: simulation k always uses substream k of an mrg32k3a
%   stream seeded with 'Seed', independent of the order of evaluation, so
%   results are identical in serial and parallel runs.
%
%   Descriptors are assumed statistically independent (as in the paper,
%   Sec. 6.2); see docs/REVIEW.md for the copula option.

ip = inputParser;
ip.addParameter('Seed', 0);
ip.addParameter('simulation', []);
if ischar(nSim) || isstring(nSim)
    varargin = [{nSim}, varargin];
    nSim = [];
end
ip.parse(varargin{:});
seed = ip.Results.Seed;
pops = normalise(pops);

if ~isempty(ip.Results.simulation)
    sets = drawOne(pops, geom, seed, ip.Results.simulation);
    return
end
sets = cell(nSim, 1);
for k = 1:nSim
    sets{k} = drawOne(pops, geom, seed, k);
end
end

function pops = normalise(pops)
if isstruct(pops)
    pops = num2cell(pops);
end
if ~iscell(pops)
    error('flpm:sample:pops', 'Populations must be a struct or cell array.');
end
end

function set = drawOne(pops, geom, seed, k)
rs = RandStream('mrg32k3a', 'Seed', seed);
rs.Substream = k;
Rb = geom.D / 2;
parts = cell(numel(pops), 1);
for j = 1:numel(pops)
    p = pops{j};
    n = defectCount(p, geom, rs);
    sz = flpm.internal.sampleDescriptor(p.size, n, p.sizeBounds, rs) + p.sizeOffset;
    ang = flpm.internal.sampleDescriptor(p.angle, n, [-pi pi], rs);
    switch p.location
        case 'radius'
            b = p.positionBounds; if isempty(b), b = [0 Rb]; end
            rad = flpm.internal.sampleDescriptor(p.position, n, b, rs);
        case 'depth'
            b = p.positionBounds; if isempty(b), b = [0 Rb]; end
            rad = Rb - flpm.internal.sampleDescriptor(p.position, n, b, rs);
        case 'surface'
            rad = repmat(Rb, n, 1);       % centred on the free surface
    end
    asp = flpm.internal.sampleDescriptor(p.aspect, n, p.aspectBounds, rs);
    if ~isempty(p.threshold)
        keep = sz > p.threshold;
        sz = sz(keep); rad = rad(keep); ang = ang(keep); asp = asp(keep);
        n = numel(sz);
    end
    if strcmp(p.type, 'surface') && ~isempty(p.shielding)
        yf = repmat(p.shielding, n, 1);
    else
        yf = repmat(geom.surfaceY, n, 1);
    end
    parts{j} = [sz, rad, ang, asp, yf, repmat(j, n, 1)];
end
M = vertcat(parts{:});
set = struct('size', M(:,1), 'radius', M(:,2), 'angle', M(:,3), ...
    'aspect', M(:,4), 'yFloor', M(:,5), 'population', M(:,6));
end

function n = defectCount(p, geom, rs)
if ~isempty(p.countPerLength)
    mu = p.countPerLength * pi * geom.D;
else
    mu = p.count;
end
switch p.countModel
    case 'fixed'
        n = round(mu);
    case 'poisson'
        n = poissonFromStream(mu, rs);    % reproducible: uses rs
end
end

function n = poissonFromStream(mu, rs)
% Inverse-CDF Poisson draw from a given RandStream.
u = rand(rs);
n = poissinv(u, mu);
end
