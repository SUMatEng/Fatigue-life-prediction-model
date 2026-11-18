function M = arrestMap(mat, opts, varargin)
%FLPM.ARRESTMAP  Fatigue strength of a single defect over size and depth.
%
%   M = flpm.arrestMap(mat)
%   M = flpm.arrestMap(mat, opts, 'SizeRange', [0.005 2], 'UseParallel', true)
%   s = M.eval(sizeMm, depthMm)          % strength of isolated defects, MPa
%
%   With crack arrest and no interaction, the strength of a set of defects
%   is the smallest strength of its members (weakest link). This map
%   tabulates that single-defect strength, dSigma_w(sqrt(area), z), by
%   flpm.fatigueStrength on a grid of size and relative depth
%   zeta = z/r (r = sqrt(area)/sqrt(pi)), and interpolates it linearly in
%   (log size, zeta, log dSigma_w). The grid is split at Murakami's
%   proximity limit zeta_c = 1 + 1/0.8, where dSigma_w jumps; each side is
%   interpolated separately. Defects deeper than the grid take the value
%   at its deepest row (an internal defect far from the surface).
%
%   Name-value pairs
%     'SizeRange'    [min max] sqrt(area), mm                   [0.005 2]
%     'NumSizes'     grid points in size (log spaced)                  [36]
%     'MaxRelDepth'  deepest zeta tabulated                             [8]
%     'UseParallel'  parfor if a pool is open                       [true]
%
%   Output struct: size (mm), zetaNear, zetaDeep, sNear, sDeep (MPa),
%   eval (function handle), and the inputs.
%
%   See also flpm.weakestLink, flpm.fatigueStrength.

if nargin < 2 || isempty(opts), opts = flpm.options(); end
ip = inputParser;
ip.addParameter('SizeRange', [0.005 2]);
ip.addParameter('NumSizes', 36);
ip.addParameter('MaxRelDepth', 8);
ip.addParameter('UseParallel', true);
ip.parse(varargin{:});
o = ip.Results;

k = flpm.internal.constants();
zc = 1 + 1/k.proximity;
sz = logspace(log10(o.SizeRange(1)), log10(o.SizeRange(2)), o.NumSizes)';
zN = [linspace(0, 1, 11), linspace(1.05, zc, 13)];
zN(end) = zc * (1 - 1e-9);                 % last near-surface point
zD = [zc * (1 + 1e-9), linspace(2.5, 4, 4), linspace(5, o.MaxRelDepth, 4)];
[A1, Z1] = ndgrid(sz, zN);
[A2, Z2] = ndgrid(sz, zD);
a = [A1(:); A2(:)];
z = [Z1(:); Z2(:)] .* a / sqrt(pi);
s = zeros(size(a));
n = numel(a);
if o.UseParallel && canUseParallelPool()
    parfor j = 1:n
        s(j) = one(a(j), z(j), mat, opts);
    end
else
    for j = 1:n
        s(j) = one(a(j), z(j), mat, opts);
    end
end
M.size = sz;
M.zetaNear = zN;
M.zetaDeep = zD;
M.sNear = reshape(s(1:numel(A1)), size(A1));
M.sDeep = reshape(s(numel(A1)+1:end), size(A2));
M.material = mat;
M.options = opts;
M.zetaCritical = zc;
gN = griddedInterpolant({log(sz), zN}, log(M.sNear), 'linear', 'nearest');
gD = griddedInterpolant({log(sz), zD}, log(M.sDeep), 'linear', 'nearest');
M.eval = @(size, depth) evalMap(gN, gD, zc, sz, size, depth);
end

function s = one(a, z, mat, opts)
r = flpm.fatigueStrength(flpm.defectSet(a, z, 0), mat, opts);
s = r.stress;
end

function s = evalMap(gN, gD, zc, sz, a, z)
if any(a(:) < sz(1) * (1 - 1e-9) | a(:) > sz(end) * (1 + 1e-9))
    error('flpm:arrestMap:range', ...
        'Size outside the map (%.4g-%.4g mm); rebuild it with a wider ''SizeRange''.', sz(1), sz(end));
end
zeta = z * sqrt(pi) ./ a;
s = zeros(size(a));
near = zeta < zc;
s(near) = exp(gN(log(a(near)), zeta(near)));
s(~near) = exp(gD(log(a(~near)), zeta(~near)));
end

function tf = canUseParallelPool()
tf = false;
try
    tf = ~isempty(ver('parallel')) && ~isempty(gcp('nocreate'));
catch
end
end
