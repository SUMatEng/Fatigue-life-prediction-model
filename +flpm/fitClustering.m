function fit = fitClustering(xyz, varargin)
%FLPM.FITCLUSTERING  Fit a Thomas cluster process to measured defect locations.
%
%   fit = flpm.fitClustering(xy,  'Disc', R)                 % 2-D, disc window
%   fit = flpm.fitClustering(xyz, 'Box', [lo; hi])           % 2-D or 3-D box
%   fit = flpm.fitClustering(xyz, 'Cylinder', [R zmin zmax]) % 3-D, axis = column 3
%   fit = flpm.fitClustering(..., 'Sizes', s, 'Depth', z)    % mark tests
%   pop = flpm.defectPopulation('volumetric', ..., 'Clustering', fit)
%
%   xyz  defect centres (mm), n x 2 or n x 3. In 3-D, column 'AxialDim'
%        (default 3) is the load axis / build direction; the other two
%        span the cross-section. Disc and cylinder windows are centred on
%        the origin in the cross-section.
%
%   Steps
%   1. Ripley's K with translation edge correction,
%        K(r) = |W|/(n(n-1)) sum_{i~=j} 1(d_ij <= r) |W|/|W cap W_{x_i-x_j}|.
%      In 3-D the cylindrical form K(r, t) (in-plane distance <= r, axial
%      separation <= t) is used, so that the in-plane and axial spreads are
%      identified separately.
%   2. Minimum contrast (Diggle): (kappa, sigma) minimise
%      sum (Khat^c - K^c)^2, c = 0.25, for r (and t) up to 'RMax', with
%        2-D   K(r)    = pi r^2 + (1 - exp(-r^2/4s^2)) / kappa
%        3-D   K(r, t) = 2 t pi r^2
%                        + (1 - exp(-r^2/4sp^2)) erf(t/(2 sa)) / kappa
%      mu = lambda/kappa, lambda = n/|W|.
%   3. Validation on the nearest-neighbour distance distribution G(r)
%      (border correction): Monte Carlo deviation test (maximum absolute
%      deviation) against 'NumSimulations' realisations of the fitted
%      Thomas process and of a Poisson process of the same intensity.
%   4. Marks (if 'Sizes' are given): random-labelling permutation test of
%      the correlation between the size of a defect and that of its nearest
%      neighbour; Spearman correlation of size with 'Depth' (if given).
%
%   Name-value pairs
%     'Disc' | 'Box' | 'Cylinder'  window (default: bounding box of xyz)
%     'AxialDim'        axial column in 3-D                              [3]
%     'Isotropic'       3-D: one spread for all directions           [false]
%     'RMax'            largest r (and t) in the contrast; default a
%                       quarter of the smallest window dimension
%     'NumR'            number of r values                              [40]
%     'NumSimulations'  realisations for the G test                     [99]
%     'NumPermutations' permutations for the mark test                 [999]
%     'Sizes', 'Depth'  marks, n x 1
%     'Seed'                                                              [0]
%     'Intensity'       [] homogeneous; 'radial' (disc/cylinder: kernel
%                       estimate of lambda against distance from the axis,
%                       e.g. near-surface porosity); or a handle lambda(x).
%                       The K function is then the inhomogeneous
%                       (intensity-reweighted) K of Baddeley, Moller &
%                       Waagepetersen (Stat. Neerl. 54 (2000) 329), the
%                       Thomas fit is the inhomogeneous Thomas process of
%                       Waagepetersen (Biometrics 63 (2007) 252), and the
%                       G-test simulations are thinned to lambda(x). The
%                       cluster term then measures clustering beyond what
%                       the intensity gradient explains.
%
%   Output struct (lengths in mm)
%     dim, n, volume, lambda    data and intensity (per mm^2 or mm^3)
%     kappa, mu, sigma          parent intensity, mean cluster size,
%                               spread ([plane axial] in 3-D)
%     clustering                flpm.clustering(...) for defectPopulation
%     K                         table r (t), Khat, Kthomas, Kpoisson
%     G                         table r, Gobs, mean and 2.5/97.5% envelopes
%                               of the Thomas and Poisson simulations
%     pThomas, pPoisson         p-values of the G deviation test
%     marks                     nnCorr, pNN, depthRho, pDepth
%
%   See also flpm.clustering, flpm.defectPopulation.

ip = inputParser;
ip.FunctionName = 'flpm.fitClustering';
ip.addParameter('Disc', []);
ip.addParameter('Box', []);
ip.addParameter('Cylinder', []);
ip.addParameter('AxialDim', 3);
ip.addParameter('Isotropic', false);
ip.addParameter('RMax', []);
ip.addParameter('NumR', 40);
ip.addParameter('NumSimulations', 99);
ip.addParameter('NumPermutations', 999);
ip.addParameter('Sizes', []);
ip.addParameter('Depth', []);
ip.addParameter('Seed', 0);
ip.addParameter('Intensity', []);
ip.parse(varargin{:});
o = ip.Results;

[n, dim] = size(xyz);
if n < 10 || ~any(dim == [2 3])
    error('flpm:fitClustering:input', 'Give at least 10 points as an n x 2 or n x 3 array.');
end
if dim == 3                     % put the axial column last
    xyz = xyz(:, [setdiff(1:3, o.AxialDim), o.AxialDim]);
end
W = makeWindow(xyz, o);
lambda = n / W.volume;
lamFun = intensityFunction(o.Intensity, xyz, W);
li = [];
if ~isempty(lamFun)
    li = lamFun(xyz);
    c = sum(1 ./ li) / W.volume;          % so that sum(1/lambda_i) = |W|
    lamFun = @(x) c * lamFun(x);
    li = c * li;
end
rmax = o.RMax;
if isempty(rmax), rmax = 0.25 * W.minExtent; end
r = linspace(rmax / o.NumR, rmax, o.NumR)';
rs = RandStream('mrg32k3a', 'Seed', o.Seed);

% ---- 1. K function ------------------------------------------------------
Khat = kEstimate(xyz, W, r, li);
% ---- 2. minimum contrast -----------------------------------------------
cpow = 0.25;
iso = dim == 2 || o.Isotropic;
Kfun = @(q) kThomas(r, dim, exp(q(1)), exp(q(2)), exp(q(end)));
sMax = W.minExtent;                               % spreads beyond the window
obj = @(q) sum((Khat.^cpow - max(Kfun(q), 0).^cpow).^2, 'all') + ...  % are not identifiable
    1e3 * max(0, q(1) - log(lambda))^2 + ...      % kappa <= lambda (mu >= 1)
    1e3 * sum(max(0, q(2:end) - log(sMax)).^2) + ...
    1e3 * max(0, log(lambda) - q(1) - log(n))^2;  % mu <= n
q0s = log([lambda/5, rmax/4, rmax/4; lambda/20, rmax/10, rmax/10; lambda/2, rmax/2, rmax/2]);
if iso, q0s = q0s(:, 1:2); end
best = Inf;
for k = 1:size(q0s, 1)
    [q, f] = fminsearch(obj, q0s(k, :), optimset('MaxFunEvals', 4000, 'MaxIter', 4000, ...
        'TolX', 1e-8, 'TolFun', 1e-12, 'Display', 'off'));
    if f < best, best = f; qb = q; end
end
kappa = exp(qb(1));
sigma = min(exp(qb(2:end)), sMax);
if iso && dim == 3, sigma = [sigma sigma]; end
mu = max(lambda / kappa, 1);

fit.dim = dim;
fit.n = n;
fit.volume = W.volume;
fit.lambda = lambda;
fit.kappa = kappa;
fit.mu = mu;
fit.sigma = sigma;
fit.contrast = best;
fit.clustering = flpm.clustering('MeanClusterSize', mu, 'Spread', sigma);
if dim == 2
    fit.K = table(r, Khat, Kfun(qb), pi * r.^2, ...
        'VariableNames', {'r', 'Khat', 'Kthomas', 'Kpoisson'});
else  % report the diagonal t = r of the cylindrical K
    Kt = Kfun(qb);
    d = sub2ind(size(Khat), 1:numel(r), 1:numel(r));
    fit.K = table(r, Khat(d)', Kt(d)', 2 * pi * r.^3, ...
        'VariableNames', {'r', 'Khat', 'Kthomas', 'Kpoisson'});
end

% ---- 3. nearest-neighbour validation -----------------------------------
rg = linspace(0, rmax, 60)';
if ~isempty(lamFun)
    lamMax = 1.05 * max(lamFun(uniformIn(W, 20000, rs)));
end
Gobs = gEstimate(xyz, W, rg);
nS = o.NumSimulations;
GT = zeros(numel(rg), nS); GP = GT;
for k = 1:nS
    if isempty(lamFun)
        GT(:, k) = gEstimate(simThomas(kappa, mu, sigma, W, rs), W, rg);
        GP(:, k) = gEstimate(simPoisson(lambda, W, rs), W, rg);
    else            % independent thinning of homogeneous processes at lambda_max
        GT(:, k) = gEstimate(thin(simThomas(kappa, lamMax / kappa, sigma, W, rs), lamFun, lamMax, rs), W, rg);
        GP(:, k) = gEstimate(thin(simPoisson(lamMax, W, rs), lamFun, lamMax, rs), W, rg);
    end
end
[fit.pThomas, mT] = devTest(Gobs, GT);
[fit.pPoisson, mP] = devTest(Gobs, GP);
fit.G = table(rg, Gobs, mT, prctile(GT, 2.5, 2), prctile(GT, 97.5, 2), ...
    mP, prctile(GP, 2.5, 2), prctile(GP, 97.5, 2), 'VariableNames', ...
    {'r', 'Gobs', 'Gthomas', 'GthomasLo', 'GthomasHi', 'Gpoisson', 'GpoissonLo', 'GpoissonHi'});

% ---- 4. marks -----------------------------------------------------------
fit.marks = struct('nnCorr', NaN, 'pNN', NaN, 'depthRho', NaN, 'pDepth', NaN);
if ~isempty(o.Sizes)
    s = o.Sizes(:);
    j = knnsearch(xyz, xyz, 'K', 2);
    j = j(:, 2);
    c0 = corr(s, s(j), 'Type', 'Spearman');
    cp = zeros(o.NumPermutations, 1);
    for k = 1:o.NumPermutations
        sp = s(randperm(rs, n));
        cp(k) = corr(sp, sp(j), 'Type', 'Spearman');
    end
    fit.marks.nnCorr = c0;
    fit.marks.pNN = (1 + sum(abs(cp) >= abs(c0))) / (o.NumPermutations + 1);
    if ~isempty(o.Depth)
        [fit.marks.depthRho, fit.marks.pDepth] = corr(s, o.Depth(:), 'Type', 'Spearman');
    end
end
end

% =========================================================================
function W = makeWindow(xyz, o)
dim = size(xyz, 2);
if ~isempty(o.Disc)
    if dim ~= 2, error('flpm:fitClustering:window', '''Disc'' needs 2-D points.'); end
    R = o.Disc;
    W = struct('type', 'disc', 'R', R, 'volume', pi * R^2, 'minExtent', 2 * R, ...
        'lo', [-R -R], 'hi', [R R]);
elseif ~isempty(o.Cylinder)
    if dim ~= 3, error('flpm:fitClustering:window', '''Cylinder'' needs 3-D points.'); end
    c = o.Cylinder; R = c(1); h = c(3) - c(2);
    W = struct('type', 'cylinder', 'R', R, 'volume', pi * R^2 * h, ...
        'minExtent', min(2 * R, h), 'lo', [-R -R c(2)], 'hi', [R R c(3)]);
else
    b = o.Box;
    if isempty(b), b = [min(xyz); max(xyz)]; end
    W = struct('type', 'box', 'R', NaN, 'volume', prod(diff(b)), ...
        'minExtent', min(diff(b)), 'lo', b(1, :), 'hi', b(2, :));
end
end

function a = overlap(W, dx)
% |W cap W_dx| for pair separations dx (m x dim)
switch W.type
    case 'box'
        a = prod(max((W.hi - W.lo) - abs(dx), 0), 2);
    case 'disc'
        a = discOverlap(W.R, hypot(dx(:, 1), dx(:, 2)));
    case 'cylinder'
        a = discOverlap(W.R, hypot(dx(:, 1), dx(:, 2))) .* ...
            max((W.hi(3) - W.lo(3)) - abs(dx(:, 3)), 0);
end
end

function a = discOverlap(R, d)
d = min(d, 2 * R);
a = 2 * R^2 * acos(d / (2 * R)) - 0.5 * d .* sqrt(4 * R^2 - d.^2);
end

function [I, J, dx] = closePairs(xyz, rlim)
idx = rangesearch(xyz, xyz, rlim);
m = cellfun(@numel, idx);
I = repelem((1:size(xyz, 1))', m);
J = [idx{:}]';
keep = I ~= J;
I = I(keep); J = J(keep);
dx = xyz(J, :) - xyz(I, :);
end

function K = kEstimate(xyz, W, r, li)
[n, dim] = size(xyz);
rmax = r(end);
if dim == 2
    [I, J, dx] = closePairs(xyz, rmax);
    w = W.volume ./ overlap(W, dx);
    if ~isempty(li), w = w ./ (li(I) .* li(J)); end
    d = hypot(dx(:, 1), dx(:, 2));
    [~, b] = histc(d, [0; r]); %#ok<HISTC>
    ok = b > 0;
    K = cumsum(accumarray(b(ok), w(ok), [numel(r) 1]));
else
    [I, J, dx] = closePairs(xyz, sqrt(2) * rmax);
    w = W.volume ./ overlap(W, dx);
    if ~isempty(li), w = w ./ (li(I) .* li(J)); end
    dp = hypot(dx(:, 1), dx(:, 2));
    dz = abs(dx(:, 3));
    [~, bp] = histc(dp, [0; r]); %#ok<HISTC>
    [~, bz] = histc(dz, [0; r]); %#ok<HISTC>
    ok = bp > 0 & bz > 0;
    H = accumarray([bp(ok) bz(ok)], w(ok), [numel(r) numel(r)]);
    K = cumsum(cumsum(H, 1), 2);          % K(r_i, t_j)
end
if isempty(li)
    K = K * W.volume / (n * (n - 1));
else
    K = K / W.volume;                     % intensity-reweighted (inhomogeneous) K
end
end

function f = intensityFunction(spec, xyz, W)
% [] (homogeneous), a handle lambda(x), or 'radial': kernel estimate of
% the intensity as a function of the distance from the axis (disc or
% cylinder windows), lambda(r) = n f(r) / (2 pi r [h])
f = [];
if isempty(spec), return, end
if isa(spec, 'function_handle'), f = spec; return, end
if ~strcmpi(spec, 'radial') || ~any(strcmp(W.type, {'disc', 'cylinder'}))
    error('flpm:fitClustering:intensity', '''Intensity'', ''radial'' needs a disc or cylinder window.');
end
rr = hypot(xyz(:, 1), xyz(:, 2));
n = numel(rr);
[~, ~, bw] = ksdensity(rr, 'Support', [0 W.R * (1 + 1e-9)], 'BoundaryCorrection', 'reflection');
h = 1;
if strcmp(W.type, 'cylinder'), h = W.hi(3) - W.lo(3); end
f = @(x) n * ksdensity(rr, min(hypot(x(:, 1), x(:, 2)), W.R), 'Support', [0 W.R * (1 + 1e-9)], ...
    'BoundaryCorrection', 'reflection', 'Bandwidth', bw) ./ (2 * pi * max(hypot(x(:, 1), x(:, 2)), bw) * h);
end

function x = thin(x, lamFun, lamMax, rs)
if isempty(x), return, end
x = x(rand(rs, size(x, 1), 1) < lamFun(x) / lamMax, :);
end

function K = kThomas(r, dim, kappa, sp, sa)
if dim == 2
    K = pi * r.^2 + (1 - exp(-r.^2 / (4 * sp^2))) / kappa;
else
    t = r.';
    K = 2 * t .* (pi * r.^2) + (1 - exp(-r.^2 / (4 * sp^2))) .* erf(t / (2 * sa)) / kappa;
end
end

function G = gEstimate(xyz, W, rg)
% nearest-neighbour distance distribution, border (reduced-sample) corrected
n = size(xyz, 1);
G = nan(size(rg));
if n < 2, return, end
[~, d] = knnsearch(xyz, xyz, 'K', 2);
d = d(:, 2);
switch W.type
    case 'box'
        b = min(min(xyz - W.lo, W.hi - xyz), [], 2);
    case 'disc'
        b = W.R - hypot(xyz(:, 1), xyz(:, 2));
    case 'cylinder'
        b = min([W.R - hypot(xyz(:, 1), xyz(:, 2)), xyz(:, 3) - W.lo(3), W.hi(3) - xyz(:, 3)], [], 2);
end
for k = 1:numel(rg)
    in = b >= rg(k);
    if any(in), G(k) = mean(d(in) <= rg(k)); end
end
end

function [p, m] = devTest(Gobs, Gsim)
m = mean(Gsim, 2, 'omitnan');
u0 = max(abs(Gobs - m), [], 'omitnan');
us = max(abs(Gsim - m), [], 1, 'omitnan');
p = (1 + sum(us >= u0)) / (numel(us) + 1);
end

function x = simPoisson(lambda, W, rs)
x = uniformIn(W, poissinv(rand(rs), lambda * W.volume), rs);
end

function x = simThomas(kappa, mu, sigma, W, rs)
dim = numel(W.lo);
s = repmat(sigma(1), 1, dim);
if dim == 3, s(3) = sigma(end); end
s = min(s, W.minExtent);
lo = W.lo - 4 * s; hi = W.hi + 4 * s;
np = poissinv(rand(rs), kappa * prod(hi - lo));
P = lo + (hi - lo) .* rand(rs, np, dim);
k = poissinv(rand(rs, np, 1), mu);
x = repelem(P, k, 1) + s .* randn(rs, sum(k), dim);
x = x(inWindow(x, W), :);
end

function x = uniformIn(W, n, rs)
dim = numel(W.lo);
x = zeros(0, dim);
while size(x, 1) < n
    y = W.lo + (W.hi - W.lo) .* rand(rs, 2 * n + 10, dim);
    x = [x; y(inWindow(y, W), :)]; %#ok<AGROW>
end
x = x(1:n, :);
end

function tf = inWindow(x, W)
switch W.type
    case 'box'
        tf = all(x >= W.lo & x <= W.hi, 2);
    case 'disc'
        tf = hypot(x(:, 1), x(:, 2)) <= W.R;
    case 'cylinder'
        tf = hypot(x(:, 1), x(:, 2)) <= W.R & x(:, 3) >= W.lo(3) & x(:, 3) <= W.hi(3);
end
end
