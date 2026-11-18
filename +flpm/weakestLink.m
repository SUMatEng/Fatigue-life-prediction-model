function W = weakestLink(pops, mat, varargin)
%FLPM.WEAKESTLINK  Fatigue-strength distribution in closed form (weakest link).
%
%   W = flpm.weakestLink(pops, mat, 'SurfaceLength', L)
%   W = flpm.weakestLink(pops, mat, 'SurfaceLength', L, 'AxialLength', Lz, ...
%                        'StressField', @(z, s, y) f, 'Map', M)
%   W.median, W.quantile(0.1), W.Pf(stress), W.killers
%
%   With crack arrest (fatigue limit), a defect fails the part at stress
%   range dSigma when dSigma f(x) >= dSigma_w(a, z), where dSigma_w is
%   the strength of the defect alone (flpm.arrestMap) and f(x) the local
%   stress relative to the nominal one. Without interaction, a set fails
%   at the smallest of these critical stresses. For a Poisson number of
%   defects, with marks (size, depth) and positions independent,
%
%     P_s(dSigma) = exp(-Lambda(dSigma)),
%     Lambda(dSigma) = lambda E[ 1(dSigma_w(a, z) / f(z, s, y) <= dSigma) ],
%
%   i.e. Lambda is the expected number of critical defects. For a fixed
%   count n, P_s = (1 - Lambda/n)^n. Populations combine as a product of
%   P_s. The expectation is evaluated by deterministic quadrature over the
%   size and depth descriptors (truncated inverse CDF at mid-point
%   quantiles) and over position and axial location when a stress field
%   is given. Thresholds and roughness valleys are handled as in
%   flpm.samplePopulation.
%
%   Interaction (BS 7910) and clustering are outside this closed form;
%   flpm.monteCarlo includes both. The closed form is exact in the limit
%   of no interaction, and its difference from flpm.monteCarlo is the
%   effect of interaction.
%
%   Name-value pairs
%     'SurfaceLength'  mm, as for flpm.samplePopulation          (required)
%     'AxialLength'    mm                                               [0]
%     'StressField'    f(depth, position, axial) in mm, the local stress
%                      range relative to the nominal one, evaluated at the
%                      defect centre                                  [1]
%     'Map'            flpm.arrestMap(...) to reuse; built if absent
%     'Options'        flpm.options(...)
%     'NumQuad'        quadrature points [size depth position axial]
%                                                       [200 80 16 16]
%     'Stress'         stresses (MPa) at which to report Pf; default the
%                      critical stresses of the quadrature nodes
%
%   Output struct
%     stress, Lambda, Pf           on the reported stresses
%     Pf(s), quantile(p)           function handles
%     median                       median strength, MPa
%     killers                      table of the quadrature nodes with
%                                  population, size, depth, position,
%                                  axial, sigmaC (MPa) and weight = P(this
%                                  node is the killer); sum(weight) ~ 1
%     map                          the arrest map used
%
%   See also flpm.arrestMap, flpm.monteCarlo, flpm.killerDistribution.

ip = inputParser;
ip.FunctionName = 'flpm.weakestLink';
ip.addParameter('SurfaceLength', []);
ip.addParameter('AxialLength', 0);
ip.addParameter('StressField', []);
ip.addParameter('Map', []);
ip.addParameter('Options', []);
ip.addParameter('NumQuad', [200 80 16 16]);
ip.addParameter('Stress', []);
ip.parse(varargin{:});
o = ip.Results;
if isempty(o.SurfaceLength)
    error('flpm:weakestLink:length', '''SurfaceLength'' (mm) is required.');
end
opts = o.Options;
if isempty(opts), opts = flpm.options(); end
if isstruct(pops), pops = num2cell(pops); end
L = o.SurfaceLength; Lz = o.AxialLength;
nq = o.NumQuad;
Ys = flpm.internal.constants().Ysurface;

% ---- quadrature nodes ----------------------------------------------------
P = cell(numel(pops), 1);
for i = 1:numel(pops)
    p = pops{i};
    if isfield(p, 'clustering') && ~isempty(p.clustering)
        warning('flpm:weakestLink:clustering', ...
            'Population %d is clustered; the closed form treats it as Poisson.', i);
    end
    [sz, ws] = nodes(p.size, nq(1), p.sizeBounds);
    sz = sz + p.sizeOffset;
    if ~isempty(p.threshold), ws(sz <= p.threshold) = 0; end
    if strcmp(p.type, 'surface')
        sz = sqrt(2) * (p.shielding / Ys)^2 * sz;
        dp = 0; wd = 1;
    else
        [dp, wd] = nodes(p.depth, nq(2), p.depthBounds);
    end
    if isempty(o.StressField)
        ps = 0; wp = 1; ya = 0; wy = 1;
    else
        ps = ((1:nq(3))' - 0.5) / nq(3) * L; wp = ones(nq(3), 1) / nq(3);
        if Lz > 0
            ya = ((1:nq(4))' - 0.5) / nq(4) * Lz; wy = ones(nq(4), 1) / nq(4);
        else
            ya = 0; wy = 1;
        end
    end
    [S1, D1, P1, Y1] = ndgrid(sz, dp, ps, ya);
    [w1, w2, w3, w4] = ndgrid(ws, wd, wp, wy);
    w = w1 .* w2 .* w3 .* w4;
    keep = w(:) > 0;
    T = table(repmat(i, nnz(keep), 1), S1(keep), D1(keep), P1(keep), Y1(keep), w(keep), ...
        'VariableNames', {'population', 'size', 'depth', 'position', 'axial', 'w'});
    % count
    if ~isempty(p.countPerLength)
        mu = p.countPerLength * L;
    elseif isfield(p, 'countPerArea') && ~isempty(p.countPerArea)
        if Lz <= 0
            error('flpm:weakestLink:axial', '''CountPerArea'' needs ''AxialLength'' > 0.');
        end
        mu = p.countPerArea * L * Lz;
    else
        mu = p.count;
    end
    poisson = strcmp(p.countModel, 'poisson');
    if ~poisson, mu = round(mu); end
    P{i} = struct('T', T, 'mu', mu, 'poisson', poisson);
end

% ---- arrest map and critical stresses ----------------------------------
allT = vertcat(P{:});
allT = vertcat(allT.T);
M = o.Map;
if isempty(M)
    rng = [min(allT.size), max(allT.size)];
    M = flpm.arrestMap(mat, opts, 'SizeRange', [0.9 1.1] .* rng);
end
for i = 1:numel(P)
    T = P{i}.T;
    sw = M.eval(T.size, T.depth);
    f = ones(height(T), 1);
    if ~isempty(o.StressField)
        f = o.StressField(T.depth, T.position, T.axial);
        f = f(:) .* ones(height(T), 1);
    end
    T.sigmaC = sw ./ max(f, eps);
    T.sigmaC(f <= 0) = Inf;
    P{i}.T = sortrows(T, 'sigmaC');
end

% ---- survival functions ---------------------------------------------
sc = unique(cell2mat(cellfun(@(q) q.T.sigmaC(isfinite(q.T.sigmaC)), P, 'UniformOutput', false)));
stress = o.Stress;
if isempty(stress), stress = sc; end
stress = stress(:);
[Lam, logS] = deal(zeros(size(stress)));
for i = 1:numel(P)
    [li, lsi] = survival(P{i}, stress);
    Lam = Lam + li;
    logS = logS + lsi;
end
W.stress = stress;
W.Lambda = Lam;
W.Pf = 1 - exp(logS);
% exact step functions on the node stresses
[~, logSn] = deal(zeros(size(sc)));
for i = 1:numel(P)
    [~, lsi] = survival(P{i}, sc);
    logSn = logSn + lsi;
end
Pfn = 1 - exp(logSn);
W.PfFun = @(s) stepEval(sc, Pfn, s);
W.quantile = @(q) quantileFromStep(sc, Pfn, q);
W.median = W.quantile(0.5);

% ---- killer probabilities --------------------------------------------
K = cell(numel(P), 1);
for i = 1:numel(P)
    T = P{i}.T;
    logOther = zeros(height(T), 1);
    for k = 1:numel(P)
        if k == i, continue, end
        [~, lsk] = survival(P{k}, T.sigmaC, true);
        logOther = logOther + lsk;
    end
    [~, ~, Fi] = survival(P{i}, T.sigmaC, true);     % mass strictly below
    if P{i}.poisson
        h = P{i}.mu * T.w .* exp(-P{i}.mu * Fi);
    else
        h = P{i}.mu * T.w .* max(1 - Fi, 0).^(P{i}.mu - 1);
    end
    T.weight = h .* exp(logOther);
    K{i} = T;
end
W.killers = removevars(vertcat(K{:}), 'w');
W.map = M;
W.surfaceLength = L;
W.axialLength = Lz;
end

% =========================================================================
function [x, w] = nodes(d, n, bounds)
% mid-point quantiles of a truncated descriptor; equal weights
if strcmp(d.kind, 'const')
    x = d.value; w = 1; return
end
lo = bounds(1); hi = bounds(2);
Flo = 0; Fhi = 1;
if isfinite(lo), Flo = cdf(d.value, lo); end
if isfinite(hi), Fhi = cdf(d.value, hi); end
u = ((1:n)' - 0.5) / n;
x = icdf(d.value, Flo + u * (Fhi - Flo));
x = min(max(x, lo), hi);
w = ones(n, 1) / n;
end

function [Lam, logS, F] = survival(Pi, s, strict)
% Lambda, log P_s of one population at stresses s; F = node mass with
% sigmaC <= s (strictly below if 'strict')
if nargin < 3, strict = false; end
T = Pi.T;
cw = [0; cumsum(T.w)];
c = T.sigmaC(isfinite(T.sigmaC));            % sorted ascending
s = s(:);
idx = zeros(size(s));
if ~isempty(c)
    if strict          % number of nodes with sigmaC < s
        [cu, ia] = unique(c, 'first');
        if numel(cu) == 1
            j = ia + 0 * s; j(s > cu) = numel(c) + 1;
        else
            j = interp1(cu, ia, s, 'next');
            j(s > cu(end)) = numel(c) + 1;
        end
        j(s <= cu(1)) = 1;
        idx = j - 1;
    else               % number of nodes with sigmaC <= s
        [cu, ia] = unique(c, 'last');
        if numel(cu) == 1
            idx = ia * (s >= cu);
        else
            idx = interp1(cu, ia, s, 'previous');
            idx(s >= cu(end)) = numel(c);
        end
        idx(s < cu(1)) = 0;
    end
end
F = cw(idx + 1);
Lam = Pi.mu * F;
if Pi.poisson
    logS = -Lam;
else
    logS = Pi.mu * log(max(1 - F, 0));
end
end

function y = stepEval(x, v, s)
y = zeros(size(s));
for k = 1:numel(s)
    j = find(x <= s(k), 1, 'last');
    if ~isempty(j), y(k) = v(j); end
end
end

function q = quantileFromStep(x, v, p)
q = nan(size(p));
for k = 1:numel(p)
    j = find(v >= p(k), 1, 'first');
    if ~isempty(j), q(k) = x(j); end
end
end
