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
%   Scatter of the threshold between specimens (mat.dKth_cov > 0) is
%   integrated over by Gauss-Hermite quadrature (7 nodes):
%   P_s = sum_k w_k P_s(dSigma | dKth_lc,k). A size-depth copula
%   (defectPopulation 'Copula') enters through joint quadrature weights.
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
%     'Map'            flpm.arrestMap(...) to reuse (a cell, one per
%                      threshold node, if mat.dKth_cov > 0); built if absent
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
    [sz, ws, us] = nodes(p.size, nq(1), p.sizeBounds);
    sz = sz + p.sizeOffset;
    if ~isempty(p.threshold), ws(sz <= p.threshold) = 0; end
    if strcmp(p.type, 'surface')
        sz = sqrt(2) * (p.shielding / Ys)^2 * sz;
        dp = 0; wd = 1;
    else
        [dp, wd, ud] = nodes(p.depth, nq(2), p.depthBounds);
    end
    cop = isfield(p, 'copula') && ~isempty(p.copula) && strcmp(p.type, 'volumetric') ...
        && numel(sz) > 1 && numel(dp) > 1;
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
    if cop                     % joint size-depth weights from the copula density
        [U1, U2] = ndgrid(us, ud);
        c = flpm.internal.copulaPdf(p.copula, [U1(:) U2(:)]);
        c = reshape(c, size(U1));
        c = c / sum(c(:));
        w12 = repmat(c, 1, 1, numel(wp), numel(wy));
        w12(ws <= 0, :, :, :) = 0;
        w = w12 .* w3 .* w4;
    else
        w = w1 .* w2 .* w3 .* w4;
    end
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

% ---- arrest maps (one per threshold node) -------------------------------
allT = vertcat(P{:});
allT = vertcat(allT.T);
[mats, wt] = flpm.internal.thresholdNodes(mat);
Ms = o.Map;
if isempty(Ms)
    rngS = [min(allT.size), max(allT.size)];
    Ms = cellfun(@(m) flpm.arrestMap(m, opts, 'SizeRange', [0.9 1.1] .* rngS), mats, ...
        'UniformOutput', false);
elseif isstruct(Ms)
    Ms = {Ms};
end
if numel(Ms) ~= numel(mats)
    error('flpm:weakestLink:map', ...
        'Give one map per threshold node (%d; flpm.internal.thresholdNodes).', numel(mats));
end
f = [];
if ~isempty(o.StressField)
    f = cellfun(@(q) o.StressField(q.T.depth, q.T.position, q.T.axial), P, 'UniformOutput', false);
end

% ---- weakest link at each node, then mixed over the threshold ---------
R = cell(numel(Ms), 1);
for k = 1:numel(Ms)
    R{k} = core(P, Ms{k}, f);
end
sc = unique(cell2mat(cellfun(@(r) r.sc, R, 'UniformOutput', false)));
stress = o.Stress;
if isempty(stress), stress = sc; end
stress = stress(:);
[Lam, S, Sn] = deal(zeros(size(stress)), zeros(size(stress)), zeros(size(sc)));
K = cell(numel(R), 1);
for k = 1:numel(R)
    [Lk, lSk] = R{k}.at(stress);
    Lam = Lam + wt(k) * Lk;
    S = S + wt(k) * exp(lSk);
    [~, lSn] = R{k}.at(sc);
    Sn = Sn + wt(k) * exp(lSn);
    Kk = R{k}.killers;
    Kk.weight = wt(k) * Kk.weight;
    Kk.thresholdNode = repmat(k, height(Kk), 1);
    K{k} = Kk;
end
W.stress = stress;
W.Lambda = Lam;
W.Pf = 1 - S;
Pfn = 1 - Sn;
W.PfFun = @(s) stepEval(sc, Pfn, s);
W.quantile = @(q) quantileFromStep(sc, Pfn, q);
W.median = W.quantile(0.5);
W.killers = vertcat(K{:});
W.map = Ms;
if numel(Ms) == 1, W.map = Ms{1}; end
W.thresholdWeights = wt;
W.surfaceLength = L;
W.axialLength = Lz;
end

function R = core(P, M, f)
% critical stresses, survival functions and killer probabilities for one map
for i = 1:numel(P)
    T = P{i}.T;
    sw = M.eval(T.size, T.depth);
    fi = ones(height(T), 1);
    if ~isempty(f), fi = f{i}(:) .* ones(height(T), 1); end
    T.sigmaC = sw ./ max(fi, eps);
    T.sigmaC(fi <= 0) = Inf;
    P{i}.T = sortrows(T, 'sigmaC');
end
R.sc = unique(cell2mat(cellfun(@(q) q.T.sigmaC(isfinite(q.T.sigmaC)), P, 'UniformOutput', false)));
R.at = @(s) atStress(P, s);
K = cell(numel(P), 1);
for i = 1:numel(P)
    T = P{i}.T;
    logOther = zeros(height(T), 1);
    for k = 1:numel(P)
        if k == i, continue, end
        [~, lsk] = survival(P{k}, T.sigmaC, true);
        logOther = logOther + lsk;
    end
    % probability that the minimum of population i falls on node j
    % (nodes in sorted order; ties split in that order)
    cw = cumsum(T.w);
    Fb = [0; cw(1:end-1)];
    if P{i}.poisson
        Sfun = @(F) exp(-P{i}.mu * F);
    else
        Sfun = @(F) max(1 - F, 0).^P{i}.mu;
    end
    h = Sfun(Fb) - Sfun(cw);
    h(~isfinite(T.sigmaC)) = 0;
    T.weight = h .* exp(logOther);
    K{i} = T;
end
R.killers = removevars(vertcat(K{:}), 'w');
end

function [Lam, logS] = atStress(P, s)
[Lam, logS] = deal(zeros(numel(s), 1));
for i = 1:numel(P)
    [li, lsi] = survival(P{i}, s);
    Lam = Lam + li;
    logS = logS + lsi;
end
end

% =========================================================================
function [x, w, u] = nodes(d, n, bounds)
% mid-point quantiles of a truncated descriptor; equal weights
if strcmp(d.kind, 'const')
    x = d.value; w = 1; u = 0.5; return
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
