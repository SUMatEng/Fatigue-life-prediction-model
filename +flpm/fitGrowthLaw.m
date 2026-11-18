function G = fitGrowthLaw(dK, dadN, varargin)
%FLPM.FITGROWTHLAW  Calibrate the near-threshold NASGRO law from da/dN data.
%
%   G = flpm.fitGrowthLaw(dK, dadN)                      % long cracks, R = 0.1
%   G = flpm.fitGrowthLaw(dK, dadN, 'LoadRatio', R, 'Material', mat)
%   G = flpm.fitGrowthLaw(dK, dadN, 'Closure', false)   % vacuum / internal
%   mat.sealed = G.sealed;                               % use for internal cracks
%
%   Fits the growth law of the model to measured (dK, da/dN) pairs
%   (MPa sqrt(m), m/cycle) by least squares in log(da/dN):
%     da/dN = C1 F dK^(m-p) (dK - dKth)^p,   dK > dKth
%   with the closure of a long crack, F = U(R)^m and dKth = dKth_lc
%   ('Closure', true, for air data on long cracks), or F = 1 and
%   dKth = dKth_eff ('Closure', false). The latter is for vacuum and
%   internal-crack data, for which the threshold is reported to be nearly
%   independent of R (Vasudevan & Kujawski, FFEMS 45 (2022), doi:
%   10.1111/ffe.13696), i.e. closure is small. The exponents m and p are
%   fitted unless fixed.
%
%   Name-value pairs
%     'LoadRatio'  R of the data (scalar or one per point)         [0.1]
%     'Material'   flpm.material, for alpha and start values  [Ti64-AB]
%     'Closure'    true: long-crack air data; false: no closure   [true]
%     'FixExponents' keep m and p of 'Material'                  [false]
%
%   Output struct
%     C1, m, p, dKth   fitted constants (dKth = dKth_lc or dKth_eff)
%     rmsLog           r.m.s. residual of log10(da/dN)
%     sealed           struct for mat.sealed: with 'Closure', false,
%                      C1, dKth_eff = dKth_lc = dKth and closure = false
%                      (no closure, no R-curve); otherwise C1 and dKth_lc
%     rate             fitted law, @(dK, R) da/dN
%
%   See also flpm.material, flpm.fatigueLife.

ip = inputParser;
ip.FunctionName = 'flpm.fitGrowthLaw';
ip.addParameter('LoadRatio', 0.1);
ip.addParameter('Material', []);
ip.addParameter('Closure', true);
ip.addParameter('FixExponents', false);
ip.parse(varargin{:});
o = ip.Results;
mat = o.Material;
if isempty(mat), mat = flpm.material('Ti64-AB'); end
x = dK(:); y = dadN(:);
ok = isfinite(x) & isfinite(y) & x > 0 & y > 0;
x = x(ok); y = y(ok);
R = o.LoadRatio(:) .* ones(numel(dK), 1);
R = R(ok);
if numel(x) < 4
    error('flpm:fitGrowthLaw:data', 'At least 4 valid (dK, da/dN) points are needed.');
end
if o.Closure
    U = flpm.internal.newmanU(R, mat.alpha);
else
    U = ones(size(x));
end
ly = log(y);
model = @(q, xx, UU) log(exp(q(1))) + log(UU) * q(2) + (q(2) - q(3)) * log(xx) + ...
    q(3) * log(max(xx - q(4), realmin));
    function f = obj(q)
        if q(4) >= min(x) || q(4) <= 0 || q(3) <= 0 || q(2) <= q(3) * 0
            f = 1e10; return
        end
        f = sum((ly - model(q, x, U)).^2);
    end
if o.FixExponents
    full = @(v) [v(1), mat.m, mat.p, v(2)];
    f = @(v) obj(full(v));
    v0 = [log(mat.C1), 0.8 * min(x)];
else
    full = @(v) v;
    f = @obj;
    v0 = [log(mat.C1), mat.m, mat.p, 0.8 * min(x)];
end
best = Inf;
for th = [0.5 0.8 0.95]
    v0(end) = th * min(x);
    [v, fv] = fminsearch(f, v0, optimset('MaxFunEvals', 2e4, 'MaxIter', 2e4, ...
        'TolX', 1e-10, 'TolFun', 1e-12, 'Display', 'off'));
    if fv < best, best = fv; vb = v; end
end
q = full(vb);
G.C1 = exp(q(1)); G.m = q(2); G.p = q(3); G.dKth = q(4);
G.rmsLog = sqrt(best / numel(x)) / log(10);
G.closure = logical(o.Closure);
if o.Closure
    G.sealed = struct('C1', G.C1, 'dKth_lc', G.dKth);
else
    G.sealed = struct('C1', G.C1, 'dKth_eff', G.dKth, 'dKth_lc', G.dKth, 'closure', false);
end
alpha = mat.alpha; cl = o.Closure;
G.rate = @(k, r) exp(model(q, k, uOf(r, alpha, cl, size(k)))) .* (k > q(4));
end

function U = uOf(r, alpha, cl, sz)
if cl
    U = flpm.internal.newmanU(r, alpha) .* ones(sz);
else
    U = ones(sz);
end
end
