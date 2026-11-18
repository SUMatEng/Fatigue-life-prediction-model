function out = fatigueLife(stressRange, defects, mat, opts)
%FLPM.FATIGUELIFE  Cycles to failure from a set of interacting defects.
%
%   out = flpm.fatigueLife(dSigma, defects, mat)
%   out = flpm.fatigueLife(dSigma, defects, mat, flpm.options(...))
%
%   dSigma   applied stress range (MPa) at the load ratio opts.LoadRatio
%   defects  defect set (flpm.defectSet / flpm.samplePopulation):
%            size (sqrt(area), mm), depth (centre below surface, mm),
%            position (along the surface, mm)
%   mat      flpm.material(...)
%
%   Model (after Acta Mater. 237 (2022) 118189 and Mater. Des. 266 (2026)
%   116162; single form, see docs/MODEL.md):
%     da/dN = C1 F dK^(m-p) (dK - dKth)^p        a = sqrt(area)
%     F     = 1 - (1 - U^m)(1 - exp(-Da/l))
%     dKth  = dKth,eff + (dKth,lc - dKth,eff)(1 - exp(-Da/l))
%     dK    = Y dSigma sqrt(pi (A + rp)),  rp = (1/6pi)(Kmax/sigma0)^2
%   with Y and the effective sqrt(area) A from Murakami's half-space rules
%   (flpm.internal.stressIntensity), U from Newman's closure function and
%   l the difference of the El Haddad lengths of dKth,lc and dKth,eff
%   (Y = 0.65). A crack grows once F2*dK >= dKth,eff with
%   F2 = 1 - (1 - U)(1 - exp(-Da/l)); below dKth it is arrested.
%   Neighbouring cracks multiply K by a continuous interaction factor
%   (Kachanov's method, flpm.internal.interactionFactor); cracks whose
%   projections touch link up and are re-characterised after BS 7910
%   (flpm.internal.coalesce). The run ends with
%     failure   Kmax >= KIC, with K not amplified by neighbours
%     run-out   no crack can grow, or N > 1.02 * TargetLife
%
%   Output struct
%     N        cycles to failure (Inf for a run-out)
%     status   'Kc' | 'runout' | 'arrest' | 'noGrowth'
%     killer   original index of the defect with the largest dK at the end
%     nSteps   number of increments
%     final    state at the end (mm)
%     history  (opts.RecordHistory) N, a_max (mm), number of defects

if nargin < 4 || isempty(opts)
    opts = flpm.options();
end
if ~(isscalar(stressRange) && isfinite(stressRange) && stressRange > 0)
    error('flpm:fatigueLife:stress', 'Stress range must be a positive scalar.');
end

S = flpm.internal.prepareState(defects);
c = constants(mat, opts);
s = stressRange;

% initial dK picks the representative of defects that touch from the start
[Y0, A0] = flpm.internal.stressIntensity(S);
S.dK = s * Y0 .* sqrt(pi * A0);

N = 0;
nSteps = 0;
gRef = struct('id', [], 'a', [], 'gam', [], 'tol', []);
hist = zeros(opts.RecordHistory * 1024, 3);
while true
    nSteps = nSteps + 1;
    S = flpm.internal.coalesce(S);
    [Y, A, cls] = flpm.internal.stressIntensity(S);
    % interaction factor, re-evaluated when a crack has grown by more than
    % 2 Eta of its size or of its gap to the nearest neighbour since the
    % last evaluation
    if ~isequal(gRef.id, S.id) || any(abs(S.a - gRef.a) > gRef.tol)
        gRef = struct('id', S.id, 'a', S.a, 'gam', flpm.internal.interactionFactor(S), ...
            'tol', 2 * opts.Eta * min(S.a, sqrt(pi) * nearestGap(S)));
    end
    gam = gRef.gam;
    Y = Y .* gam;                      % interaction amplifies or shields K
    S.dK = s * Y .* sqrt(pi * A);

    % fracture of the specimen: K of each crack without the amplification
    % by its neighbours (amplified K >= KIC on the side facing a neighbour
    % breaks the ligament, i.e. link-up, which growth then reaches)
    if max(S.dK ./ max(gam, 1) .* local(S, c, s)) >= mat.KIC
        status = 'Kc';
        break
    end
    [rate, active] = growthRates(S, Y, A, s, c, mat);
    if nSteps == 1 && ~any(active)
        status = 'noGrowth'; N = Inf;
        break
    end
    if ~any(rate > 0)
        status = 'arrest'; N = Inf;
        break
    end
    if N > opts.MaxLife
        status = 'runout'; N = Inf;
        break
    end

    dN = stepSize(S, rate, cls, c, opts);

    % Heun predictor-corrector
    T = S;
    T.a = S.a + rate * dN;
    T.ext = S.ext + rate * dN;
    [YT, AT] = flpm.internal.stressIntensity(T);
    YT = YT .* gam;
    rateT = growthRates(T, YT, AT, s, c, mat);
    rate = 0.5 * (rate + rateT .* (rate > 0));

    S.a = S.a + rate * dN;
    S.ext = S.ext + rate * dN;
    N = N + dN;
    if opts.RecordHistory
        if nSteps > size(hist, 1), hist = [hist; zeros(size(hist))]; end %#ok<AGROW>
        hist(nSteps, :) = [N, max(S.a) * 1e3, numel(S.a)];
    end
end

[~, k] = max(S.dK);
out.N = N;
out.status = status;
out.killer = S.id(k);
out.nSteps = nSteps;
out.final = struct('size', S.a * 1e3, 'depth', S.z * 1e3, ...
    'position', S.s * 1e3, 'axial', S.y * 1e3, 'extension', S.ext * 1e3, 'id', S.id);
if opts.RecordHistory
    h = hist(1:min(nSteps, size(hist, 1)), :);
    h = h(h(:,1) > 0, :);
    out.history = struct('N', h(:,1), 'aMax', h(:,2), 'nDefects', h(:,3));
end
end

% -------------------------------------------------------------------------
function c = constants(mat, opts)
k = flpm.internal.constants();
c.R = opts.LoadRatio;
c.U = flpm.internal.newmanU(c.R, mat.alpha);
c.alpha = mat.alpha;
c.res = opts.ResidualStress;
a0  = (1/pi) * (mat.dKth_lc  / (k.Ysurface * mat.dSigma_e))^2;   % El Haddad
a0e = (1/pi) * (mat.dKth_eff / (k.Ysurface * mat.dSigma_e))^2;
c.l = a0 - a0e;
c.rpFactor = k.rpFactor;
c.proximity = k.proximity;
end

function [rate, active] = growthRates(S, Y, A, s, c, mat)
% NASGRO / Maierhofer growth rate of every defect (m/cycle).
dK = s * Y .* sqrt(pi * A);
[kmax, U, isOpen] = local(S, c, s);
rp = c.rpFactor * (dK .* kmax / mat.sigma0).^2;
dKg = s * Y .* sqrt(pi * (A + rp));
decay = exp(-S.ext / c.l);
[C1, dKe, dKl, noCl] = flpm.internal.growthParameters(S, mat);
U = U .* ones(size(S.a));
U(noCl) = 1;
dKth = dKe + (dKl - dKe) .* (1 - decay);
F1 = 1 - (1 - U.^mat.m) .* (1 - decay);
F2 = 1 - (1 - U) .* (1 - decay);
active = F2 .* dK >= dKe & isOpen;
rate = zeros(size(S.a));
g = active & dKg > dKth;
rate(g) = C1(g) .* F1(g) .* dKg(g).^(mat.m - mat.p) .* (dKg(g) - dKth(g)).^mat.p;
end

function g = nearestGap(S)
% edge gap of every crack to its nearest neighbour (Inf if alone)
n = numel(S.a);
if n < 2, g = Inf(n, 1); return, end
r = S.a / sqrt(pi);
D = sqrt((S.s - S.s.').^2 + (S.z - S.z.').^2 + (S.y - S.y.').^2) - (r + r.');
D(1:n+1:end) = Inf;
g = max(min(D, [], 2), 0);
end

function [kmax, U, isOpen] = local(S, c, s)
% Kmax/dK, closure ratio and open flag of every defect. Without residual
% stress these are the nominal 1/(1-R) and U(R).
if isempty(c.res)
    kmax = ones(size(S.a)) / (1 - c.R); U = c.U; isOpen = true(size(S.a));
    return
end
smax = s / (1 - c.R);
sr = c.res(S.z * 1e3);
sr = sr(:) .* ones(size(S.a));
isOpen = smax + sr > 0;
kmax = max(smax + sr, 0) / s;
Rl = (smax * c.R + sr) ./ max(smax + sr, eps);
U = flpm.internal.newmanU(min(Rl, 0.99), c.alpha);
end

function dN = stepSize(S, rate, cls, c, opts)
% Largest increment such that no crack advances by more than Eta of its
% governing length (its size; the closure length while the R-curve builds
% up) and that contact between discs and the near-surface and surface
% transitions are resolved.
eta = opts.Eta;
g = rate > 0;
len = S.a;
building = S.ext < 5 * c.l;
len(building) = min(len(building), c.l);
dN = min(eta * len(g) ./ rate(g));

r = S.a / sqrt(pi);
vr = rate / sqrt(pi);
small = 0.01 * eta * r;
n = numel(S.a);
if n > 1
    % distance to contact of the projections (link-up), which closes at
    % most at (dr1 + dr2)/dN; aligned planes only (flpm.internal.coalesce)
    rr = r + r.';
    gap = hypot(S.s - S.s.', S.z - S.z.') - rr;
    aligned = abs(S.y - S.y.') <= rr;
    v = vr + vr.';
    allow = 0.5 * gap + min(small, small.');
    mask = v > 0 & ~eye(n, 'logical') & gap > 0 & aligned;
    if any(mask(:))
        dN = min(dN, min(allow(mask) ./ v(mask)));
    end
end
% internal -> near-surface (r = z/(1 + 1/proximity)) and -> surface (r = z)
for j = find(g & cls < 2).'
    if cls(j) == 0
        rStar = S.z(j) / (1 + 1/c.proximity);
    else
        rStar = S.z(j);
    end
    dN = min(dN, (0.5 * (rStar - r(j)) + small(j)) / vr(j));
end
dN = min(max(dN, opts.MinStep), opts.MaxStep);
end
