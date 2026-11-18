function out = fatigueLife(stressRange, defects, mat, geom, opts)
%FLPM.FATIGUELIFE  Cycles to failure of a bar containing multiple defects.
%
%   out = flpm.fatigueLife(dSigma, defects, mat, geom, opts)
%
%   dSigma   applied stress range (MPa) at load ratio geom.R
%   defects  defect set (see flpm.samplePopulation / flpm.defectSet):
%            size (sqrt(area), mm), radius (mm), angle (rad), aspect (a/c)
%   mat      flpm.material(...)       geom  flpm.geometry(...)
%   opts     flpm.modelOptions(...)   (default: 'paper' preset)
%
%   Model (Acta Mater. 237 (2022) 118189; Mater. Des. 266 (2026) 116162):
%     da/dN  = C1 F dK^(m-p) (dK - dKth)^p                       (Eq. 1)
%     F      = 1 - (1 - U^m)(1 - exp(-Da/l))                     (Eq. 2)
%     dKth   = dKth,eff + (dKth,lc - dKth,eff)(1 - exp(-Da/l))   (Eq. 3)
%     dK     = Y sigma_net sqrt(pi (a + rp)),  a = sqrt(area)    (Eq. 4/8/9)
%   with U from Newman's closure function, l the difference of El Haddad
%   lengths for dKth,lc and dKth,eff, rp the cyclic plastic-zone size, Y from
%   flpm.internal.geometryFactors, and sigma_net the growing net-section
%   stress of the largest defect. Defects grow independently, coalesce on
%   contact (flpm.internal.coalesce) and the run ends on
%     Kmax >= KIC, sigma_net,max >= flow stress, a > D       (failure)
%     no defect can grow / N > opts.MaxLife                   (run-out)
%
%   Output struct
%     N         cycles to failure (Inf for a run-out)
%     status    'Kc' | 'netSection' | 'through' | 'runout' | 'noGrowth'
%     killer    original index of the defect with the largest dK at the
%               end of the run (the 'killer' defect, Mater. Des. Sec. 7.2)
%     nSteps    number of increments
%     final     state at the end (sizes in mm)
%     history   (opts.RecordHistory) N, a_max (mm), sigma_net, n defects

if nargin < 5 || isempty(opts)
    opts = flpm.modelOptions();
end
if ~(isscalar(stressRange) && isfinite(stressRange) && stressRange > 0)
    error('flpm:fatigueLife:stress', 'Stress range must be a positive scalar.');
end

S = flpm.internal.prepareState(defects, geom);
c = constants(mat, geom, opts);
s = stressRange;
D = c.D;

% Initial dK, used only to pick the representative of a group that is
% already in contact before the first increment (as in the published code).
[Y0, ~] = flpm.internal.geometryFactors(S, D, opts, c.rScale);
S.dK = c.U * s * Y0 .* sqrt(pi * S.a);

N = 0;
nSteps = 0;
status = '';
sigG = [];
fixed = strcmpi(opts.Stepping, 'fixed');
hist = initHistory(opts.RecordHistory);

while true
    nSteps = nSteps + 1;

    % 1. coalescence, geometry factors, nominal dK
    S = flpm.internal.coalesce(S, c.rScale);
    [Y, cls] = flpm.internal.geometryFactors(S, D, opts, c.rScale);
    dK = s * Y .* sqrt(pi * S.a);
    S.dK = dK;
    [amax, ind] = max(S.a);

    if fixed
        N = N + opts.FixedStep;                 % code counts before growth
    end
    if amax > D
        status = 'through';
        break
    end

    % 2. driving force and thresholds
    if opts.NetSectionInDeltaK
        if isempty(sigG)
            sigG = netStress(amax, s, c);       % initialise (code: undefined)
        elseif amax < D/2 && max(c.U * dK) < opts.NetSectionKCap
            sigG = netStress(amax, s, c);
        end
    else
        sigG = s;
    end
    [rate, active] = growthRates(S, Y, dK, s, sigG, c, mat, opts);

    if nSteps == 1 && ~any(active)
        status = 'noGrowth';
        N = Inf;
        break
    end

    if ~any(rate > 0)
        % Nothing can change any more: the state is frozen -> run-out.
        status = 'runout';
        N = Inf;
        break
    end

    % 3. increment size
    if fixed
        dN = opts.FixedStep;
    else
        dN = adaptiveStep(S, rate, cls, D, c, opts);
    end

    % 4. growth (adaptive: Heun predictor-corrector, 2nd order in dN)
    if ~fixed
        T = S;
        T.a = S.a + rate * dN;
        T.ext = S.ext + rate * dN;
        YT = flpm.internal.geometryFactors(T, D, opts, c.rScale);
        rateT = growthRates(T, YT, s * YT .* sqrt(pi * T.a), s, sigG, c, mat, opts);
        rate = 0.5 * (rate + rateT .* (rate > 0));
    end
    da = rate * dN;
    aOld = S.a;
    S.a = S.a + da;
    S.ext = S.ext + da;
    g = da > 0;
    switch opts.AspectUpdate
        case 'legacy'   % increment applied twice, as in the published code
            S.asp(g) = (S.a(g) + da(g)) ./ (S.a(g) ./ S.asp(g) + da(g));
        otherwise       % a/c of an ellipse whose semi-axes both grow by da
            S.asp(g) = (aOld(g) + da(g)) ./ (aOld(g) ./ S.asp(g) + da(g));
    end
    S.asp = min(max(S.asp, opts.AspectBounds(1)), opts.AspectBounds(2));
    if ~fixed
        N = N + dN;
    end
    if opts.RecordHistory
        hist = pushHistory(hist, N, max(S.a)*1e3, sigG, numel(S.a));
    end

    % 5. failure / run-out checks (order as in the published code)
    if S.a(ind) < D/2 && netStress(S.a(ind), s/(1 - c.R), c) > mat.flowStress
        status = 'netSection';
        break
    end
    if max(dK) / (1 - c.R) > mat.KIC
        status = 'Kc';
        break
    end
    if N > opts.MaxLife
        status = 'runout';
        N = Inf;
        break
    end
end

[~, k] = max(S.dK);
out.N = N;
out.status = status;
out.killer = S.id(k);
out.nSteps = nSteps;
out.final = struct('size', S.a*1e3, 'radius', hypot(S.x, S.y)*1e3, ...
    'angle', atan2(S.y, S.x), 'aspect', S.asp, 'extension', S.ext*1e3, ...
    'id', S.id);
if opts.RecordHistory
    out.history = finaliseHistory(hist);
end
end

% -------------------------------------------------------------------------
function c = constants(mat, geom, opts)
c.D = geom.D / 1000;
c.R = geom.R;
c.e = geom.e;
c.U = flpm.internal.newmanU(geom.R, mat.alpha, mat.newmanSmaxSigma0);
a0  = (1/pi) * (mat.dKth_lc  / (mat.Yw * mat.dSigma_e))^2;   % El Haddad
a0e = (1/pi) * (mat.dKth_eff / (mat.Yw * mat.dSigma_e))^2;
c.l = a0 - a0e;                                % closure build-up length l1
c.rpFactor = opts.PlasticZoneFactor;
switch lower(opts.InteractionRadius)
    case 'sqrtarea',         c.rScale = 1;
    case 'equivalentcircle', c.rScale = 1/sqrt(pi);
    otherwise
        error('flpm:opts:InteractionRadius', 'Unknown InteractionRadius.');
end
end

function [rate, active] = growthRates(S, Y, dK, s, sigG, c, mat, opts)
% NASGRO / Maierhofer growth rate of every defect (m/cycle).
rp = c.rpFactor * (dK / (1 - c.R) / mat.sigma0).^2;
dKg = sigG * Y .* sqrt(pi * (S.a + rp));
decay = exp(-S.ext / c.l);
dKth = mat.dKth_eff + (mat.dKth_lc - mat.dKth_eff) * (1 - decay);
F1 = 1 - (1 - c.U^mat.m) * (1 - decay);
F2 = 1 - (1 - c.U) * (1 - decay);
if strcmpi(opts.InitiationStress, 'net')
    dKcheck = sigG * Y .* sqrt(pi * S.a);
else
    dKcheck = dK;
end
active = F2 .* dKcheck >= mat.dKth_eff;
rate = zeros(size(S.a));
above = active & dKg > dKth;
rate(above) = mat.C1 * F1(above) .* dKg(above).^(mat.m - mat.p) ...
    .* (dKg(above) - dKth(above)).^mat.p;
creep = active & ~above & dKg < dKth;      % dKg == dKth: no growth (code)
rate(creep) = opts.CreepRate;
end

function S = netStress(a, s, c)
if a < c.D/2
    S = flpm.internal.netSectionStress(a, s, c.D, c.e);
else
    S = flpm.internal.netSectionStress(c.D/2*(1 - 1e-9), s, c.D, c.e);
end
end

function dN = adaptiveStep(S, rate, cls, D, c, opts)
% Largest increment such that no defect advances by more than Eta of its
% governing length (size, closure length while the R-curve builds up) and
% that coalescence and surface/near-surface transitions are resolved.
eta = opts.Eta;
g = rate > 0;
len = S.a;
building = S.ext < 5 * c.l;
len(building) = min(len(building), c.l);
dN = min(eta * len(g) ./ rate(g));

r = c.rScale * S.a;
vr = c.rScale * rate;
Rpos = hypot(S.x, S.y);
small = 0.01 * eta * r;                         % lets an event be crossed

% gaps between defects
n = numel(S.a);
if n > 1
    dx = S.x - S.x.';
    dy = S.y - S.y.';
    gap = sqrt(dx.^2 + dy.^2) - (r + r.');
    v = vr + vr.';
    allow = 0.5 * gap + min(small, small.');
    mask = v > 0 & ~eye(n, 'logical');
    if any(mask(:))
        dN = min(dN, min(allow(mask) ./ v(mask)));
    end
end
% approach to the free surface and to the near-surface criterion
emb = g & cls < 2;
if any(emb)
    lig = D/2 - Rpos(emb) - r(emb);
    dN = min(dN, min((0.5 * lig + small(emb)) ./ vr(emb)));
    q = opts.NearSurfaceRatio;
    rStar = q * (D/2 - Rpos(emb)) / (1 + q);
    internal = cls(emb) == 0;
    if any(internal)
        e = find(emb);
        e = e(internal);
        dN = min(dN, min((0.5 * (rStar(internal) - r(e)) + small(e)) ./ vr(e)));
    end
end
dN = min(max(dN, opts.MinStep), opts.MaxStep);
end

function h = initHistory(on)
h = struct('on', on, 'k', 0, 'data', zeros(on*1024, 4));
end

function h = pushHistory(h, N, amax, sigG, n)
h.k = h.k + 1;
if h.k > size(h.data, 1)
    h.data = [h.data; zeros(size(h.data))];
end
h.data(h.k, :) = [N, amax, sigG, n];
end

function out = finaliseHistory(h)
d = h.data(1:h.k, :);
out = struct('N', d(:,1), 'aMax', d(:,2), 'sigmaNet', d(:,3), 'nDefects', d(:,4));
end
