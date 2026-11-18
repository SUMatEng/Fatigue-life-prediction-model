function res = fatigueStrength(defects, mat, geom, opts)
%FLPM.FATIGUESTRENGTH  Stress range giving failure at the target life.
%
%   res = flpm.fatigueStrength(defects, mat, geom, opts)
%
%   Finds dSigma_w such that N_f(dSigma_w) = opts.TargetLife (1e7 by
%   default) for the given defect set.
%
%   Solver 'bracket' (default, 'paper' preset)
%     N_f(dSigma) decreases monotonically with dSigma and is infinite below
%     the stress at which no defect can grow. The lower end of the bracket is
%     that growth-onset stress, computed in closed form from the initial
%     geometry factors; the upper end is found by geometric expansion. The
%     root of g = log N_f - log N_target is then located by bisection,
%     accelerated by secant steps on (log dSigma, log N_f) once both ends of
%     the bracket are finite, to an absolute stress tolerance opts.StressTol.
%     Unlike the Nelder-Mead objective of the published code, this cannot
%     return a point outside the tolerance band.
%   Solver 'fminsearch' ('legacy' preset)
%     Reproduces Command_script.m: fminsearch on |TargetLife - N|, with
%     N = 0 for a run-out, error set to zero within 1e5 cycles,
%     optimset('TolFun',1e4,'TolX',2), start value opts.InitialStress.
%
%   Output struct
%     stress     dSigma_w (MPa); NaN if not found
%     converged  true if the target was bracketed to tolerance (bracket) or
%                the legacy error is zero (fminsearch)
%     N          cycles to failure at 'stress'
%     killer     original index of the killer defect at 'stress'
%     status     termination status of the life calculation at 'stress'
%     onset      growth-onset stress (MPa) (bracket solver)
%     nEval      number of life evaluations

if nargin < 4 || isempty(opts)
    opts = flpm.modelOptions();
end
switch lower(opts.Solver)
    case 'bracket'
        res = bracketSolve(defects, mat, geom, opts);
    case 'fminsearch'
        res = legacySolve(defects, mat, geom, opts);
    otherwise
        error('flpm:fatigueStrength:solver', 'Unknown solver "%s".', opts.Solver);
end
end

% -------------------------------------------------------------------------
function res = bracketSolve(defects, mat, geom, opts)
Nt = opts.TargetLife;
nEval = 0;
    function [g, out] = G(s)
        out = flpm.fatigueLife(s, defects, mat, geom, opts);
        nEval = nEval + 1;
        g = log(out.N) - log(Nt);          % +Inf for a run-out
    end

onset = growthOnset(defects, mat, geom, opts);
lo = onset;                                % N(lo-) = Inf for 'nominal'
[gLo, outLo] = G(lo);
k = 0;
while gLo < 0 && strcmpi(opts.InitiationStress, 'net') && k < 20
    lo = 0.9 * lo;                         % net-section onset is lower
    [gLo, outLo] = G(lo);
    k = k + 1;
end
if gLo < 0
    % N_f drops from infinity to below the target as soon as a defect can
    % grow (large or severely located defect): the S-N curve has a step at
    % the growth-onset stress, which is then the fatigue strength.
    res = pack(lo, true, outLo, onset, nEval);
    res.onsetLimited = true;
    return
end
hi = max(1.2 * lo, lo + 1);
[gHi, outHi] = G(hi);
k = 0;
while gHi > 0
    lo = hi; gLo = gHi; outLo = outHi;
    hi = 1.5 * hi;
    [gHi, outHi] = G(hi);
    k = k + 1;
    if k > 40
        res = pack(NaN, false, outHi, onset, nEval);
        return
    end
end

useSecant = true;
while hi - lo > opts.StressTol
    if useSecant && isfinite(gLo)
        % secant in (log s, g)
        x = exp(log(lo) + (log(hi) - log(lo)) * gLo / (gLo - gHi));
        if ~(x > lo && x < hi) || min(x - lo, hi - x) < 0.01 * (hi - lo)
            x = 0.5 * (lo + hi);
        end
    else
        x = 0.5 * (lo + hi);
    end
    useSecant = ~useSecant;                % alternate: guaranteed shrinkage
    [gx, outx] = G(x);
    if gx > 0
        lo = x; gLo = gx; outLo = outx;
    else
        hi = x; gHi = gx; outHi = outx;
    end
    if isfinite(gx) && abs(gx) < 1e-4    % |N/Nt - 1| < 0.01 %
        res = pack(x, true, outx, onset, nEval);
        return
    end
end
% Report the side that fails (finite N) at the target tolerance.
res = pack(hi, true, outHi, onset, nEval);
end

function res = pack(s, ok, out, onset, nEval)
res = struct('stress', s, 'converged', ok, 'N', out.N, ...
    'killer', out.killer, 'status', out.status, ...
    'onset', onset, 'nEval', nEval, 'onsetLimited', false);
end

function s = growthOnset(defects, mat, geom, opts)
% Smallest stress range at which any defect satisfies the growth check
% F2*dK >= dKth,eff on the first increment (F2 = 1 there). Coalescence of
% initially touching defects is taken into account.
S = flpm.internal.prepareState(defects, geom);
switch lower(opts.InteractionRadius)
    case 'sqrtarea',         rS = 1;
    otherwise,               rS = 1/sqrt(pi);
end
D = geom.D / 1000;
Y0 = flpm.internal.geometryFactors(S, D, opts, rS);
S.dK = Y0 .* sqrt(pi * S.a);
S = flpm.internal.coalesce(S, rS);
Y = flpm.internal.geometryFactors(S, D, opts, rS);
s = min(mat.dKth_eff ./ (Y .* sqrt(pi * S.a)));
s = s * (1 + 1e-9);
end

% -------------------------------------------------------------------------
function res = legacySolve(defects, mat, geom, opts)
Nt = opts.TargetLife;
nEval = 0;
    function e = err(s)
        nEval = nEval + 1;
        if s <= 0
            e = Nt; return
        end
        out = flpm.fatigueLife(s, defects, mat, geom, opts);
        N = out.N;
        if ~isfinite(N), N = 0; end
        e = abs(Nt - N);
        if e < 1e5, e = 0; end
    end
o = optimset('TolFun', 1e4, 'TolX', 2);
[s, e] = fminsearch(@err, opts.InitialStress, o);
out = flpm.fatigueLife(s, defects, mat, geom, opts);
res = struct('stress', s, 'converged', e == 0, 'N', out.N, ...
    'killer', out.killer, 'status', out.status, 'onset', NaN, ...
    'nEval', nEval, 'onsetLimited', false, 'legacyError', e);
end
