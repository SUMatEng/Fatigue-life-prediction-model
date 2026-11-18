function res = fatigueStrength(defects, mat, opts)
%FLPM.FATIGUESTRENGTH  Stress range giving failure at the target life.
%
%   res = flpm.fatigueStrength(defects, mat)
%   res = flpm.fatigueStrength(defects, mat, flpm.options(...))
%
%   Finds dSigma_w with N_f(dSigma_w) = opts.TargetLife. N_f decreases with
%   dSigma and is infinite below the stress at which no crack can grow. The
%   lower end of the bracket is that growth-onset stress (closed form); the
%   upper end is found by expansion. The root of log N_f - log N_target is
%   then located by bisection alternated with secant steps, to opts.StressTol.
%
%   Output struct
%     stress        dSigma_w (MPa); NaN if not found
%     converged     true if bracketed to tolerance
%     N, killer, status   life calculation at 'stress'
%     onset         growth-onset stress (MPa)
%     nEval         number of life evaluations
%     onsetLimited  true if N_f < target as soon as a crack can grow
%     stepInSN      true if N_f jumps from run-out to below the target at
%                   'stress' (cracks arrest just below it)

if nargin < 3 || isempty(opts)
    opts = flpm.options();
end
Nt = opts.TargetLife;
nEval = 0;
    function [g, out] = G(s)
        out = flpm.fatigueLife(s, defects, mat, opts);
        nEval = nEval + 1;
        g = log(out.N) - log(Nt);          % +Inf for a run-out
    end

onset = growthOnset(defects, mat);
lo = onset;
[gLo, outLo] = G(lo);
if gLo < 0
    res = pack(lo, true, outLo, onset, nEval);
    res.onsetLimited = true;
    return
end
hi = max(1.2 * lo, lo + 1);
[gHi, outHi] = G(hi);
k = 0;
while gHi > 0
    lo = hi; gLo = gHi;
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
        lo = x; gLo = gx;
    else
        hi = x; gHi = gx; outHi = outx;
    end
    if isfinite(gx) && abs(gx) < 1e-4
        res = pack(x, true, outx, onset, nEval);
        return
    end
end
% Report the failing end. If N there is far from the target, N_f is
% discontinuous at the solution (arrest just below, failure just above).
res = pack(hi, true, outHi, onset, nEval);
res.stepInSN = abs(gHi) > log(1.02);
end

function res = pack(s, ok, out, onset, nEval)
res = struct('stress', s, 'converged', ok, 'N', out.N, ...
    'killer', out.killer, 'status', out.status, ...
    'onset', onset, 'nEval', nEval, 'onsetLimited', false, 'stepInSN', false);
end

function s = growthOnset(defects, mat)
% Smallest stress range at which a crack passes the growth check
% dK >= dKth,eff on the first increment, after the initial coalescence.
S = flpm.internal.prepareState(defects);
[Y, A] = flpm.internal.stressIntensity(S);
S.dK = Y .* sqrt(pi * A);
S = flpm.internal.coalesce(S);
[Y, A] = flpm.internal.stressIntensity(S);
Y = Y .* flpm.internal.interactionFactor(S);
s = min(mat.dKth_eff ./ (Y .* sqrt(pi * A))) * (1 + 1e-9);
end
