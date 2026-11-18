function F = fitSizeDistribution(x, distName, varargin)
%FLPM.FITSIZEDISTRIBUTION  Fit a defect-size distribution truncated at the XCT detection limit.
%
%   F = flpm.fitSizeDistribution(sqrtArea, 'Lognormal', 'DetectionLimit', u)
%   pop = flpm.defectPopulation('volumetric', 'Size', F.dist, ...
%           'SizeBounds', [u Inf], 'CountPerArea', nDetected/area, ...)
%
%   Defects smaller than the detection limit u are not seen by XCT, so the
%   measured sizes are a sample of the distribution left-truncated at u.
%   The parameters are found by maximising
%       sum log f(x_i) - n log(1 - F(u)),
%   which removes the bias of an ordinary fit to truncated data. The fitted
%   distribution describes all defects; the detectable ones (density =
%   measured count per unit area or volume) follow it restricted to
%   x > u, which is what 'SizeBounds', [u Inf] does when sampling.
%   The untruncated fit is returned for comparison.
%
%   distName   any two- or three-parameter makedist name ('Lognormal',
%              'Weibull', 'Gamma', 'GeneralizedExtremeValue', ...)
%   'DetectionLimit'  u, mm (default: no truncation)
%
%   Output: dist (truncation-corrected), naive (ordinary fit), logL,
%   detectionLimit, fractionBelow (estimated share of defects below u).

ip = inputParser;
ip.addParameter('DetectionLimit', 0);
ip.parse(varargin{:});
u = ip.Results.DetectionLimit;
x = x(:);
x = x(isfinite(x));
if any(x < u)
    error('flpm:fitSize:limit', 'Observations below the detection limit.');
end
naive = fitdist(x, distName);
names = naive.ParameterNames;
q0 = naive.ParameterValues;
F.naive = naive;
F.detectionLimit = u;
if u <= 0
    F.dist = naive; F.logL = sum(log(pdf(naive, x))); F.fractionBelow = 0;
    return
end
nll = @(q) negLL(q, distName, names, x, u);
q = fminsearch(nll, q0, optimset('MaxFunEvals', 1e4, 'MaxIter', 1e4, 'Display', 'off'));
F.dist = mk(distName, names, q);
F.logL = -nll(q);
F.fractionBelow = cdf(F.dist, u);
end

function f = negLL(q, distName, names, x, u)
try
    d = mk(distName, names, q);
    S = 1 - cdf(d, u);
    l = pdf(d, x);
    if S <= 1e-12 || any(~isfinite(l)) || any(l <= 0)
        f = 1e10; return
    end
    f = -(sum(log(l)) - numel(x) * log(S));
catch
    f = 1e10;
end
end

function d = mk(distName, names, q)
args = [names; num2cell(q)];
d = makedist(distName, args{:});
end
