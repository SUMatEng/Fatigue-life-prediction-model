function K = killerDistribution(src, varargin)
%FLPM.KILLERDISTRIBUTION  Predicted killer-defect size and depth, vs fractography.
%
%   K = flpm.killerDistribution(W)                 % flpm.weakestLink output
%   K = flpm.killerDistribution(T)                 % flpm.monteCarlo table
%   K = flpm.killerDistribution(..., 'Observed', [sqrtArea depth])  % mm
%
%   The predicted distribution of the defect that sets the strength:
%   from the killer probabilities of the weakest-link quadrature nodes, or
%   from the killers of the Monte Carlo simulations. With fractographic
%   measurements, one-sample Kolmogorov-Smirnov tests compare the observed
%   killer sizes and depths with the prediction.
%
%   Output struct
%     size, depth      tables x, cdf (prediction)
%     sizeQuantiles, depthQuantiles   5/25/50/75/95 %
%     ksSize, ksDepth  [D p] of the KS tests (if 'Observed' is given)
%
%   The prediction is unconditional on the test stress, which is exact for
%   a fatigue-limit model: the killer is the defect with the lowest
%   arrest stress. Fractography of specimens that failed at finite life is
%   comparable only as far as that holds.

ip = inputParser;
ip.addParameter('Observed', []);
ip.parse(varargin{:});
obs = ip.Results.Observed;
if istable(src)
    ok = isfinite(src.killerSize);
    a = src.killerSize(ok); z = src.killerDepth(ok); w = ones(nnz(ok), 1);
else
    a = src.killers.size; z = src.killers.depth; w = src.killers.weight;
end
w = w / sum(w);
K.size = wecdf(a, w);
K.depth = wecdf(z, w);
q = [0.05 0.25 0.5 0.75 0.95];
K.sizeQuantiles = wq(K.size, q);
K.depthQuantiles = wq(K.depth, q);
K.ksSize = [NaN NaN]; K.ksDepth = [NaN NaN];
if ~isempty(obs)
    K.ksSize = ks(obs(:, 1), K.size);
    K.ksDepth = ks(obs(:, 2), K.depth);
end
end

function T = wecdf(x, w)
[x, i] = sort(x(:)); w = w(i);
[xu, ~, j] = unique(x);
T = table(xu, cumsum(accumarray(j, w)), 'VariableNames', {'x', 'cdf'});
T.cdf = min(T.cdf, 1);
end

function v = wq(T, q)
v = arrayfun(@(p) T.x(find(T.cdf >= p - 1e-12, 1, 'first')), q);
end

function r = ks(xo, T)
xo = sort(xo(:)); n = numel(xo);
Fp = @(x) interp1([-Inf; T.x], [0; T.cdf], x, 'previous', 1);
Fo = (1:n)' / n;
Fm = Fp(xo);
Fl = arrayfun(@(x) interp1([-Inf; T.x], [0; T.cdf], x - eps(x), 'previous', 1), xo);
D = max([Fo - Fm; Fm - Fo + 1/n; Fl - (0:n-1)'/n]);
D = max(D, 0);
% asymptotic Kolmogorov p-value with Stephens' small-sample correction
lam = (sqrt(n) + 0.12 + 0.11/sqrt(n)) * D;
k = (1:100)';
p = min(max(2 * sum((-1).^(k-1) .* exp(-2 * k.^2 * lam^2)), 0), 1);
r = [D p];
end
