function F = fitSN(stress, N, runout, varargin)
%FLPM.FITSN  Maximum-likelihood S-N fit with a random fatigue limit.
%
%   F = flpm.fitSN(dSigma, N, runout)
%   F = flpm.fitSN(dSigma, N, runout, 'FatigueLimit', W)   % limit from the model
%   p = F.cdf(dSigma, n)            P(N <= n | dSigma)
%   n = F.quantile(dSigma, p)       life quantile
%
%   Random fatigue-limit model (Pascual & Meeker, Technometrics 41 (1999)
%   277-289):
%     ln N = b0 + b1 ln(dSigma - gamma) + e,   e ~ N(0, s^2),   dSigma > gamma
%   and N = Inf if dSigma <= gamma. The fatigue limit gamma of each
%   specimen is random. Run-outs are right-censored at their N:
%     failure  L_i = E_gamma[ 1(gamma < S) phi(z_i)/s ]
%     run-out  L_i = E_gamma[ 1(gamma >= S) + 1(gamma < S)(1 - Phi(z_i)) ]
%   with z_i = (ln N_i - b0 - b1 ln(S_i - gamma))/s. The expectation over
%   gamma is a mid-point quantile rule.
%
%   gamma distribution
%     default             lognormal, parameters (mu_g, s_g) fitted
%     'FatigueLimit', W   the strength distribution of the defect model
%                         (flpm.weakestLink output, or a vector of
%                         simulated strengths, e.g. flpm.monteCarlo's
%                         T.stress). Only b0, b1 and s are fitted, so the
%                         fit tests the model's prediction of the
%                         fatigue-limit scatter against the S-N data.
%                         Specimens the model gives no failing defect
%                         (P_f < 1) have gamma = 'MatrixLimit'.
%
%   Name-value pairs
%     'FatigueLimit'  W struct or strength vector                      [[]]
%     'MatrixLimit'   gamma for defect-free specimens, MPa             [Inf]
%     'NumNodes'      quadrature nodes over gamma                      [400]
%
%   Output struct
%     params        b0, b1, s, (mu_g, s_g)
%     se            standard errors (observed information)
%     logL, k, aic  log-likelihood (of ln N), number of parameters, AIC
%     gammaMedian   median fatigue limit, MPa
%     cdf, quantile function handles
%     data          the input
%
%   See also flpm.compareSN, flpm.weakestLink, flpm.killerDistribution.

ip = inputParser;
ip.FunctionName = 'flpm.fitSN';
ip.addParameter('FatigueLimit', []);
ip.addParameter('MatrixLimit', Inf);
ip.addParameter('NumNodes', 400);
ip.parse(varargin{:});
o = ip.Results;
S = stress(:); n = N(:); ro = logical(runout(:));
if ~(numel(S) == numel(n) && numel(n) == numel(ro)) || any(S <= 0) || any(n <= 0)
    error('flpm:fitSN:input', 'stress, N and runout must be equal-length, positive.');
end
u = ((1:o.NumNodes)' - 0.5) / o.NumNodes;
fixed = ~isempty(o.FatigueLimit);
if fixed
    g = modelQuantiles(o.FatigueLimit, u, o.MatrixLimit);
    nll = @(q) -logLik(S, n, ro, q(1), q(2), exp(q(3)), g);
else
    nll = @(q) -logLik(S, n, ro, q(1), q(2), exp(q(3)), exp(q(4) + exp(q(5)) * sqrt(2) * erfinv(2*u - 1)));
end

% start: Basquin through the failures with gamma at 0.6 min(S)
fl = ~ro;
g0 = 0.6 * min(S);
if fixed, g0 = min(0.9 * min(S), median(g(isfinite(g)))); end
X = [ones(nnz(fl), 1), log(max(S(fl) - g0, 1))];
b = X \ log(n(fl));
s0 = max(std(log(n(fl)) - X * b), 0.2);
starts = [b(1), b(2), log(s0)];
if ~fixed
    starts = [starts, log(g0), log(0.1)];
    starts = [starts; starts + [0 0 0 log(1.3) 0]; starts + [0 0 0 log(0.75) log(2)]];
end
best = Inf;
for k = 1:size(starts, 1)
    [q, f] = fminsearch(nll, starts(k, :), optimset('MaxFunEvals', 2e4, 'MaxIter', 2e4, ...
        'TolX', 1e-9, 'TolFun', 1e-10, 'Display', 'off'));
    if f < best, best = f; qb = q; end
end
F.params.b0 = qb(1);
F.params.b1 = qb(2);
F.params.s = exp(qb(3));
if fixed
    gam = g;
else
    F.params.mu_g = qb(4);
    F.params.s_g = exp(qb(5));
    gam = exp(qb(4) + exp(qb(5)) * sqrt(2) * erfinv(2*u - 1));
end
F.logL = -best;
F.k = numel(qb);
F.aic = 2 * F.k + 2 * best;
F.se = sqrt(abs(diag(pinv(numHessian(nll, qb)))))';
F.gammaMedian = median(gam);
F.gamma = gam;
F.cdf = @(s, nn) lifeCdf(s, nn, qb(1), qb(2), exp(qb(3)), gam);
F.quantile = @(s, p) lifeQuantile(s, p, qb(1), qb(2), exp(qb(3)), gam);
F.fixedLimit = fixed;
F.data = table(S, n, ro, 'VariableNames', {'stress', 'N', 'runout'});
end

% =========================================================================
function g = modelQuantiles(W, u, gMatrix)
if isstruct(W)
    g = W.quantile(u);
else
    x = sort(W(:));
    g = x(min(numel(x), max(1, ceil(u * numel(x)))));
end
g = g(:);
g(isnan(g)) = gMatrix;
end

function L = logLik(S, n, ro, b0, b1, s, g)
g = g(:).';                              % 1 x nodes
live = S > g;                            % specimens x nodes
d = max(S - g, realmin);
z = (log(n) - b0 - b1 * log(d)) / s;
Lf = mean(live .* exp(-0.5 * z.^2) / (sqrt(2*pi) * s), 2);
Lr = mean(~live + live .* 0.5 .* erfc(z / sqrt(2)), 2);
Li = Lf;
Li(ro) = Lr(ro);
if b1 >= 0 || s <= 0
    L = -1e10; return
end
L = sum(log(max(Li, 1e-300)));
end

function P = lifeCdf(S, nn, b0, b1, s, g)
g = g(:).';
P = zeros(size(S));
for k = 1:numel(S)
    live = S(k) > g;
    z = (log(nn(min(k, numel(nn)))) - b0 - b1 * log(max(S(k) - g, realmin))) / s;
    P(k) = mean(live .* 0.5 .* erfc(-z / sqrt(2)));
end
end

function n = lifeQuantile(S, p, b0, b1, s, g)
n = zeros(size(S));
for k = 1:numel(S)
    f = @(ln) lifeCdf(S(k), exp(ln), b0, b1, s, g) - p;
    if lifeCdf(S(k), 1e30, b0, b1, s, g) < p
        n(k) = Inf;
    else
        n(k) = exp(fzero(f, [0 log(1e30)]));
    end
end
end

function H = numHessian(f, x)
k = numel(x); H = zeros(k);
h = 1e-4 * max(1, abs(x));
for i = 1:k
    for j = i:k
        ei = zeros(1, k); ej = ei; ei(i) = h(i); ej(j) = h(j);
        H(i, j) = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * h(i) * h(j));
        H(j, i) = H(i, j);
    end
end
end
