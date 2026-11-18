function D = testDependence(sz, depth)
%FLPM.TESTDEPENDENCE  Test size-depth dependence and fit a copula.
%
%   D = flpm.testDependence(sqrtArea, depth)
%   pop = flpm.defectPopulation('volumetric', ..., 'Copula', D)   % if warranted
%
%   Rank correlations (Spearman, Kendall) with their p-values, and maximum-
%   likelihood Gaussian and t copulas fitted to the pseudo-observations
%   (ranks/(n+1)), compared with independence by AIC and a likelihood-
%   ratio test. D.copula is the copula with the lowest AIC if it beats
%   independence by AIC and its likelihood-ratio test at the 5% level
%   (D.warranted); otherwise [].
%
%   Output struct
%     n, spearman [rho p], kendall [tau p]
%     gaussian, t   structs with rho (nu), logL, aic, pLR
%     aicIndependent (= 0), copula, warranted

x = sz(:); y = depth(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
n = numel(x);
if n < 10
    error('flpm:testDependence:input', 'Need at least 10 paired observations.');
end
D.n = n;
[r, p] = corr(x, y, 'Type', 'Spearman'); D.spearman = [r p];
[r, p] = corr(x, y, 'Type', 'Kendall');  D.kendall = [r p];
U = [tiedrank(x), tiedrank(y)] / (n + 1);

rg = copulafit('Gaussian', U);
lg = sum(log(copulapdf('Gaussian', U, rg)));
D.gaussian = struct('rho', rg(1, 2), 'logL', lg, 'aic', 2 - 2*lg, ...
    'pLR', 1 - chi2cdf(2 * max(lg, 0), 1));
try
    [rt, nu] = copulafit('t', U);
    lt = sum(log(copulapdf('t', U, rt, nu)));
    D.t = struct('rho', rt(1, 2), 'nu', nu, 'logL', lt, 'aic', 4 - 2*lt, ...
        'pLR', 1 - chi2cdf(2 * max(lt, 0), 2));
catch
    D.t = struct('rho', NaN, 'nu', NaN, 'logL', -Inf, 'aic', Inf, 'pLR', 1);
end
D.aicIndependent = 0;
[best, k] = min([0, D.gaussian.aic, D.t.aic]);
pk = [1, D.gaussian.pLR, D.t.pLR];
D.warranted = k > 1 && best < 0 && pk(k) < 0.05;
if ~D.warranted, k = 1; end
switch k
    case 1, D.copula = [];
    case 2, D.copula = struct('family', 'gaussian', 'rho', D.gaussian.rho, 'nu', Inf);
    case 3, D.copula = struct('family', 't', 'rho', D.t.rho, 'nu', D.t.nu);
end
end
