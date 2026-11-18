function C = compareSN(Ffree, Fmodel)
%FLPM.COMPARESN  Compare a model-predicted fatigue-limit distribution with S-N data.
%
%   Ffree  = flpm.fitSN(S, N, runout);                       % lognormal gamma
%   Fmodel = flpm.fitSN(S, N, runout, 'FatigueLimit', W);    % gamma from the model
%   C = flpm.compareSN(Ffree, Fmodel)
%
%   Both fits share the life model (b0, b1, s); they differ only in the
%   fatigue-limit distribution. C reports their log-likelihoods, AIC, the
%   AIC difference (positive: the data favour the model-predicted limit)
%   and the likelihood-ratio statistic LR = 2(logL_free - logL_model) with
%   its chi-square p-value for the two extra parameters. The model
%   distribution is not a member of the lognormal family, so the
%   chi-square reference is approximate; read the AIC difference first.
%   Also given: the median fatigue limit of each fit.

LR = 2 * (Ffree.logL - Fmodel.logL);
df = Ffree.k - Fmodel.k;
C = struct('logLfree', Ffree.logL, 'logLmodel', Fmodel.logL, ...
    'aicFree', Ffree.aic, 'aicModel', Fmodel.aic, ...
    'deltaAIC', Ffree.aic - Fmodel.aic, 'LR', LR, 'df', df, ...
    'p', 1 - chi2cdf(max(LR, 0), df), ...
    'gammaMedianFree', Ffree.gammaMedian, 'gammaMedianModel', Fmodel.gammaMedian);
end
