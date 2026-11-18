%% Example 2 - Fatigue-strength distribution from population descriptors
% The defect population is described statistically; no defect list is read.
% SurfaceLength is the length of free surface one simulation represents
% (pi*D for a round bar). It sets how many defects a simulation contains,
% not the stress-intensity solution.

mat = flpm.material('Ti64-AB');

pores = flpm.defectPopulation('volumetric', ...
    'Size',  makedist('Lognormal', 'mu', log(0.060), 'sigma', 0.20), ... % sqrt(area), mm
    'Depth', makedist('Gamma', 'a', 2.0, 'b', 0.12), ...               % centre depth, mm
    'CountPerLength', 0.5, ...                                         % per mm of surface
    'Name',  'XCT pores');

% parpool('Processes');   % optional; results are identical without it
T = flpm.monteCarlo(pores, mat, 'SurfaceLength', pi*7.3, ...
    'NumSimulations', 400, 'Seed', 1, 'UseParallel', true);

q = quantile(T.stress, [0.25 0.5 0.75]);
fprintf('Median dSigma_w = %.0f MPa (IQR %.0f-%.0f)\n', q(2), q(1), q(3));

figure;
subplot(1,2,1); histogram(T.stress, 25); xlabel('\Delta\sigma_w (MPa)'); ylabel('count');
subplot(1,2,2); scatter(T.killerDepth, T.killerSize, 20, T.stress, 'filled');
set(gca, 'XScale', 'log');
xlabel('killer depth (mm)'); ylabel('killer \surd area (mm)'); colorbar
