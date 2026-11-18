%% Example 2 - Fatigue-strength distribution from population descriptors
% The defect population is described statistically - no defect list is
% read. Sizes, depths and aspect ratios are sampled inside flpm.monteCarlo.
% The descriptors below are illustrative (similar in scale to the PBF-LB
% Ti-6Al-4V XCT data of Mater. Des. 2026, Fig. 1 / Table 3).

mat  = flpm.material('Ti64-AB');
geom = flpm.geometry('Diameter', 7.3);

pores = flpm.defectPopulation('volumetric', ...
    'Size',   makedist('Lognormal', 'mu', log(0.060), 'sigma', 0.20), ... % sqrt(area), mm
    'Depth',  makedist('Gamma', 'a', 2.0, 'b', 0.12), ...               % centre below surface, mm
    'Aspect', makedist('Beta', 'a', 5, 'b', 1.6), ...                   % a/c
    'Count',  12, ...                                                   % per cross-section (N, Eq. 5)
    'Name',   'XCT pores');

% Optional: use a parallel pool (results are identical with or without it)
% parpool('Processes');
T = flpm.monteCarlo(pores, mat, geom, 'NumSimulations', 400, 'Seed', 1, ...
    'UseParallel', true);

q = quantile(T.stress, [0.25 0.5 0.75]);
fprintf('Median dSigma_w = %.0f MPa (IQR %.0f-%.0f), %d/%d converged\n', ...
    q(2), q(1), q(3), sum(T.converged), height(T));

figure;
subplot(1,2,1); histogram(T.stress, 25); xlabel('\Delta\sigma_w (MPa)'); ylabel('count');
subplot(1,2,2); scatter(T.killerDepth, T.killerSize, 20, T.stress, 'filled');
xlabel('killer depth (mm)'); ylabel('killer \surd area (mm)'); colorbar
