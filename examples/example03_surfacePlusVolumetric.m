%% Example 3 - Surface roughness and volumetric defects combined
% Surface valleys are sampled with the peaks-over-threshold approach: the
% exceedances of Rv above the threshold u follow an exponential (GP with
% zero shape) distribution, and V exceedances per mm of profile are scaled
% to the circumference (Mater. Des. 2026, Eq. 6 and Sec. 4). The surface
% geometry factor is the roughness shielding factor Y_R (Eqs. 7-9).

mat  = flpm.material('Ti64-AB');
geom = flpm.geometry('Diameter', 7.3);

% Shielding factor from mean valley depth and mean pitch of the profile
YR = flpm.roughnessShielding(0.020, 0.160);          % Rv, RSm in mm

u = 0.020;                                           % EVT threshold, mm
rough = flpm.defectPopulation('surface', ...
    'Size', makedist('Exponential', 'mu', 0.010), 'SizeOffset', u, ...
    'CountPerLength', 0.35, ...                      % exceedances per mm
    'ShieldingFactor', YR, 'Name', 'AB surface');

pores = flpm.defectPopulation('volumetric', ...
    'Size',   makedist('Lognormal', 'mu', log(0.060), 'sigma', 0.20), ...
    'Depth',  makedist('Gamma', 'a', 2.0, 'b', 0.12), ...
    'Aspect', makedist('Beta', 'a', 5, 'b', 1.6), ...
    'Count',  12, 'Name', 'pores');

T = flpm.monteCarlo({pores, rough}, mat, geom, 'NumSimulations', 200, ...
    'Seed', 3, 'UseParallel', true);

fromSurface = mean(T.killerPopulation == 2);
fprintf('Median dSigma_w = %.0f MPa; %.0f %% of failures from the surface\n', ...
    median(T.stress), 100*fromSurface);
