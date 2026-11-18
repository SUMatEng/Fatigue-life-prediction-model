%% Example 4 - Descriptors fitted to measured XCT data
% Passing a data vector to a descriptor fits a kernel density (as in
% Mater. Des. 2026, Sec. 3). Only the fitted distributions are used by the
% model; the measured defects themselves are not simulated.

T = readtable('CT_excel.xlsx', 'Sheet', 1);          % your XCT export
sqrtArea = sqrt(T.ProjectedArea_xy_plane__mm__);
[~, radius] = cart2pol(T.CenterX_mm_ - mean(T.CenterX_mm_), ...
                       T.CenterY_mm_ - mean(T.CenterY_mm_));
sx = T.ProjectedSizeX_mm_; sy = T.ProjectedSizeY_mm_;
aspect = min(sx./sy, sy./sx);

pop = flpm.defectPopulation('volumetric', ...
    'Size', sqrtArea, 'Radius', radius, 'Aspect', aspect, 'Count', 9);

R = flpm.monteCarlo(pop, flpm.material('Ti64-AB'), flpm.geometry(), ...
    'NumSimulations', 200, 'UseParallel', true);
disp(summary(R(:, {'stress', 'killerSize', 'killerDepth', 'killerAspect'})))
