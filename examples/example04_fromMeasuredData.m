%% Example 4 - Descriptors fitted to measured XCT data
% A data vector passed as a descriptor is fitted with a kernel density;
% only the fitted distributions enter the simulation. Depth is measured
% from the free surface: for an XCT scan of a round bar, the distance of
% each pore centre from the fitted outer contour.

T = readtable('CT_excel.xlsx', 'Sheet', 1);          % your XCT export
sqrtArea = sqrt(T.ProjectedArea_xy_plane__mm__);
depth = T.DistanceToSurface_mm_;                     % centre to free surface
scannedSurface = 22.9;                               % mm of surface in the scan section
nPerLength = numel(sqrtArea) / scannedSurface;

pop = flpm.defectPopulation('volumetric', 'Size', sqrtArea, 'Depth', depth, ...
    'CountPerLength', nPerLength);

R = flpm.monteCarlo(pop, flpm.material('Ti64-AB'), 'SurfaceLength', scannedSurface, ...
    'NumSimulations', 200, 'UseParallel', true);
disp(summary(R(:, {'stress', 'killerSize', 'killerDepth', 'killerSpacing'})))
