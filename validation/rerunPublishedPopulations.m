%% Re-run the defect sets behind Mater. Des. 2026, Table 5
% Requires the unpublished sim_XCT_{siz,rad,ang}_dist{1,2,3}.xlsx and
% Output{1,2,3}.xlsx (2024_DefectPopulation/Nic's Model). Not part of the
% repository. The polar positions of the 7.3 mm bars are converted once
% to depth below the surface and arc length along it; the model itself
% uses no diameter.
dataDir = fullfile('..', '2024_DefectPopulation', 'Nic''s Model');
Rb = 7.3/2;
mat = flpm.material('Ti64-AB');
for dist = 1:3
    rd = @(q) readmatrix(fullfile(dataDir, sprintf('sim_XCT_%s_dist%d.xlsx', q, dist)));
    S = rd('siz'); R = rd('rad'); A = rd('ang');
    sets = arrayfun(@(i) flpm.defectSet(S(i,:), Rb - R(i,:), Rb*A(i,:)), ...
        (1:size(S,1))', 'UniformOutput', false);
    pub = readtable(fullfile(dataDir, sprintf('Output%d.xlsx', dist)));
    T = flpm.monteCarlo([], mat, 'DefectSets', sets, 'UseParallel', true);
    qp = quantile(pub.stress(pub.err == 0), [0.25 0.5 0.75]);
    q  = quantile(T.stress, [0.25 0.5 0.75]);
    fprintf('XCT volume %d: published %.0f (%.0f-%.0f), flpm %.0f (%.0f-%.0f) MPa\n', ...
        dist, qp([2 1 3]), q([2 1 3]));
end
