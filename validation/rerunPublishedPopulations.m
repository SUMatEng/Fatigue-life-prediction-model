%% Re-run the defect sets behind Mater. Des. 2026, Table 5
% Requires the (unpublished) simulated populations sim_XCT_{siz,rad,ang,asp}
% _dist{1,2,3}.xlsx and Output{1,2,3}.xlsx from 2024_DefectPopulation/
% Nic's Model. Set dataDir accordingly. These data are not part of the
% repository.
dataDir = fullfile('..', '2024_DefectPopulation', 'Nic''s Model');

mat  = flpm.material('Ti64-AB');
geom = flpm.geometry('Diameter', 7.3);
presets = {
    'legacy (adaptive, bracket)', flpm.modelOptions('legacy', 'Stepping', 'adaptive', ...
        'Eta', 0.005, 'Solver', 'bracket', 'MaxLife', 1.02e7)
    'paper',                      flpm.modelOptions('paper', 'Eta', 0.005)};

for dist = 1:3
    rd = @(q) readmatrix(fullfile(dataDir, sprintf('sim_XCT_%s_dist%d.xlsx', q, dist)));
    S = rd('siz'); R = rd('rad'); A = rd('ang'); P = rd('asp');
    sets = arrayfun(@(i) flpm.defectSet(S(i,:), R(i,:), A(i,:), P(i,:)), ...
        (1:size(S,1))', 'UniformOutput', false);
    pub = readtable(fullfile(dataDir, sprintf('Output%d.xlsx', dist)));
    fprintf('\nXCT volume %d (n = %d defects)\n', dist, numel(sets{1}.size));
    q = quantile(pub.stress(pub.err == 0), [0.25 0.5 0.75]);
    fprintf('  %-28s median %.1f (%.1f-%.1f) MPa\n', 'published (converged)', q([2 1 3]));
    for k = 1:size(presets, 1)
        T = flpm.monteCarlo([], mat, geom, 'DefectSets', sets, ...
            'Options', presets{k,2}, 'UseParallel', true);
        q = quantile(T.stress, [0.25 0.5 0.75]);
        fprintf('  %-28s median %.1f (%.1f-%.1f) MPa, %.1f s/set\n', ...
            presets{k,1}, q([2 1 3]), mean(T.runtime));
    end
end
