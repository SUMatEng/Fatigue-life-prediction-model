%% Example 1 - Fatigue strength of a single defect and an S-N curve
% One defect of given size and depth below the free surface.
% Units: mm for sizes and positions, MPa for stress.

mat  = flpm.material('Ti64-AB');               % Mater. Des. 2026, Table 2
opts = flpm.options('LoadRatio', 0.1, 'TargetLife', 1e7);

% sqrt(area) = 56 um, centre 0.256 mm below the surface
d = flpm.defectSet(0.056, 0.256, 0);

res = flpm.fatigueStrength(d, mat, opts);
fprintf('Fatigue strength at 1e7 cycles: %.1f MPa (%d life evaluations)\n', ...
    res.stress, res.nEval);

% S-N curve of the same defect
s = linspace(res.stress*1.01, 2*res.stress, 25);
N = arrayfun(@(x) flpm.fatigueLife(x, d, mat, opts).N, s);
figure; semilogx(N, s, 'o-'); grid on
xlabel('Cycles to failure'); ylabel('\Delta\sigma (MPa)');
title('Single defect, PBF-LB Ti-6Al-4V (AB)');

% Strength against depth: the near-surface band is the most severe
z = logspace(-2, 0.5, 30);
sw = arrayfun(@(q) flpm.fatigueStrength(flpm.defectSet(0.056, q, 0), mat, opts).stress, z);
figure; semilogx(z, sw, 'o-'); grid on
xlabel('Depth of the defect centre (mm)'); ylabel('\Delta\sigma_w (MPa)');
