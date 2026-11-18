%% Example 1 - Fatigue strength of a single defect and an S-N curve
% Deterministic use of the solver: one defect, given size and location.
% Units: mm for sizes/positions, MPa for stress.

mat  = flpm.material('Ti64-AB');               % Mater. Des. 2026, Table 2
geom = flpm.geometry('Diameter', 7.3, 'LoadRatio', 0.1);
opts = flpm.modelOptions();                     % 'paper' preset

% sqrt(area) = 56 um, centre 0.256 mm below the surface, a/c = 1
d = flpm.defectSet(0.056, 3.394, 2.544, 1);

res = flpm.fatigueStrength(d, mat, geom, opts);
fprintf('Fatigue strength at 1e7 cycles: %.1f MPa (%d life evaluations)\n', ...
    res.stress, res.nEval);

% S-N curve of the same defect
s = linspace(res.stress*1.01, 2.5*res.stress, 25);
N = arrayfun(@(x) flpm.fatigueLife(x, d, mat, geom, opts).N, s);
figure; semilogx(N, s, 'o-'); grid on
xlabel('Cycles to failure'); ylabel('\Delta\sigma (MPa)');
title('Single defect, PBF-LB Ti-6Al-4V (AB)');
