# Fatigue-life-prediction-model

Defect-tolerant fatigue-life and fatigue-strength model for additively
manufactured (PBF-LB) bars, after

* N. Macallister et al., *Acta Materialia* 237 (2022) 118189: NASGRO
  short-crack growth from multiple interacting defects;
* N. Macallister, F. Deillon, T.H. Becker, *Materials & Design* 266 (2026)
  116162: stochastic fatigue strength from volumetric (XCT) and surface
  (roughness) defect populations.

The `+flpm` package is a refactor of the 2023 scripts (kept in the repository
root for reference). It is driven by **statistical population descriptors**
instead of measured defect lists, reproduces the original cycle counts when
run with the legacy options, and is about 40× faster per life evaluation and
50× faster per fatigue-strength solve. See [docs/REVIEW.md](docs/REVIEW.md)
for the code review, the deviations from the published model, and open
questions.

## Quick start

```matlab
addpath('path/to/Fatigue-life-prediction-model')   % the folder containing +flpm
mat  = flpm.material('Ti64-AB');                     % Mater. Des. 2026, Table 2
geom = flpm.geometry('Diameter', 7.3, 'LoadRatio', 0.1);

% one defect set -> fatigue strength at 1e7 cycles
d = flpm.defectSet(0.056, 3.39, 2.54, 1);            % sqrt(area), radius, angle, a/c
r = flpm.fatigueStrength(d, mat, geom);              % r.stress in MPa

% population descriptors -> strength distribution
pores = flpm.defectPopulation('volumetric', ...
    'Size',   makedist('Lognormal', 'mu', log(0.06), 'sigma', 0.2), ...
    'Depth',  makedist('Gamma', 'a', 2, 'b', 0.12), ...
    'Aspect', makedist('Beta', 'a', 5, 'b', 1.6), ...
    'Count',  12);
rough = flpm.defectPopulation('surface', ...
    'Size', {'GeneralizedPareto', 'k', -0.29, 'sigma', 0.0147, 'theta', 0}, ...
    'SizeOffset', 0.050, 'CountPerLength', 0.37, 'ShieldingFactor', 1.02);
T = flpm.monteCarlo({pores, rough}, mat, geom, 'NumSimulations', 400, 'Seed', 1);
median(T.stress)
```

Descriptors accept `makedist`/`fitdist` objects (including kernel fits), a
data vector (kernel fit), a constant, or a `makedist` argument cell.

## Main functions

| Function | Purpose |
|---|---|
| `flpm.material`, `flpm.geometry`, `flpm.modelOptions` | Parameters; option presets `'paper'` (default) and `'legacy'` |
| `flpm.defectPopulation` | Statistical descriptor of a volumetric or surface population |
| `flpm.samplePopulation` | Reproducible defect sets from descriptors (per-simulation RNG substreams) |
| `flpm.defectSet` | A defect set from vectors (measured data) |
| `flpm.fatigueLife` | Cycles to failure at a stress range |
| `flpm.fatigueStrength` | Stress range giving the target life (bracketed secant/bisection) |
| `flpm.monteCarlo` | Strength distribution, killer-defect statistics, optional `parfor` |
| `flpm.roughnessShielding` | Y_R and √area_R from Rv and RSm |

Examples are in `examples/`; tests in `tests/` (`runtests('tests')`).

## Legacy scripts

`Nasgro_many_01.m`, `multi_coalition_version3.m`, `net_stress.m`,
`sanbox_multi.m`: the 2023 code. `Nasgro_many_01.m` calls
`multi_coalition_version3` with seven arguments (the specimen diameter is
missing) and does not run as committed.
