# Fatigue-life-prediction-model

Defect-tolerant fatigue-life and fatigue-strength model for additively
manufactured (PBF-LB) bars, after

* N. Macallister et al., *Acta Materialia* 237 (2022) 118189: NASGRO
  short-crack growth from multiple interacting defects;
* N. Macallister, F. Deillon, T.H. Becker, *Materials & Design* 266 (2026)
  116162: stochastic fatigue strength from volumetric (XCT) and surface
  (roughness) defect populations.

The `+flpm` package replaces the 2023–2025 scripts, which remain in the
history of the `main` branch. It is driven by **statistical population
descriptors** instead of measured defect lists. It has a single model in the
Murakami √area framework, located by distance from the free surface, with no
specimen diameter ([docs/MODEL.md](docs/MODEL.md)). The review of the earlier
code is in [docs/REVIEW.md](docs/REVIEW.md).

## Quick start

```matlab
addpath('path/to/Fatigue-life-prediction-model')   % the folder containing +flpm
mat = flpm.material('Ti64-AB');                      % Mater. Des. 2026, Table 2

% one defect set: sqrt(area), centre depth below the surface, position along it (mm)
d = flpm.defectSet([0.056 0.04], [0.26 0.30], [0 0.08]);
r = flpm.fatigueStrength(d, mat);                    % r.stress in MPa

% population descriptors -> strength distribution
pores = flpm.defectPopulation('volumetric', ...
    'Size',  makedist('Lognormal', 'mu', log(0.06), 'sigma', 0.2), ...
    'Depth', makedist('Gamma', 'a', 2, 'b', 0.12), ...
    'CountPerLength', 0.5);
rough = flpm.defectPopulation('surface', ...
    'Size', {'GeneralizedPareto', 'k', -0.29, 'sigma', 0.0147, 'theta', 0}, ...
    'SizeOffset', 0.050, 'CountPerLength', 0.37, 'ShieldingFactor', 1.02);
T = flpm.monteCarlo({pores, rough}, mat, 'SurfaceLength', 22.9, ...
    'NumSimulations', 400, 'Seed', 1);
median(T.stress)
```

The model needs no specimen dimension. `SurfaceLength` is only the amount
of free surface one simulation represents, i.e. the statistical sample size.
`'AxialLength'` (gauge length, mm) spreads defects along the load axis.
Flaw interaction then follows the BS 7910 rules in 3D, and
`'CountPerArea'` gives the count per mm² of surface. A measured set takes
axial positions with `flpm.defectSet(size, depth, position, 'Axial', y)`.

## Main functions

| Function | Purpose |
|---|---|
| `flpm.material`, `flpm.options` | Material parameters; loading ratio, target life, numerics |
| `flpm.defectPopulation` | Statistical descriptor of a volumetric or surface population |
| `flpm.samplePopulation` | Reproducible defect sets from descriptors |
| `flpm.defectSet` | A defect set from vectors (measured data) |
| `flpm.fatigueLife` | Cycles to failure at a stress range |
| `flpm.fatigueStrength` | Stress range giving the target life |
| `flpm.monteCarlo` | Strength distribution and killer-defect statistics |
| `flpm.roughnessShielding` | Y_R and √area_R from Rv and RSm |
| `flpm.arrestMap`, `flpm.weakestLink` | Closed-form strength distribution P_s = exp(−Λ), stress fields, killer probabilities |
| `flpm.fitSN`, `flpm.compareSN` | Random fatigue-limit S–N fit (censored run-outs); test of the predicted limit |
| `flpm.killerDistribution` | Predicted killer size/depth, KS tests vs fractography |
| `flpm.clustering` | Thomas cluster process for defect locations (`'Clustering'` of a population) |
| `flpm.fitClustering` | Fit κ, μ, σ to measured centres (Ripley's K, minimum contrast), G-test, mark tests |

The model is described in [docs/MODEL.md](docs/MODEL.md). Examples are in
`examples/`, and tests in `tests/` (`runtests('tests')`).

