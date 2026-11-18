# Code review and rebuild notes

> **Status (18 Nov 2026).** At the owner's request the package now has a
> single, diameter-free model in the Murakami √area framework; see
> [MODEL.md](MODEL.md). Sections 1–5 below record the review of the
> 2023–2025 code that led to it. Section 0 gives the current results.

## 0. Current model on the published XCT defect sets

The 600 defect sets behind [M] Table 5 were converted once from polar
coordinates in the 7.3 mm bars to depth below the surface and arc length
along it. Results are median and interquartile range, in MPa:

| Volume | Published [M] | 2025 code, D = 7.3 mm, paper fixes | Contact-only coalescence | BS 7910, one plane | BS 7910, axial spread 2 mm |
|---|---|---|---|---|---|
| 1 (9 defects) | 222 (159–262) | 262 (245–283) | 303 (286–311) | 292 (268–309) | 305 (287–311) |
| 2 (12 defects) | 180 (149–203) | 238 (214–252) | 274 (261–284) | 257 (230–279) | 278 (265–285) |
| 3 (15 defects) | 158 (143–186) | 221 (203–246) | 278 (264–282) | 259 (239–277) | 279 (269–284) |

* The current model uses BS 7910 interaction (MODEL.md). The simulated XCT
  sets hold no axial coordinate, so "one plane" is how they are evaluated by
  default. That case is the upper bound on interaction. The last column
  gives the same sets random axial positions over 2 mm (seed 7). The
  results for 5 mm and 10 mm are the same within 1 MPa.
* With every defect in one plane, BS 7910 changes 60, 94 and 125 of the 200
  strengths per volume, by up to 133 MPa. Once the defects are spread
  along the axis, interaction almost disappears.

* All 600 sets converge, in about 0.5 s per set.
* **Every strength is a crack-arrest step.** Just below Δσ_w all cracks
  arrest. Just above it the bar fails by K_IC after 3.5–5 × 10⁵ cycles.
  Without the 1e-12 m/cycle creep of the old code, an arrested crack never
  restarts. The model therefore predicts a fatigue limit (the
  Kitagawa–El Haddad arrest stress) and gives no finite lives between about
  5 × 10⁵ and 10⁷ cycles. In the old code that creep, accumulating about
  10 µm over 10⁷ cycles, produced the lives near the target.
* The volume ranking weakens. Volumes 2 and 3 are equal within their IQRs,
  because the larger defect count of volume 3 no longer adds near-surface
  severity beyond what volume 2 already has.


### Interaction factor vs BS 7910 (18 Nov 2026)

The published XCT sets (200 per set) were rerun with the continuous
interaction factor and link-up on contact. They are compared with the
BS 7910 combination rules. Values are Δσ_w in MPa, given as median (P10).

| Set | Planar, BS 7910 | Planar, interaction factor | 2 mm axial spread, BS 7910 | 2 mm axial spread, interaction factor |
|---|---|---|---|---|
| 1 | 292.0 (245.7) | 301.3 (250.2) | 304.7 (249.9) | 304.8 (250.9) |
| 2 | 257.0 (220.4) | 272.3 (226.3) | 277.7 (227.5) | 278.4 (227.5) |
| 3 | 258.7 (214.3) | 271.9 (243.7) | 279.1 (256.5) | 279.3 (261.2) |

When all defects are in one plane, BS 7910 merges flaws across gaps.
That gives 9–15 MPa lower medians and up to 29 MPa lower P10, which is
the conservatism the ORNL report notes. With a realistic 2 mm axial
spread, the two methods agree to within 1 MPa in the median, because few
pairs are close enough for either. The interaction factor therefore
matters mainly for clustered or banded porosity. Its cost was about 18 min
for the 600 planar sets on 9 workers.

### Defect clustering (18 Nov 2026, revised)

Fits with `flpm.fitClustering`:

| Data | n | Window | κ | μ | σ (mm) | G test p (Thomas / Poisson) |
|---|---|---|---|---|---|---|
| Ti64 bars 21 / 22 / 23 behind [M] (`CT_excel.xlsx`, 3D), homogeneous null | 97 / 157 / 123 | cylinder R 3.6 mm × 9.9 mm | – | – | – | Poisson rejected (p = 0.005 each) |
| same, inhomogeneous null (`'Intensity', 'radial'`) | | | 0.23 / 0.38 / 0.27 /mm³ | 1.00 / 1.00 / 1.15 | – / – / [0.20 0.28] | 0.26 / 0.21 / 0.28 (Thomas), 0.17 / 0.21 / 0.04 (Poisson) |
| `XCTdefects.xlsx` (0.1 mm slab, 2D; not the [M] bars) | 98 | disc R 2.55 mm | = λ | 1.0 | – | 0.055 / 0.075 |
| NTU XY export (3D) | 27 541 | cylinder R 2.3 mm × 8.4 mm | 3.4 /mm³ | 58 | 0.34 in-plane, 0.22 axial | 0.025 / 0.025 |
| NTU ZX export (3D) | 9 451 | cylinder R 2.3 mm × 8.7 mm | (degenerate) | – | in-plane at the cap, 0.18 axial | 0.63 / 0.025 |

* **Ti64 bars of [M]:** the porosity is concentrated near the surface (median
  depth 0.25–0.34 mm in a 3.65 mm-radius bar). Against a homogeneous null
  that gradient looks like clustering. Against an inhomogeneous Poisson
  null, with the intensity a function of depth (which the population
  model already represents through its depth distribution), bars 21 and 22
  show no clustering (μ = 1). Bar 23 shows a weak excess of close pairs
  (Poisson p = 0.04, μ = 1.15, σ ≈ 0.2–0.3 mm). Size does not correlate
  with the neighbour's size in any bar. Size correlates with depth only in
  bar 23 (Spearman 0.22, p = 0.013).
* A first version of this section used `XCTdefects.xlsx`, a 0.1 mm slab
  that is not the data behind [M]. Its conclusion (no clustering) stands
  for the [M] bars, but for the reason given above.
* **NTU XY:** K exceeds the Poisson value by 15–25% between 0.1 and 0.5 mm,
  a mild clustering on the 0.3 mm scale. Below about 0.05 mm, G lies below
  both models (inhibition: pores do not overlap, and segmentation merges
  close pores). Neither the Thomas nor the Poisson process reproduces that
  range.
* **NTU ZX:** the excess pair density comes from banding along the
  specimen axis (σ_axial ≈ 0.18 mm), not from compact clusters. A Thomas
  process is a poor description of it, and a layered (Cox) intensity
  would fit better.
* **Thinning:** the NTU fits cover all detected pores (about 200 /mm³).
  The [M] sets keep only the 9–15 largest defects per simulation, i.e.
  about 0.2–0.3 /mm³ in the near-surface shell. With κ = 3.4 /mm³ the
  thinned mean cluster size is ≪ 1, so this clustering does not carry
  over to the largest defects. That would change if the large pores
  themselves cluster.

Strength (median and IQR, MPa), with the [M] sets resampled from their
pooled size and depth, SurfaceLength 22.9 mm, AxialLength 2 mm, BS 7910
interaction in 3D, 200 simulations:

| Volume | Poisson | Thomas μ = 3, σ = [0.34 0.22] mm | One cluster, σ = [0.34 0.22] mm |
|---|---|---|---|
| 1 (9 defects) | 304 (288–312) | 283 (250–305) | 275 (247–302) |
| 2 (12 defects) | 279 (267–285) | 234 (219–253) | 233 (212–254) |
| 3 (15 defects) | 279 (272–284) | 232 (215–253) | 226 (205–244) |

If the largest defects cluster on the scale seen in the NTU data, the
median falls by 20–50 MPa and the lower quartile by 40–65 MPa. The drop
comes almost entirely from BS 7910 interaction. Most of it appears as soon
as clusters of three form.


### Weakest link vs Monte Carlo (18 Nov 2026)

[M] sets resampled (pooled size and depth, fixed count, SurfaceLength
22.9 mm, AxialLength 2 mm). The table gives the 10/25/50/75/90% quantiles
of strength in MPa:

| Volume | Closed form (no interaction) | Monte Carlo (BS 7910, 3D) |
|---|---|---|
| 1 | 261 / 290 / 304 / 312 / 319 | 265 / 288 / 304 / 312 / 318 |
| 2 | 224 / 264 / 277 / 284 / 290 | 227 / 267 / 279 / 285 / 292 |
| 3 | 261 / 273 / 279 / 284 / 291 | 262 / 272 / 279 / 284 / 294 |

* The two agree within 3 MPa at every quantile, so with defects spread
  in 3D, interaction is negligible for these populations.
* The closed form takes 0.05 s per population after a 20 s map, against
  about 30 s for 200 Monte Carlo simulations. The predicted mean killer
  sizes are identical (0.062, 0.081 and 0.079 mm).
* Monte Carlo remains necessary when defects lie in one plane or cluster
  (Sections above).


### Scatter sources (18 Nov 2026)

* **Size–depth dependence** (`flpm.testDependence`, copula option): it is
  warranted only for bar 23 (Gaussian ρ = 0.21). Including it changes the
  volume 3 strength by at most 5 MPa (10% quantile 261 → 266 MPa, median
  unchanged). The simulated [M] sets carry no dependence (|ρ_s| < 0.02).
* **XCT truncation** (`flpm.fitSizeDistribution`): truncated and naive fits
  differ by < 3% in median and 99th percentile for these data, because the
  smallest recorded pores (26–37 µm) lie well below the median
  (50–60 µm). The estimated share of pores below the limit is 1–5%.
* **Bootstrap of the size descriptor** (`flpm.bootstrapStrength`, 200
  replicates of the 157 pores of bar 22, truncated lognormal refit, depth
  and count as for volume 2): estimate 233 / 275 / 289 MPa (10/50/90%),
  with 95% intervals [232, 235], [270, 283] and [285, 296] MPa. The
  sampling uncertainty of the size distribution is therefore about
  ±5–10 MPa on the median. That is smaller than the volume-to-volume
  differences in Table 5.
* **Threshold scatter**, specimen to specimen, CoV of ΔK_th,lc = 0.10
  (illustrative; no measured value is available). 10/50/90% quantiles in
  MPa:

| Volume | CoV 0 | CoV 0.10 |
|---|---|---|
| 1 | 261 / 304 / 319 | 265 / 303 / 324 |
| 2 | 224 / 277 / 290 | 225 / 276 / 297 |
| 3 | 261 / 279 / 291 | 259 / 279 / 298 |

  The scatter is dominated by the defect population. A 10% threshold
  scatter widens only the upper tail, by 5–7 MPa.

### Validation tools: first use (18 Nov 2026)

The S–N data of the XCT-characterised bars are not in the repository or
the project files. `S_N_data_TB.m` holds as-fabricated Ti-6Al-4V data
from this group (14 tests, amplitude given, no run-outs, R not stated)
and literature data. As a demonstration only, Wycisk et al. (as quoted
there: 14 tests, R = 0.1, four run-outs at 10⁷) was used:

* Free fit: median fatigue limit 185 MPa (s_g = 0.08).
* Fatigue limit fixed to the model prediction for the [M] volumes
  (medians 304, 277 and 279 MPa): strongly rejected (ΔAIC −52 to −2800).
  Run-outs at 175–198 MPa cannot coexist with a predicted limit near
  280 MPa.

This is a different material batch, so it does not validate or
invalidate the model. It does show that the current fatigue-limit
predictions sit well above typical as-built PBF-LB Ti-6Al-4V run-out
levels. They are also above the published [M] values (222/180/158 MPa).
The predicted killer √area (median 0.061–0.079 mm) and depth (median
0.26–0.32 mm) are ready for comparison with fractography.


Scope: the 2023 scripts in this repository, the February 2025 versions in the
2024_DefectPopulation project (`Nasgro.m`, `multi_coalition_version3.m`,
`Command_script.m`), and the two papers:

* **[A]** Macallister et al., Acta Mater. 237 (2022) 118189
* **[M]** Macallister, Deillon & Becker, Mater. Des. 266 (2026) 116162

The `+flpm` package keeps the 2023/2025 physics as the baseline. Where the
code and the papers disagree, the `'legacy'` preset follows the code and the
`'paper'` preset follows the papers **only for items 1 and 2 below**. All
other discrepancies are listed for a decision and are not changed.

## 1. Discrepancies between code and papers (decision needed)

| # | Item | Code | Paper | `+flpm` default |
|---|---|---|---|---|
| 1 | Plastic-zone size r_p | (π/6)(K_max/σ₀)² | (1/6π)(K_max/σ₀)², [A] Eq. 12 | paper (`PlasticZoneFactor`) |
| 2 | Aspect-ratio update | updated twice per increment, unbounded | single update, a/c ≤ 1 | paper (`AspectUpdate`, `AspectBounds`) |
| 3 | ΔK_th,lc and flow stress | 2.4 MPa√m, 1026 MPa (`Nasgro.m`) | 2.7 MPa√m, 1150 MPa ([M] Table 2) | paper values in `flpm.material('Ti64-AB')` |
| 4 | Surface geometry factor | max(Forman–Shivakumar, 1.0335) | 1.12 ([A] Sec. 2) | code |
| 5 | Near-surface / internal | a/h ≥ 0.8; internal Y × 0.5/0.65 | described inconsistently between [A] and [M] | code |
| 6 | Newman closure | coefficients 0.8613/0.3387/0.0465 | Newman (1984) gives different A₀–A₃ polynomial coefficients | code |
| 7 | Alegre embedded-crack Y | M₀ contains a duplicated λ² term; h = D/2 − R + a | source not verified against Alegre (2022) | code (flagged in `alegreEmbeddedY.m`) |
| 8 | Roughness √area_R | coefficient 2.97 | 0.97 ([M] Eq. 7) | 2.97 (Murakami); [M] Eq. 7 appears to be a typo |
| 9 | Initiation check | nominal stress | net-section stress ([M] Sec. 2) | code (`InitiationStress='nominal'`; `'net'` available) |
| 10 | Arrested cracks | grow at 1e-12 m/cycle | not described | code (`CreepRate`) |
| 11 | Interaction radius | √area used as circle radius | not stated explicitly | code (`InteractionRadius`) |
| 12 | Net-section stress in ΔK | σ_net (eccentric bending, e = 0.13) multiplies an already geometry-corrected Y | not discussed | code (`NetSectionInDeltaK`); risk of double-counting the finite-width correction |
| 13 | Coalesced crack extension | max of members' Δa | not described | code |

Items 1 and 2 are clear transcription errors in the code relative to [A]. They
change the strength by +0 to +5.5 MPa (median ≈ +3 MPa) in the 30-set
validation below. Items 3–13 change the physics and need your decision before
they are altered.

## 2. Published results: cause of the discrepancy found

`Output1–3.xlsx` (the source of [M] Table 5: median and IQR of the converged
rows) are **not** reproduced by the Feb 2025 `Nasgro.m`, nor by the Aug 2024
`MATLAB scripts 2024/Nasgro_many_01.m`, with either 2.4/1026 or 2.7/1150.
The current code is never lower than the published value and is markedly
higher (by more than 2 MPa in 519 of 600 rows, up to ~140 MPa).

**Cause.** The 2024 `multi_coalition_version3.m` has the specimen diameter as
its last argument, but `Nasgro_many_01.m` does not pass it. The commented line
`%D=(7.1/1000)` in that function indicates that the outputs were produced
with **D = 7.1 mm inside the coalescence/geometry-factor routine**, while the
rest of the model (net-section stress, K, failure checks) used 7.3 mm. With
D = 7.1 mm in `multi_coalition_version3` only, the 2024 code returns
N = 0.996–1.008 × 10⁷ at the published stress for every mismatched row tested
(Output1 rows 1, 5, 10, 33, 34). With D = 7.3 mm, nothing grows at those
stresses. Defects with centres 0.1–0.4 mm below the surface are then
classified and magnified against a free surface 0.1 mm closer than the
specimen's, which lowers their strength by up to ~140 MPa.

**Regenerated Table 5** (all 200 sets per volume; D = 7.3 mm throughout;
median and interquartile range, MPa):

| Volume | Published [M] (converged rows) | `legacy` physics | `paper` preset |
|---|---|---|---|
| 1 (9 defects) | 222 (159–262) | 258 (242–278) | 262 (245–283) |
| 2 (12 defects) | 180 (149–203) | 235 (212–248) | 238 (214–252) |
| 3 (15 defects) | 158 (143–186) | 220 (202–242) | 221 (203–246) |

All 600 sets converge with the bracket solver. In 52, 54 and 79 sets per
volume the strength is a step in N_f(Δσ) (arrest just below it). Runtime: about
20 s per 200 sets on 10 workers.

## 3. Single-defect check (`sanbox_multi`)

Defect √area = 0.056 mm, R = 3.394 mm, a/c = 1, target 10⁷ cycles:

| Implementation | Strength (MPa) |
|---|---|
| Original code + fminsearch, 2.7/1150 | 315.70 |
| Original code + fminsearch, 2.4/1026 | 302.34 |
| Earlier Python port | 306.3 |
| `+flpm` legacy (exact) | 315.70 (bitwise same cycle counts) |
| `+flpm` legacy physics, bracket solver | 315.6 ± 0.1 |
| `+flpm` paper preset | 321.4 |

## 4. Numerical and robustness issues fixed

* **Objective cliff in fminsearch.** The legacy objective |N_t − N| with
  run-out → N = 0 is flat and discontinuous. fminsearch stops on whichever
  plateau it first lands on, and the result depends on the start value. It
  was replaced with a bracketed secant/bisection on log N vs log Δσ that
  starts from the closed-form growth-onset stress. It needs 12–15 life
  evaluations (about 1 s per set, compared with about 50 s for the legacy fminsearch on the same engine).
* **S–N steps.** N_f(Δσ) can jump from run-out to well below the target life
  (arrest just below a stress, growth just above it). The solver reports the
  jump stress and flags it (`stepInSN`), or `onsetLimited` when the step is at
  the growth onset. In the 30-set validation, 7–14 of 30 sets per volume (23–47 %) show this.
* **Fixed 100-cycle increments** were replaced with adaptive increments (crack
  advance ≤ η·min(a, l, gaps), η = 0.005, Heun corrector). This gives the same
  life to < 0.5 % near the strength with 35× fewer steps. Near the threshold
  the life depends to first order on the step size because the active/arrest
  switch is discontinuous; fixed 100-cycle stepping is itself about 0.3 % off.
* **`Nasgro_many_01.m` does not run** (seven arguments passed, eight needed).
* Coalescence uses connected components, so chains A–B–C merge in one
  step regardless of defect ordering.
* Units: all internal quantities are SI; inputs are mm and MPa.
* Reproducibility: each simulation has its own RNG substream, so results do
  not depend on `parfor` scheduling or on the number of simulations drawn.
* Input validation for defects outside the bar, negative stresses, and empty
  sets (empty sets are assigned Δσ_e and flagged `noDefects`).

## 5. Speed

| | Original | `+flpm` |
|---|---|---|
| One life evaluation near 10⁷ cycles | ~39 s | ~0.02–1 s |
| One fatigue-strength solve | ~20 min | ~1 s |
| 30 XCT sets, both presets, 10 workers | n/a | 20 s |

## 6. Population descriptors

`flpm.defectPopulation` replaces the defect tables. Volumetric populations
take Size (√area), Radius or Depth, Aspect, Angle and Count. Surface
populations take Size (Rv, optionally an EVT exceedance plus `SizeOffset`),
`CountPerLength` V (count = V·πD, [M] Eq. 6) and `ShieldingFactor` Y_R. An EVT
`Threshold` thins sub-threshold defects. Sampling uses an exact truncated
inverse CDF.

## 7. Other notes

* An earlier configuration-file refactor exists locally
  (`~/Documents/MATLAB/defectTolerance`, March 2026). It was not merged.
* The unpublished XCT/profilometry data, `sim_XCT_*` and `Output*` files are
  **not** in this public repository. `validation/rerunPublishedPopulations.m`
  reads them from a local folder.
