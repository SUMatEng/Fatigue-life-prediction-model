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
