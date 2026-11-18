# Code review and rebuild notes

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

## 2. Published results are not reproduced by the current code

`Output1–3.xlsx` (the source of [M] Table 5) were rerun with the Feb 2025
`Nasgro.m` on the same `sim_XCT_*` defect sets, with both 2.4/1026 and
2.7/1150:

* About half of the rows match to < 1 MPa (rows 3, 4, 7, 8, 13, 14, 19 of
  Output1, …). The rest give **2–82 MPa higher** strengths (e.g. row 1:
  published 201 MPa, current code 254 MPa). At the published stress, rows 1–2
  give no growth or a run-out, so the published value is not a solution of
  this code.
* No pattern by surface/near-surface classification was found.

Conclusion: Output1–3 were generated with a different code version (or
different inputs) from the one in the repository. Before the `+flpm` results
are compared with [M], the version used for the paper should be located.

Validation, first 30 defect sets of each XCT volume (`'legacy'` physics with
adaptive stepping and bracket solver, and `'paper'` preset):

| Volume | Published median, 30 sets | `legacy` median | `paper` median | Rows within 2 MPa of published |
|---|---|---|---|---|
| dist1 | 229.0 | 254.0 | 257.9 | 11/30 |
| dist2 | 191.5 | 237.3 | 240.5 | 4/30 |
| dist3 | 172.2 | 218.7 | 218.8 | 1/30 |

The full sets (200 per volume) were not run, on request; at about 1 s per
set they take about 3–4 min per volume per preset on 10 workers.

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
