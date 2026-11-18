# Model (single version)

The `+flpm` package implements one model. There are no presets or
switches for alternative physics; `flpm.options` holds only the loading
(R, target life) and numerical settings. No specimen dimension enters.

## Defects

Each defect is a disc in the plane normal to the load:

* `size`: √area of the whole defect (mm). The equivalent-circle radius is
  r = √area/√π.
* `depth`: depth of its centre below the free surface (mm, ≥ 0).
* `position`: position of its centre along the free surface (mm). Only
  relative positions matter.

## Stress intensity (Murakami √area, half-space)

K = Y σ √(π A), where A is the effective √area.

| Class | Condition | A | Y |
|---|---|---|---|
| Internal | ligament h = z − r > r/0.8 | √area | 0.5 |
| Near-surface | 0 < h ≤ r/0.8 | √(π r (z + r)/2), the half-ellipse enclosing the defect and its ligament | 0.65 |
| Surface-breaking | z ≤ r | √(area of the disc below the surface) | 0.65 |

A is continuous where a disc touches the surface. Murakami's proximity rule
gives a step at h = r/0.8. That step is what makes a pore just inside the
proximity band the most severe position.

Roughness valleys of depth Rv with shielding factor Y_R, so that
K = Y_R σ √(π Rv) (Mater. Des. Eqs. 7–9), become surface-breaking discs
centred on the surface with A = (Y_R/0.65)² Rv. This is equivalent to the
same K. With no shielding (Y_R = 1.12), A ≈ 2.97 Rv, which is Murakami's
long shallow-notch limit.

## Crack growth (NASGRO with closure build-up)

```
da/dN = C1 F ΔK^(m−p) (ΔK − ΔKth)^p,     a = √area
F     = 1 − (1 − U^m)(1 − e^(−Δa/l))
ΔKth  = ΔKth,eff + (ΔKth,lc − ΔKth,eff)(1 − e^(−Δa/l))
ΔK    = Y Δσ √(π (A + rp)),   rp = (1/6π)(Kmax/σ0)²     [Acta Eq. 12]
```

* U is Newman's closure ratio. l is the difference of the El Haddad
  lengths for ΔKth,lc and ΔKth,eff with Y = 0.65.
* A crack starts to grow when F₂ΔK ≥ ΔKth,eff, with
  F₂ = 1 − (1 − U)(1 − e^(−Δa/l)).
* Once ΔK ≤ ΔKth the crack is arrested: it does not grow. The 1e-12 m/cycle
  creep of the 2023–2025 code is not used.

## Interaction and termination

* Discs that touch coalesce. Areas add, the centre is the area-weighted
  centroid, and the merged crack takes the largest member extension.
* Failure is Kmax = ΔK/(1 − R) ≥ K_IC. There is no net-section criterion,
  since that needs a cross-section.
* A run-out is either arrest of every crack, or N > 1.02 × the target life.

## Removed relative to the 2023–2025 code

* Specimen diameter, the Forman–Shivakumar and Alegre round-bar solutions,
  and the net-section stress with eccentric bending.
* The 1.0335 surface value and the 0.5/0.65 internal factor applied to the
  radius-based Alegre solution.
* 'legacy' and 'paper' presets, the fminsearch solver, fixed 100-cycle
  stepping, and creep growth of arrested cracks.

## Statistical sample size

`SurfaceLength` (mm) is the free-surface length that one Monte Carlo
simulation represents (πD for a round bar). It sets the defect count of
`CountPerLength` populations and the range of positions along the
surface. It is a sample size, not a geometry input to K.
