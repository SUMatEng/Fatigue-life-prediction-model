# Model (single version)

The `+flpm` package implements one model. There are no presets or
switches for alternative physics; `flpm.options` holds only the loading
(R, target life) and numerical settings. No specimen dimension enters.

## Defects

Each defect is idealised as a circular disc normal to the load (a pore
seen in projection on that plane):

* `size`: √area of the whole defect (mm). The equivalent-circle radius is
  r = √area/√π.
* `depth`: depth of its centre below the free surface (mm, ≥ 0).
* `position`: position of its centre along the free surface (mm). Only
  relative positions matter.
* `axial` (optional): position along the load axis (mm). If it is absent,
  all defects lie in one plane, which is the most severe case for
  interaction.

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

## Flaw interaction (BS 7910 rules)

Interaction is decided in 3D with the BS 7910 flaw-interaction rules,
applied at every step, so that cracks that grow into interaction are
combined. Each defect is treated as a circular flaw of length l = 2r and
height 2a = 2r (r = √area/√π).

1. **Alignment (non-coplanar flaws).** Two flaws in planes H apart, with
   in-plane gap g between their projections, are treated as coplanar when
   √(H² + g²) ≤ 0.5 (l₁ + l₂). Otherwise they do not interact.
2. **Combination (in the plane of projection).** With the containment
   rectangles of the projections:
   * overlapping projections interact;
   * side by side (overlapping in depth, gap s along the surface):
     s ≤ min(l₁, l₂);
   * stacked (overlapping along the surface, gap s in depth):
     s ≤ a₁ + a₂ (sum of half-heights);
   * flaws that are offset diagonally do not interact.
3. **Re-characterisation.** Interacting flaws are replaced by the rectangle
   that contains them, idealised as its inscribed ellipse (half-axes c, h).
   This is stored as an equal-area flaw, √area = √(π c h), at the centre
   of the rectangle, so that Murakami's √area and Y rules
   (above) apply to it unchanged. Its axial position is the area-weighted
   mean of the members. The merged crack takes the largest member
   extension, and the member with the largest ΔK names it.

Interaction is transitive: chains of pairwise-interacting flaws are
combined in one step, and the result is checked again until no pair
interacts.

**Sources and status.** The BS 7910 text itself was not consulted. The
rules above follow secondary sources:

* Kim (2007), Key Eng. Mater. 345–346, 411,
  doi:10.4028/www.scientific.net/KEM.345-346.411: alignment
  D ≤ 0.5(l₁ + l₂) and combination limits s ≤ min(l₁, l₂), s ≤ a₁ + a₂.
* TWI, overview of BS 7910:2013 and the history of its interaction
  criteria (twi-global.com). For surface flaws, the 2013 edition replaced
  "s = 0" with s ≤ 0.5 max(a₁, a₂), which is not used here because the
  half-space treatment of near-surface flaws is Murakami's.
* ORNL report on AM round bars, doi:10.2172/1985290. For porosity in AM
  bars, BS 7910 / ASME XI combination is conservative, with ~40% higher K
  from the bounding flaw. The report suggests s ≤ a_min instead.

Points to confirm against BS 7910:2019: the exact form of the alignment
and combination limits, and the rule for an embedded flaw near the free
surface. The latter is handled by Murakami's proximity rule (above), not
by BS 7910.

## Termination

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

## Defect clustering (Thomas process)

Defect locations are by default independent (a Poisson process). A
population can instead be clustered with `'Clustering', flpm.clustering(...)`
(a Thomas / Neyman–Scott process):

* Cluster centres form a Poisson process of intensity κ. Their depth comes
  from `ParentDepth` (default: the population's `Depth`). Their positions
  are uniform along the surface and along the axis.
* Each centre has Poisson(μ) defects (for a `'fixed'` count, the fixed
  total is shared multinomially among round(λ/μ) centres).
* Offsets are Gaussian with spread σ_plane in the cross-section (depth and
  position along the surface) and σ_axial along the load axis. Positions
  wrap around SurfaceLength and AxialLength, and depths are reflected at
  the surface.
* The population count stays the mean total λ, so κ = λ/μ. Sizes and
  other marks are drawn per defect, independently of location.

`flpm.fitClustering` estimates the process from measured centres:

* Ripley's K with translation edge correction for disc, box or cylinder
  windows. In 3D it uses the cylindrical K(r, t), which separates the
  in-plane and axial spreads:
  K(r, t) = 2tπr² + (1 − e^(−r²/4σp²)) erf(t/2σa)/κ.
* Minimum-contrast fit (Diggle, exponent 0.25) of κ and σ, with
  μ = λ/κ.
* Validation by a Monte Carlo deviation test on the border-corrected
  nearest-neighbour distribution G(r), against the fitted Thomas process
  and against Poisson.
* Mark tests: a random-labelling permutation test of nearest-neighbour size
  correlation, and the Spearman correlation of size with depth.

**Thinning.** Keeping only defects above a size threshold, with
independent sizes, leaves κ and σ unchanged and scales μ by the retained
fraction. A fit made on all detected pores must therefore be thinned
before it is used for a population of only the largest defects.

## Statistical sample size

`SurfaceLength` (mm) is the free-surface length that one Monte Carlo
simulation represents (πD for a round bar). It sets the defect count of
`CountPerLength` populations and the range of positions along the
surface. It is a sample size, not a geometry input to K. `AxialLength`
(mm) spreads defects along the load axis (3D sampling). `CountPerArea`
populations take CountPerArea × SurfaceLength × AxialLength defects.
