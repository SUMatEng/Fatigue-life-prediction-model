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

## Flaw interaction (continuous interaction factor, link-up on contact)

Each defect is the equivalent penny crack of radius r = √area/√π in the
plane normal to the load, centred at (position, depth, axial). Two
mechanisms act at every increment.

**1. Interaction (K amplification or shielding).** Every crack's ΔK is
multiplied by its interaction factor γ ≥ 0, which follows continuously
from the spacing, size ratio and offset of its neighbours. γ is computed
with Kachanov's method of average tractions (Kachanov, Int. J. Solids
Struct. 23 (1987) 23; Kachanov & Laures, Int. J. Fract. 41 (1989) 289):

```
p_i = 1 + Σ_j Λ_ij p_j,     Λ_ij = ⟨σ̂_j⟩ over crack i
t_i(P) = 1 + Σ_j p_j σ̂_j(P)          (traction on crack i)
K_i(Q) = 1/(π√(πr)) ∫ t_i(P) √(r² − ρ_P²)/|PQ|² dA      (Galin; Kassir & Sih)
γ_i = max_Q K_i(Q) / K_i,alone
```

σ̂_j is the exact normal stress of a penny crack under unit pressure
(Sneddon's closed form, `flpm.internal.pennyStress`, checked against
numerical Hankel integration to 5 digits). The weight-function integral is
evaluated with its singularity subtracted (the kernel integrates to 2π r)
on a 32 × 64 area-uniform grid at 24 front points.

| Configuration (equal cracks) | γ (this code) | Note |
|---|---|---|
| coplanar, edge gap 0.1 r | 1.29 | converged reference 1.292 |
| coplanar, edge gap 0.5 r | 1.06 | converged reference 1.064 |
| coaxial (stacked), 0.25 r apart | 0.74 | shielding; reference 0.744 |
| coaxial, 1 r apart | 0.83 | shielding |
| edge gap > 4 (r₁ + r₂) | 1 | cut-off, amplification < ~1 % |

The "converged reference" values are from the same formulation on a
4-times finer grid.

**Check against 3D FE (19 Nov 2026).** The cracks were modelled as thin
oblate cavities (half-thickness 0.03 r) in a block under remote tension,
using the MATLAB PDE Toolbox with linear tetrahedra
(`validation/fe/feCrackK.m`). K was taken from the crack-opening
displacement by displacement extrapolation. γ_FE is K of the pair divided
by K of the same crack alone, with the same mesh settings.

| Configuration (equal cracks) | γ, FE | γ, this code |
|---|---|---|
| single crack, K / Sneddon | 0.97 (±2.5 % around the front) | 1 |
| coplanar, edge gap 0.1 r | 1.27 / 1.25 | 1.29 |
| coplanar, edge gap 0.5 r | 1.07 / 1.06 | 1.06 |
| stacked, 0.25 r apart | 0.72 | 0.74 |

The factor agrees with FE to within 3 %, which is inside the FE accuracy.
It is kept as it is, with no correction. Gaps below 0.1 r could not be
meshed within the 16 GB of memory available. Such pairs link up after a
very small amount of growth, so they matter little for the life. Published SIF tables (Kachanov & Laures 1989) were not accessible; the
factor was instead checked against 3D FE (below). The method
is known to be accurate to a few percent down to gaps of about 0.1 r and
to under-predict the amplification at smaller gaps.

γ is re-evaluated whenever any crack has grown by more than 1 % since the
last evaluation.

**2. Link-up (coalescence).** Cracks are combined only when their
projections on the plane normal to the load touch and their planes are
aligned (the BS 7910 alignment rule, H ≤ 0.5 (l₁ + l₂) with l = 2r). The
linked crack is re-characterised as in BS 7910: the rectangle containing
both projections, idealised as its inscribed ellipse (half-axes c, h), and
stored as the equal-area crack √area = √(π c h) at the centre of the
rectangle. It takes the largest member extension; the member with the
largest ΔK names it. Link-up is transitive.

This replaces the earlier binary BS 7910 combination rules (s ≤ min(l₁, l₂)
side by side, s ≤ a₁ + a₂ stacked), which merge cracks across a gap and
are known to be conservative for AM porosity (ORNL, doi:10.2172/1985290,
~40 % higher K from the bounding flaw). The interaction across a gap is now
the computed γ.

**Limitations (to confirm).**

* The free surface is not included in the interaction. It enters only
  through Murakami's Y and A of each crack. A near-surface pair is
  therefore treated as in an infinite body, and the front points above
  the surface are not excluded from the maximum.
* Surface roughness valleys are surface discs in the same set, so the
  roughness–pore interaction enters through γ and link-up, with the same
  caveat.
* Cracks are taken as flat discs normal to the load. Tilted pores enter
  through their projected area (`flpm.projectedArea`), not their tilt.

## Residual stress (local mean-stress shift)

`flpm.options('ResidualStress', @(z) ...)` gives the residual stress
normal to the crack plane (MPa, tensile positive) as a function of depth
z (mm). At a defect centred at depth z:

```
Kmax = (σ_max + σ_res(z)) Y √(π A),   R_loc = (σ_min + σ_res)/(σ_max + σ_res)
```

The range ΔK is unchanged. Closure uses U(R_loc). The plastic zone uses
the local Kmax, and the K_IC check also uses it. A crack with
σ_max + σ_res ≤ 0 is closed throughout the cycle and does not grow. The
residual stress is evaluated at the defect centre: there is no gradient
over the crack and no relaxation. It is a hook for measured profiles
(e.g. hole drilling or XRD after machining or shot peening). Whether a
superposition through R alone is adequate for LPBF Ti-6Al-4V needs your
confirmation.

## Projected area from XCT morphology

`flpm.projectedArea` returns Murakami's √area on the plane normal to the
maximum principal stress n. It has two forms:

* From the three XCT projections (Ayz, Axz, Axy), taking the aligned
  ellipsoid with those projections:
  A(n) = √((n_x Ayz)² + (n_y Axz)² + (n_z Axy)²). This is exact for any
  shape when n is a scan axis.
* From the equivalent ellipsoid (semi-axes and principal directions):
  A(n) = π a b c √(nᵀ M n), with M = Q diag(a⁻², b⁻², c⁻²) Qᵀ.

## Finite life and the internal environment

**Growth-law calibration.** `flpm.fitGrowthLaw(dK, dadN, ...)` fits the
model's own NASGRO form to measured (ΔK, da/dN) pairs by least squares in
log da/dN:

* `'Closure', true` is for long-crack data in air. The law is
  F = U(R)^m with ΔKth = ΔKth,lc, which returns C1, m, p and ΔKth,lc.
* `'Closure', false` is for vacuum or internal-crack data. The law is
  F = 1 with ΔKth = ΔKth,eff, which returns C1 (with m and p) and the
  threshold.

The finite-life branch of the S–N curve (above the arrest step) follows
from this law. It is the part that the removed 1e-12 m/cycle creep used
to fill below the step. The model stays a fatigue-limit model: below the
step every crack arrests.

**Sealed (internal) cracks.** `mat.sealed` gives the constants of cracks
whose disc does not reach the free surface (z > r). Any of C1, ΔKth,eff
and ΔKth,lc can be given, and `closure = false` sets U = 1 for those
cracks. When a crack breaks the surface it switches to the air
constants. The default is empty, so there is no environment effect.

**Literature check (abstracts only; full texts not accessible here).**

* Yoshinaka et al., FFEMS 42 (2019), doi:10.1111/ffe.13085, used SR-µCT on
  wrought (α+β) Ti-6Al-4V in the VHCF regime. Internal crack growth rates
  were below 1e-10 m/cycle and fell to 1e-13 to 1e-11 m/cycle at
  ΔK < 5 MPa√m. These rates closely matched surface cracks in high
  vacuum, and the authors attribute the long internal-fracture lives to
  the vacuum-like environment inside the crack.
* Xue et al., FFEMS 45 (2022), doi:10.1111/ffe.13765, used nano-CT. The
  internal crack tip was blunted, as in vacuum, and its CTOD was almost
  independent of the load, unlike surface cracks.
* Yoshinaka et al., Int. J. Fatigue 142 (2021), doi:10.1016/j.ijfatigue.2020.105961,
  studied the effect of vacuum pressure on small-crack growth in
  Ti-6Al-4V. Only the title was seen.
* Vasudevan & Kujawski, FFEMS 45 (2022), doi:10.1111/ffe.13696, report that in
  vacuum ΔKth is nearly independent of R and that the small-crack effect
  is diminished. This motivates `closure = false` for sealed cracks.

**Flagged as unverified for this material.** The evidence is for wrought
Ti-6Al-4V with cracks initiated at α facets. LPBF pores contain the build
gas (argon), not vacuum. Whether a pore-initiated internal crack in LPBF
Ti-6Al-4V grows at vacuum rates has not been established here. No
vacuum constants are set by default. To use the effect, fit `mat.sealed`
from vacuum or internal-crack data, for example digitised curves from the
papers above or your own SR-CT data. Note that C1 changes only the life.
The fatigue strength at 10⁷ cycles changes only through the sealed
thresholds, because the strength is an arrest step. The closure length l
keeps its air value.

## Termination

* Failure is Kmax ≥ K_IC (Kmax = ΔK/(1 − R), shifted by any residual
  stress), evaluated with the K of each crack without amplification by its
  neighbours (γ clipped at 1). An amplified K ≥ K_IC on the side facing a
  neighbour breaks the ligament. That is link-up, which growth reaches
  within a few increments, and not fracture of the part. Without this
  rule, the life jumped by about 8 % with the step size, depending on
  whether a crack reached K_IC just before or just after linking. There is no net-section criterion,
  since that needs a cross-section.
* A run-out is either arrest of every crack, or N > 1.02 × the target life.

## Removed relative to the 2023–2025 code

* Specimen diameter, the Forman–Shivakumar and Alegre round-bar solutions,
  and the net-section stress with eccentric bending.
* The 1.0335 surface value and the 0.5/0.65 internal factor applied to the
  radius-based Alegre solution.
* 'legacy' and 'paper' presets, the fminsearch solver, fixed 100-cycle
  stepping, and creep growth of arrested cracks.

## Weakest link (closed form)

With crack arrest, a part fails at the smallest stress at which any defect
cannot arrest. `flpm.arrestMap` tabulates that stress for an isolated
defect, Δσ_w(√area, ζ = z/r), on a log-size × relative-depth grid. The
grid is split at Murakami's proximity limit, where Δσ_w jumps.
`flpm.weakestLink` then gives the strength distribution without sampling:

```
P_s(Δσ) = exp(−Λ(Δσ))              Poisson count
P_s(Δσ) = (1 − Λ(Δσ)/n)^n          fixed count n
Λ(Δσ)   = λ ∫ 1(Δσ_w(a, z)/f(x) ≤ Δσ) dF(a) dF(z) dx/|x|
```

f(x) = f(depth, position, axial) is the local stress range relative to
the nominal one, taken at the defect centre. That makes the result
transferable to notches, gradients and components. Several populations
multiply their P_s. The expectation is evaluated by quadrature at
mid-point quantiles of the truncated descriptors. The probability that
each quadrature node is the killer follows from the same expression, so
killer-size and killer-depth distributions come for free.

The closed form leaves out interaction and clustering, which only
`flpm.monteCarlo` includes. It is exact in the limit of no interaction,
and its difference from Monte Carlo measures the effect of interaction.

* **Flag (needs confirmation):** with a stress gradient, f is evaluated
  at the defect centre, and the crack-growth integration assumes uniform
  stress over the defect. For gradients that vary over a defect size, use
  a weight-function K. That is not implemented.

## Validation against S–N data and fractography

`flpm.fitSN` fits the random fatigue-limit model of Pascual & Meeker
(Technometrics 41 (1999) 277) by maximum likelihood, with run-outs
right-censored:

```
ln N = b0 + b1 ln(Δσ − γ) + ε,   ε ~ N(0, s²),   N = ∞ if Δσ ≤ γ
```

γ is either lognormal with fitted parameters, or fixed to the defect
model's fatigue-strength distribution (`'FatigueLimit', W`). In the second
case only b0, b1 and s are free, so `flpm.compareSN` (AIC, likelihood
ratio) tests whether the S–N data support the predicted fatigue-limit
scatter. The life model here is statistical; the defect model's own
finite-life prediction is stage 5.

`flpm.killerDistribution` gives the predicted killer √area and depth
distributions, from the weakest-link killer probabilities or the Monte
Carlo killers. With fractographic data it adds one-sample
Kolmogorov–Smirnov tests.

## Scatter sources

* **Threshold:** `mat.dKth_cov` is the coefficient of variation of
  ΔK_th,lc between specimens (lognormal, median ΔK_th,lc). Through the
  El Haddad length, it also scatters the closure length l. Monte Carlo
  draws it once per simulation from its own random stream. The weakest
  link integrates over it by 7-point Gauss–Hermite quadrature, with one
  arrest map per node. ΔK_th,eff is kept fixed. **Flag:** whether the
  scatter is better treated per specimen (as here) or per defect, as
  microstructural variation would suggest, needs your judgement.
* **Size–depth dependence:** `'Copula'` (Gaussian or t) on a volumetric
  population couples the size and depth draws. In the closed form it
  enters as joint quadrature weights. `flpm.testDependence` decides
  whether it is warranted (AIC and likelihood-ratio test against
  independence).
* **Detection limit:** `flpm.fitSizeDistribution` fits a size distribution
  left-truncated at the XCT detection limit u. Use it with
  `'SizeBounds', [u Inf]` and the count of detected defects.
* **Sampling uncertainty:** `flpm.bootstrapStrength` resamples the defect
  data, refits the descriptors with a user function, and returns
  confidence intervals of the strength quantiles. It uses the closed form.

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
* `'Intensity', 'radial'` (or a handle λ(x)): inhomogeneous K
  (Baddeley, Møller & Waagepetersen, Stat. Neerl. 54 (2000) 329), the
  inhomogeneous Thomas process (Waagepetersen, Biometrics 63 (2007) 252),
  and nulls thinned to λ(x). Use it for near-surface porosity, which is
  otherwise mistaken for clustering.
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
