# Proposal: image-based FE of real XCT defect clusters

Status (19 Nov 2026): **superseded.** With no real geometries, the check
was reduced to idealised penny-crack pairs, which agree with the factor
to within 3 % (docs/MODEL.md). The real-geometry campaign below is not
pursued.

## Purpose

The closed-form model treats each pore as a flat penny crack normal to the
load. It interacts through Kachanov's average-traction factor γ, with no
free surface in the interaction. Image-based FE on segmented XCT clusters
would do three things:

1. **Calibrate γ.** It gives a correction γ_FE/γ_Kachanov for (a) small
   gaps (< 0.1 r), where Kachanov under-predicts, (b) near-surface pairs
   (free surface plus interaction), and (c) real 3D pore shapes (tilt,
   lack-of-fusion flatness).
2. **Validate the projected-area idealisation.** It tests √area on the
   plane normal to σ₁ against the K of the real pore with an initiated
   crack ring.
3. **Train an offline surrogate.** The surrogate would be
   γ = f(spacing/r, size ratio, axial offset/r, depth/r, aspect), used in
   place of, or as a correction to, the Kachanov factor. No FE would run
   inside the Monte Carlo.

## Selection of clusters (from data you already hold)

* CT_excel bars 21/22/23: the 30–50 pairs and triplets with the largest
  Kachanov γ, and the 10 near-surface pores with the smallest
  ligament/r. The selection is automatic from the existing code, by
  ranking `flpm.internal.interactionFactor` over all sets.
* NTU XY/ZX sets: 10 clusters from the dense banded regions (μ ≈ 58).
* Controls: 10 isolated pores spanning the size range, where γ should be 1.

That is about 100 geometries in total.

## Workflow

1. **Segmentation.** Crop sub-volumes of about 6 r around each cluster from
   the XCT data. Use the existing threshold, or Otsu with a local check.
   Write STL surfaces with Avizo/Dragonfly or MATLAB `isosurface`.
2. **Crack insertion.** Insert a thin annular crack of depth Δa ≈ 0.05 r
   (e.g. 5 µm) at the pore equator normal to σ₁, so that the pore is
   treated as the crack it becomes. Compute K along the ring with a
   domain (J) integral, either in Abaqus through `*CONTOUR INTEGRAL` on
   a seam or with an XFEM enriched ring in a submodel. Zencrack or
   FRANC3D would be the alternatives for tetrahedral meshes.
   - A pore-only stress-concentration analysis (K_t) is cheaper. It is
     useful for initiation but does **not** give K, so it is not the
     primary output.
3. **Boundary conditions.** Apply uniaxial σ in a block of about 10 r around
   the cluster. For near-surface clusters, use one traction-free face at
   the measured depth. Linear elastic analysis: E = 114 GPa, ν = 0.34.
4. **Output.** Record K_I(θ) per pore, giving γ_FE = max K_I/K_I,alone,
   where K_I,alone is from the same pore with its neighbours removed (the
   same mesh, so discretisation errors cancel).
5. **Verification first.** Before using any real geometry, reproduce the
   two-coplanar and two-coaxial penny cases in the MODEL.md table, and a
   surface semi-circular crack against Newman–Raju.

## Effort estimate (to be confirmed with the analyst)

| Step | Effort |
|---|---|
| Verification cases | 1 week |
| Segmentation and meshing pipeline, scripted | 2 weeks |
| ~100 clusters at ~1–3 h CPU each (Abaqus Standard, 1–3 M dof) | about 1 week of runs on a workstation, or 1–2 days on the UCT HPC |
| Surrogate fit and its use in `interactionFactor` | 1 week |

## Decisions needed from you

* Go or no-go, and who runs it (student, postdoc, or scripted here with an
  Abaqus licence on your machine or HPC).
* Whether pores are to be treated as cracks from the outset (the current
  model assumption), or whether a short initiation phase is wanted (K_t
  based).
* Access to the segmented XCT volumes (not only the Excel exports). These
  stay off the public repository.
