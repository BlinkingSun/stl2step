# How stl2step Reconstructs Analytic Geometry from Triangle Meshes

![An STL mesh (left) converts to a faceted B-Rep solid (center), and then to a
fully analytic solid (right).](assets/conversion-stages.png)

*Figure 1 — the two conversion stages. Left: the input STL, a triangle mesh.
Center: the Stage 1 result, a watertight faceted B-Rep solid. Right: the
Stage 2 result, the same part with every curved wall recovered as a true
analytic cylinder.*

## 1. The problem

An STL file describes a part as a bag of triangles. CAD and CAM systems,
however, reason about *boundary representation* (B-Rep) solids: faces carried
by analytic surfaces — planes, cylinders — bounded by analytic edges — lines,
circles. A triangle mesh of a machined part is unambiguous to the eye but
useless to a toolpath planner that needs to know "this wall is an arc of
radius 20 mm," not "this wall is 28 little rectangles."

Converting mesh to B-Rep therefore has two distinct levels of ambition, and
the engine keeps them as two distinct stages a user chooses between:

| Stage | Mode | What it produces | When to stop here |
|---|---|---|---|
| **1** | `verbatim` (default) | A byte-faithful, watertight B-Rep of the mesh exactly as tessellated. Coplanar triangles merge into single planar faces; nothing else is reinterpreted. | When geometric fidelity to the mesh is the requirement — inspection, archival, downstream repair, or any case where reinterpretation is unwanted. |
| **2** | `trueform` (`--engine trueform`) | The Stage 1 solid *plus* recognition and reconstruction: tessellated curved regions are identified and replaced by true analytic surfaces. | When the part must be machinable or editable — arcs as arcs. |

Stage 1 is always safe: it never invents geometry. Stage 2 is where the
interesting problem lives, and the rest of this document describes how it is
solved.

## 2. A tessellated arc is not noise — it is a pattern with a law

The key observation is that a mesh exported from a CAD system is not a noisy
sampling of the original surfaces. It is the *deterministic output of a
tessellation algorithm*, and that algorithm follows a law. Measured against
original CAD geometry, exported meshes show:

- Every vertex lies **exactly on the source surface** (measured deviations at
  the 10⁻⁵ mm scale — pure floating-point residue).
- A cylindrical face of radius *R* and angular extent *Φ* becomes a band of
  **equal-width flat strips joined corner to corner**: the tessellator caps
  the chordal deviation at some *d* and the facet-to-facet normal change at
  some angle *α*, then divides the arc into *N = ⌈Φ/θ_max⌉* equal steps of
  angle *θ = Φ/N*.

This yields a closed-form **inverse**. Given a chain of strips of chord width
*w* meeting at dihedral steps *θ*:

```
R = w / (2 · sin(θ/2))
```

Applied to correctly grouped strip chains, this recovers the source radius
essentially exactly — not "fitted within a few percent," but reproducing the
CAD value to the precision of the mesh itself. The axis follows from the
shared strip-edge directions, and the angular extent from the sum of the
steps.

Two disciplines make this robust in practice:

1. **The law's parameters are never assumed.** Different exporters and
   presets use different *d* and *α*. The engine estimates the law from the
   mesh's own band statistics, so recognition is parameter-free from the
   user's point of view.
2. **The law is also the validity test.** A chain is accepted as a
   tessellated arc only if its strips actually satisfy an equal-step law
   within tolerance — which is what separates true arcs from coincidental
   near-flat geometry.

## 3. Grouping, not fitting, is the hard part

A subtle failure mode dominates naive approaches: the radius mathematics can
be perfect while the *grouping* of triangles into candidate regions is wrong.
Generic region-growing segmentation makes two characteristic errors on
tessellated arcs:

- **Absorption.** A shallow arc — say 10° of a large radius, only a few
  strips wide — differs so little from a plane that it is swallowed into an
  adjacent planar region and never fitted at all.
- **Chimeras.** Strips from two different cylinders that happen to be
  adjacent get merged into one region, which is then "fitted" to a blended
  radius belonging to neither.

Stage 2 therefore inverts the usual order: the law-based chain detector runs
*first* and claims its strips before general segmentation can absorb them;
merges are permitted only between chains that share a generator edge, are
coaxial, agree in radius, and remain equal-step when concatenated; and any
chain whose per-strip radii disagree is split at the discontinuity rather
than averaged. Grouping quality is scored against labeled
triangle-to-surface correspondence data in the test suite, so regressions in
grouping are caught directly rather than inferred from radius coincidences.

## 4. Prismatic parts: solve the arcs in 2D, rebuild in 3D

A large class of machined parts is *prismatic* (2.5D): parallel top and
bottom faces, with every lateral wall parallel to one axis. For these parts
there is a strictly better strategy than repairing curved faces one at a
time in 3D:

1. **Detect** prismaticity from the recognition results themselves — all
   cylinder axes parallel, planes cleanly split into caps (normal along the
   axis) and laterals (normal perpendicular to it). Tolerances are computed
   from the mesh, not hard-coded.
2. **Slice** the part at its cap levels. Each level contributes a closed 2D
   profile — an outer loop and any hole loops — in the sketch plane.
3. **Fit the profile** as lines and arcs using the inverse law from §2. This
   is where the arcs are actually solved: in two dimensions, where chord
   chains are unambiguous and closure is checkable (every loop must close
   exactly, and its enclosed area must reconcile with the measured cap
   area). The fitted profiles can be emitted directly as DXF drawings
   (`--dxf <dir>`, one file per level), which doubles as a human-checkable
   view of exactly what was recognized — see Figure 1's center-to-right
   transition.
4. **Rebuild** the solid by extruding each level's profile and uniting the
   levels. Every curved wall now comes into existence as a true cylindrical
   face with correct shared edges *by construction* — there is nothing to
   repair afterward.

Parts that are not prismatic are routed, byte-identically, through the
general Stage 2 path or served by Stage 1; the detection gate makes this
choice explicit and testable.

## 5. Verification: every claim is a gate

Reconstruction earns trust only if it is checked by machinery, not by
inspection. Three families of checks run on every conversion and in the test
suite:

- **Volume accounting.** A true-arc solid legitimately differs in volume
  from its chord mesh by the arc-versus-chord defect, which is itself
  computable from the law: per strip, ΔV = (R²/2)(θ − sin θ)·h. Budgets are
  derived from that formula, so real reconstructions pass while phantom
  volume — the classic failure where an invalid shell reports a plausible
  face count with wildly wrong volume — is rejected.
- **Topological validity.** Shipped solids must be watertight and pass the
  kernel's boundary-representation checks; a closed-but-invalid shell is
  treated as failure regardless of how good its census looks.
- **Determinism.** Stage 1 output is byte-reproducible, and every optional
  feature (including DXF emission) is inert when unset. The conversion path
  contains no nondeterministic components.

## 6. Summary

The engine treats a tessellated mesh as *evidence about the CAD model that
produced it*. Stage 1 preserves that evidence faithfully as a solid. Stage 2
reads the evidence: it identifies the tessellation law the exporter followed,
inverts the law to recover exact radii, groups strips by the pattern rather
than by generic similarity, and — for prismatic parts — re-derives the solid
from 2D profiles so that curved walls are analytic by construction. The
result, verified by volume, validity, and supervised-recognition gates, is a
STEP file whose arcs are arcs.
