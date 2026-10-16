# Design notes

## 1. Current coupling (phase A: ppiclF on CPU, host-staged)

nekRS owns the fluid on the device; ppiclF owns the particles on the host.
Per fluid step, in `UDF_ExecuteStep(time, tstep)` (called after the step, so
fields are at t^{n+1}):

1. `ppiclf::interpField(mesh, jp, o_field)` — device → host copy of each field
   to interpolate, handed to `ppiclf_solve_InterpFieldUser`.
2. `ppiclc_solve_IntegrateParticle` — ppiclF's RK3 over one fluid dt.
3. (two-way, next) `ppiclf::projectedField(mesh, m, o_out)` — host → device copy
   of projected fields, added in `nrs->userSource` to `fluid->o_explicitTerms()`.

nekRS's element layout (`Nelements` blocks of `Nq^3`, i fastest) equals
ppiclF's `(LEX,LEY,LEZ,NEE)`, so the local nekRS mesh is passed directly as
ppiclF's overlap mesh (`ppiclf::initOverlapMesh`). `PPICLF_LEX` must equal
`polynomialOrder+1`, and `PPICLF_LEE` must cover the local element count *and*
every element ppiclF ships to a rank's bin.

Momentum source units: `fluid->o_explicitTerms()` (`FU`) enters `sumMakef`
as `JW * FU`, i.e. a force per unit volume. For the volume-filtered equations
the particle feedback is `F_proj / (1 - phi_p)` (as in ppiclF's Nek5000 example,
which divides by `rho` because Nek5000's `userf` is an acceleration).

## 2. ppiclF findings

Found while verifying `cases/rotation_tracer`; fixes live in `patches/ppiclF/`.

| # | Issue | Effect | Status |
|---|---|---|---|
| 1 | `-DPPICLC` renames routines to their C labels, breaking ppiclF's own internal calls (`InitParticle` → `AddParticles`) | UDF fails to `dlopen` | worked around: build without it, wrap in `src/ppiclf_cbind.f` |
| 2 | `IntegrateRK3` sets the loop bound `NPART*LRS` before `SetYdot`, which migrates particles | a particle arriving on a rank mid-step skips that RK stage (moved an extra `dt*u/3`) — multi-rank only | `0001-rk3-recount-particles-after-migration.patch` |
| 3 | `MapOverlapMesh` ships an element to a bin only if one of its GLL nodes lies in the bin | skewed/curved elements overlapping a bin with no node in it are missing; particles there are "not found" and get extrapolated values — silently | `0002-overlap-mesh-element-bbox.patch` (element bounding box, padded 1%) |
| 4 | Fluid fields are interpolated only at the first RK stage (`FinalizeInterp` resets the counter) | particle RK3 is first order in the fluid velocity seen along the path; for tracers it collapses to forward Euler | by design; fine for drag-dominated inertial particles, revisit for tracer-like particles |

Verification (`rotation_tracer`, frozen solid-body rotation on the unit-sphere
`ethier.re2` mesh, N=7, 1000 tracers, one revolution in 600 steps):
positions match ppiclF's discrete (Euler-equivalent) solution to 2e-14 on
1, 2, 4 and 8 ranks with no particle lost; error against the exact rotation is
1.7e-2, the O(dt) of item 4. Before patches 0001/0002, multi-rank runs were off
by 1e-2.

## 3. GPU roadmap (open)

The bottleneck of a GPU port of ppiclF as-is is not compute but the overlap-mesh
exchange: ppiclF's bins (particle bounding box, one bin per rank) are unrelated
to the fluid partition, so every RHS evaluation ships element coordinates,
fluid fields and projected fields (whole GLL elements) between fluid owners and
bin owners, and re-bins every stage.

Direction under consideration: keep nekRS's `lpm_t` ownership (particle lives
on the rank owning its element, so interpolation and projection are local) and
port only what `lpm_t` lacks from ppiclF:

- ghost particles within the filter width, sent to ranks whose element AABBs
  are within reach (host or GPU-aware MPI exchange of particles, not meshes);
- a device cell list (counting sort) for collisions and Gaussian projection;
- projection either as ppiclF's direct Gaussian sum (needs ghosts) or as a
  transpose-interpolation deposit + `oogs` dssum + Helmholtz (diffusion) filter
  (needs no ghosts; watch negative high-order weights for `phi_p`).

Open questions: load imbalance when particles cluster (element ownership has
no bin balancing), `lpm_t::maxEntriesPerParticleMigration = 50`, and how IBM
marker → body reductions (NekIBM queen/worker) map onto device segmented
reductions.
