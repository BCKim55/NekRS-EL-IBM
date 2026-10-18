# Design notes

## 1. Current coupling (phase A: ppiclF on CPU, host-staged)

nekRS owns the fluid on the device; ppiclF owns the particles on the host.
Per fluid step, in `UDF_ExecuteStep(time, tstep)` (called after the step, so
fields are at t^{n+1}):

1. `ppiclf::interpField(mesh, jp, o_field)` — device → host copy of each field
   to interpolate, handed to `ppiclf_solve_InterpFieldUser`.
2. `ppiclc_solve_IntegrateParticle` — ppiclF's RK3 over one fluid dt.
3. Two-way: `ppiclc_solve_Project` re-projects with the final particle state
   (ppiclF's own last projection is one RK stage stale), `ppiclf::projectedField`
   reads the force and `phi_p` fields, and `FU = F_proj / (1 - phi_p)` goes to
   the device. `nrs->userSource` copies it into `fluid->o_explicitTerms()` for
   the next fluid step.

The particle side keeps the force it exerts on the fluid in its own `rprop`
slots (set in `SetYdot`, read in `MapProjPart`), not in `ydotc`, which ppiclF
zeroes every stage and which collision forces would also land in.

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

Found while verifying `cases/rotation_tracer` and `cases/settling_twoway`; fixes live
in `patches/ppiclF/`.

| # | Issue | Effect | Status |
|---|---|---|---|
| 1 | `-DPPICLC` renames routines to their C labels, breaking ppiclF's own internal calls (`InitParticle` → `AddParticles`) | UDF fails to `dlopen` | worked around: build without it, wrap in `src/ppiclf_cbind.f` |
| 2 | `IntegrateRK3` sets the loop bound `NPART*LRS` before `SetYdot`, which migrates particles | a particle arriving on a rank mid-step skips that RK stage (moved an extra `dt*u/3`) — multi-rank only | `0001-rk3-recount-particles-after-migration.patch` |
| 3 | `MapOverlapMesh` ships an element to a bin only if one of its GLL nodes lies in the bin | skewed/curved elements overlapping a bin with no node in it are missing; particles there are "not found" and get extrapolated values — silently | `0002-overlap-mesh-element-bbox.patch` (element bounding box, padded 1%) |
| 4 | Fluid fields are interpolated only at the first RK stage (`FinalizeInterp` resets the counter) | particle RK3 is first order in the fluid velocity seen along the path; for tracers it collapses to forward Euler | by design; fine for drag-dominated inertial particles, revisit for tracer-like particles |
| 5 | `ProjectParticleGrid` tests every particle against every GLL node of every element in its bin range | with one bin per rank, cost ~ particles × all local nodes (13 s/step for 10k particles, 216 N=5 elements) | `0003-projection-element-loop.patch`: skip elements whose bounding box is beyond the cutoff (bit-identical, 6× faster) |
| 6 | `ProjectParticleGrid` leaves the element loop (`exit`) at the first element whose lowest bin exceeds the particle's | assumes every element on a rank has a node in that rank's bin; false after 0002, so most elements were skipped (79 % of the force lost on 4 ranks) | `0003`: `cycle` instead of `exit` |
| 7 | `CreateGhost` drops a ghost when one was already sent to the same rank, unless the new one is periodic | with periodic bins narrower than twice the cutoff a rank needs both the unshifted and the periodic image; the unshifted one was dropped (1e-4 of the projected force on 4 ranks) | `0004-ghost-periodic-images.patch`: deduplicate on (rank, image position) |

Verification, one-way (`rotation_tracer`, frozen solid-body rotation on the unit-sphere
`ethier.re2` mesh, N=7, 1000 tracers, one revolution in 600 steps):
positions match ppiclF's discrete (Euler-equivalent) solution to 2e-14 on
1, 2, 4 and 8 ranks with no particle lost; error against the exact rotation is
1.7e-2, the O(dt) of item 4. Before patches 0001/0002, multi-rank runs were off
by 1e-2.

Verification, two-way (`settling_twoway`, 10 000 particles, rho_p/rho_f = 100,
d_p = 0.01, Stokes time 0.056, mean phi_p = 5.2e-3, mass loading 0.52, fully
periodic unit box, 6^3 N=5 elements, Gaussian filter FWHM 0.15 truncated at
1e-4 of its peak, dt = 1e-3, 200 steps; 1 and 4 ranks agree to all printed
digits):

- projection: `int F_proj dV` vs `sum f_p` differs by 3.6–4.6e-4. A 3D Gaussian
  cut where it falls to 1e-4 of its peak loses `P(chi2_3 > 18.42) = 3.6e-4`
  of its mass, so this is the truncation plus GLL quadrature, nothing else.
- fluid budget: `P_f - int_0^t int FU dV dt` is -1.6e-3 of the gravity impulse
  at t = 0.2, the explicit extrapolation of `FU` inside nekRS (EXT2 telescopes
  to O(dt (F_n - F_0))).
- global momentum: `P_f + P_p - sum(m_p) g t` is -9.6e-4 of the gravity impulse
  at t = 0.2 (1/(1 - phi_p) and the explicit lag).
- physics: the whole periodic system accelerates at `a = m_p g / (m_p + m_f)
  = 0.344`, so drag tends to `m_p (g - a) = 0.344`; the run reaches 0.341 at
  3.6 Stokes times with slip velocity 0.036 = `v_t (1 - a/g)`.

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
