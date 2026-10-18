# Design notes

## 1. Coupling (ppiclF on CPU, host-staged)

nekRS owns the fluid on the device; ppiclF owns the particles on the host.
Two layers:

- `src/ppiclfCoupler.hpp` — thin helpers: pass the local nekRS mesh as
  ppiclF's overlap mesh, hand fields over for interpolation, read projected
  fields and particles. nekRS's element layout (`Nelements` blocks of `Nq^3`,
  i fastest) equals ppiclF's `(LEX,LEY,LEZ,NEE)`, so nothing is reordered.
  `PPICLF_LEX` must equal `polynomialOrder+1`; `PPICLF_LEE` must cover the
  local element count *and* every element ppiclF ships to a rank's bin.
- `src/elibm.hpp` (`elibm::Coupling`) with the ppiclF user routines
  `src/ppiclf_user_elibm.f` — the EL + IBM model below. A case only sets array
  sizes in its `PPICLF_USER.h` and includes `src/PPICLF_ELIBM.h`.

### Particles

ppiclF carries two kinds of particles (`rprop(JTYPE)`):

- **EL point particles**: Schiller–Naumann drag on the interpolated fluid
  velocity, gravity, and soft-sphere contact with other EL particles (ppiclF's
  binned neighbour search, `InitNeighborBin`) and with IBM bodies (direct loop;
  bodies are few). The drag reaction goes to the fluid.
- **IBM surface markers**: Fibonacci points on each sphere, moved rigidly with
  their body (snapped to the exact rigid position after every step). Each
  exerts the explicit direct-forcing force of NekIBM (`ibm_f_marker`),
  `f = rho_f dV (u_body - u)/dt`, on the fluid.

The body state (centre, velocity, angular velocity, rotation matrix) is
replicated on every rank in C++ and handed to Fortran each step. Bodies follow
Uhlmann (2005) with the internal fluid taken as rigid:
`(rho_b - rho_f) V dU/dt = -sum f + F_contact + (rho_b - rho_f) V g` and
`0.4 (rho_b - rho_f) V R^2 dw/dt = -sum r x f`.

Contact is a linear spring–dashpot along the normal; the dashpot is set from
the restitution coefficient for an isolated binary impact. The contact time
`pi sqrt(m_eff / k_n)` usually needs a smaller dt than the fluid, so particles
can take `substeps` per fluid step (fluid frozen). The EL particle stores the
contact force it receives from bodies (`rprop(JC*)`); summed over sub-steps it
is the impulse the bodies receive back.

### Fluid side

After the particle step `ppiclc_solve_Project` re-projects with the final state
(ppiclF's own last projection is one RK stage stale) onto seven fields:
EL force, IBM force, `phi_p`. Then

- momentum: `FU = F_EL / (1 - phi_p) + F_IBM`, a force per unit volume
  (`sumMakef` uses `JW * FU`), through `nrs->userSource`;
- continuity (`qtl`): `div u = -(1/phi_f) (D phi_f/Dt)`, `phi_f = 1 - phi_p`,
  with a backward difference in time and `opSEM::strongGrad`, through
  `nrs->userDivergence`. In a closed or fully periodic domain its volume mean is
  removed (it must integrate to zero).

IBM markers project with `phi_p = 0`: bodies are resolved, not volume-filtered.

### Time integration of the source

nekRS extrapolates explicit terms to t^{n+1} (EXT2: `2 E^n - E^{n-1}`) and its
BDF2 integrates them half a step late; together that is second order for the
smooth EL drag, which is passed through unchanged. For the IBM direct forcing,
a proportional controller on the marker slip, the extrapolation doubles the
gain: with gain `G` the slip obeys `z^2 - (1 - 2G) z - G = 0`, unstable at
`G = 1` (root -1.618). So the IBM part is written as
`E^n = (F - sum_{s>=1} c_s E^{n-s}) / c_0` with nekRS's current `c_s`
(`fluid->o_coeffEXT`), which makes the step apply exactly `F`; the slip then
decays as `(1 - G)^n`, stable for `0 < G < 2`.

### IBM gain and stability

The controller gain is what spreading then interpolating returns at a marker.
A Gaussian sheet integrates across the surface to `1 / (sqrt(2 pi) sigma)`,
so the marker volume is (surface share) x `sqrt(2 pi) sigma` (gain 1). Using the
grid spacing `h` instead, as Uhlmann does for his narrower delta, gives gain
`h / (sqrt(2 pi) sigma)` = 0.39 here and lingering slip.

The explicit body update adds a second constraint: the fluid shell the markers
drive, `rho_f A sqrt(2 pi) sigma`, must be lighter than the body's excess mass
`(rho_b - rho_f) V`, otherwise a change in body velocity is fed back amplified
each step. For a sphere: `rho_b / rho_f > 1 + 3 sqrt(2 pi) sigma / R`.
`elibm::Coupling::initialize` aborts below it (the collision demo at
`rho_b / rho_f = 2` < 2.28 diverged within 30 steps; 4 runs cleanly).

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

The numbers above are with `qtl = 0`; with `qtl = 1` they change in the fifth
digit (mean `phi_p` is small and nearly uniform). The imposed divergence has
RMS 2.0e-4 at t = 0.2; the solved `div u` differs from it by 6.0e-5 RMS, the
same as nekRS's own pointwise divergence error without QTL in this run
(5.9e-5), so the constraint is enforced to the discretisation error. Dense,
strongly non-uniform `phi_p` is still untested.

Verification, EL–IBM (`el_ibm_collision`, fully periodic unit box, 8^3 N=5
elements, IBM sphere R = 0.15, rho_b/rho_f = 4, initially at rest; 200 EL
particles d_p = 0.01, rho_p = 2000 (Stokes time 1.1) fired at it at 0.5;
k_n = 500, e = 0.9, 10 particle sub-steps per fluid step dt = 5e-3, 240 steps;
1 and 4 ranks agree to the last printed digit but one):

- no EL particle enters the sphere: the closest centre stays at
  0.47 d_p or more from the surface (touching is 0.5 d_p);
- the sphere is pushed from u = 0.03 (carried by the flow the cloud drags
  along) to 0.27 while the cloud hits it between t = 0.15 and 0.35; the first
  particles rebound (up to 14 with u_x < 0), later ones only slow down since
  the sphere is already moving;
- total x momentum (fluid + EL + body excess) drifts by -1.0 % of the initial
  EL momentum by t = 1.2; -0.4 % of it builds up before any contact (the
  `1/(1 - phi_p)` factor in the dense cloud), the rest from contacts and IBM
  (lagged forces, the rigid internal-fluid model, contact impulse taken at the
  last RK stage);
- IBM slip (RMS `|u_body - u|` at the markers) is 5e-3 at the end, under 2 %
  of the body speed, and peaks at 1.8e-2 during the impacts.

## 3. Known limitations

- Body contact impulses use the force of the last RK stage of each sub-step,
  not the RK-weighted average the EL particle actually received, so momentum
  exchange in contact is first order in the sub-step.
- No IBM–IBM contact and no wall contact; ppiclF's wall support
  (`InitWall`) is wrapped but not used by the EL + IBM model.
- Gravity acts on EL particles without buoyancy (heavy-particle limit).
- One Gaussian width serves both EL projection and IBM spreading (ppiclF has
  a single filter), so it must suit both: wider than the GLL spacing and
  narrow enough for the IBM interface.
- Bodies are replicated on every rank; fine for tens of bodies, not for
  thousands.

## 4. GPU roadmap (open)

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
