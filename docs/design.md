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
| 8 | `CreateBin` loops forever (`do while (bins_dx < d2chk)` with `n_bins` already 1) when a periodic direction of the domain is shorter than the filter cutoff | a thin periodic slab with a filter wider than its thickness hangs at start-up with no message | `0005-abort-on-short-periodic-length.patch`: abort with `Periodic length < cutoff, dir`; the ghost logic holds one periodic image, so the limit itself stays (keep every periodic length ≥ 2 x cutoff) |
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
- total x momentum (fluid + EL + body excess) drifts by -0.3 % of the initial
  EL momentum by t = 1.2 (-1.0 % before the body contact impulse was made
  RK-weighted); -0.4 % of it builds up before any contact (the
  `1/(1 - phi_p)` factor in the dense cloud), the rest from the IBM coupling
  (lagged forces, the rigid internal-fluid model);
- IBM slip (RMS `|u_body - u|` at the markers) is 5e-3 at the end, under 2 %
  of the body speed, and peaks at 1.8e-2 during the impacts.

## 2b. Spray step 1: parcels and injection (`cases/spray_jet2d`)

Both paths treat an EL particle as a parcel of `N_p` droplets of diameter
`d_p`: the equation of motion is per droplet (drag + gravity, nothing else
is needed at `rho_l / rho_g ~ 1000`), the force on the fluid and the liquid
volume are per parcel. The ppiclF path stores the parcel liquid volume in
`rprop(JVOLP)` and derives `N_p` from it; the device path has a `np`
property. Parcels are created during the run (`Coupling::inject`) from a
deterministic generator shared by both cases, so the two runs see identical
input. No breakup, deformation or dense-regime drag yet.

Verification, quasi-2D slot jet: box 2 x 1 x 0.1 (40 x 20 x 1 elements,
N = 5; periodic x, z; free-slip y, so the x momentum of gas + parcels equals
the injected momentum exactly), gas `rho = 1`, `nu = 5e-4`, slot width 0.02
at x = 0.1, liquid `rho = 1000`, `d_p = 0.002` (tau_p = 0.44, Re_p0 = 4),
`u_jet = 1`, liquid volume fraction 0.01 in the slot (mass loading 10 at the
nozzle), 10 parcels per step of 2.4 droplets, dt = 2e-3, 1000 steps, no
gravity, 4 ranks:

| | ppiclF path | device path |
|---|---|---|
| single droplet vs drag ODE, max error / U0 | 7e-8 (RK3) | 7e-8 (SSP-RK3, see §4) |
| momentum budget at t = 2 (max over the run) | +8.6e-4 (2.8e-3) | +9.4e-4 (2.8e-3) |
| fluid budget `P_gas - int FU` | 1.6e-4 | 1.6e-4 |
| P_gas, P_parcels at t = 2 | 1.0720e-2, 2.9314e-2 | 1.0806e-2, 2.9232e-2 |
| centreline gas u_x, max | 0.8333 | 0.8336 |
| centreline difference | — | max 0.0055, rms 0.0013 |
| step time | 1.25 s (10 k parcels) | 0.19 s |

The two implementations (ppiclF's bins, ghosts and Gaussian sum vs `lpm_t`,
L2 projection and diffusion filter) agree to 0.8 % in momenta and 0.7 % in
the centreline velocity at the jet head, with the momentum budget closed to
1e-3 of the injected momentum in both. With lpm_t's RK2 and the end-of-step
force the device path gave 1.5e-5 on the droplet check and a centreline
difference of 0.0145; the RK-weighted force brought the `settling_gpu`
momentum drift from -5.6e-3 to -1.2e-3 of the gravity impulse (ppiclF path:
-1.0e-3). At dt = 5e-3 both runs blew up at
t = 0.5 when the induced gas velocity (0.83 u_jet) pushed the gas CFL past 1:
nekRS's own advection limit, not the coupling.

## 3. Known limitations

- No IBM–IBM contact and no wall contact; ppiclF's wall support
  (`InitWall`) is wrapped but not used by the EL + IBM model.
- Gravity acts on EL particles without buoyancy (heavy-particle limit).
- One Gaussian width serves both EL projection and IBM spreading (ppiclF has
  a single filter), so it must suit both: wider than the GLL spacing and
  narrow enough for the IBM interface.
- Bodies are replicated on every rank; fine for tens of bodies, not for
  thousands.
- The body contact impulse is the RK-weighted integral of the force the EL
  particle received (`ppiclc_elibm_BeginStep`, `rprop(JI*)`), so EL–body
  momentum exchange is now consistent with the particle update; the body
  itself is still advanced with forward Euler.

## 4. EL particles on the GPU (`src/elgpu.hpp`, prototype)

### What nekRS already provides

`lpm_t` keeps the particle state on the device and does findpts, interpolation,
RK/AB integration and migration (host crystal router). After `migrate()` every
particle lives on the rank owning its element and `interpolator().data()`
holds that element and the reference coordinates in particle order on the
host — exactly what a projection kernel needs, and two small arrays to copy.

### What was added (`src/elgpu.okl`)

- `elgpuRHS`: Schiller–Naumann drag + gravity as the `lpm_t` user RHS.
- `elgpuForce`: the force on the fluid and the particle volume per particle,
  from the *current* fluid velocity at the new positions
  (`pointInterpolation_t::eval`; `lpm_t::interpolate()` would use the state
  extrapolated to the last RK stage).
- `elgpuDeposit`: transpose interpolation, `acc_n += q_p l_j(r) l_k(s) l_l(t)`
  over the nodes of the particle's element (`@atomic`), then `oogs` gather-
  scatter and division by the assembled mass matrix: the L2 projection of the
  point sources, conservative to round-off (`proj` = 5e-16 in
  `settling_gpu`). No ghost particles are needed: the projection is assembled
  across elements by the gather-scatter.
- Gaussian filter as explicit diffusion to pseudo-time `sigma^2 / 2` with
  `L f = invM dssum(K f)`, `K` the SEM stiffness matrix of nekRS's elliptic
  solvers (`core-weakLaplacianHex3D`); the weak form integrates to zero, so
  the filter conserves. The step is set from the spectral radius of `L`,
  measured once by power iteration (4.2e4 on the 6^3 N=5 box: 57 steps of
  3.6e-5 per field), and the sign from the Rayleigh quotient. Two things
  found on the way: `opSEM::laplacian` wraps this kernel but is unused in
  nekRS and passes an uninitialised field-offset table to it (segfault), and
  the obvious substitute, the weak divergence of the dssum-averaged strong
  gradient, has a spectral radius 3x smaller and leaves element-interface
  kinks in `phi_p` that the pressure solver cannot follow (divergence error
  8e-4 instead of 2e-5).
- The QTL divergence takes `d phi_f/dt = +div(phi_p v_p)` from the projected
  particle volume flux (deposited and filtered like the other fields) rather
  than a finite difference in time: the L2 projection's slope is
  discontinuous when a particle crosses an element face, and the finite
  difference turned that into node-scale noise (imposed divergence 6x too
  large). `Config::qtlFromFlux = false` keeps the finite difference.
- `elgpuSource`: `FU = F / (1 - phi_p)` and the QTL divergence on the device.
- Periodic wrap kernel; particles leaving a non-periodic domain are deleted.

### nekRS `lpm_t` finding

`lpm_t::integrateRK3` and `integrateRK4` (nekRS v26, `src/core/plugins/lpm.cpp`)
pass the start-of-step state `o_y_` to the user RHS at the intermediate
stages instead of the stage state `o_ytmp_` (the positions for the fluid
interpolation *are* the stage positions, only the state handed to the RHS is
stale). For a velocity-dependent RHS the scheme is then first order: in the
single-droplet check of `spray_jet2d` the velocity error was 1.9e-3 U0 with
"RK3" and 1.5e-5 U0 with RK2, which is implemented correctly.
`patches/nekRS/0001-lpm-rk-stage-state.patch` is the upstream fix (not
applied to the submodule or the installed nekRS). `elgpu` does not depend on
it: it integrates with its own SSP-RK3 through `lpm_t::setUserODESolver`
(`Coupling::ssprk3`, stage positions handed to the interpolator), which gives
7e-8 U0 on the same check, and accumulates the force on the fluid with the
same stage weights (prop `fAcc`, migrates with the parcel), so the fluid
receives what the parcels lose. `Config::rkOrder = 1, 2` fall back to
lpm_t's own integrators with the end-of-step force.

### Not there yet

Particle–particle contact (needs a device cell list: counting sort by cell,
then a neighbour loop; ghosts across ranks only for particles within the
contact distance of a rank's element boundary), IBM markers (the same deposit
kernel spreads marker forces; the body reduction is a segmented sum over
marker ids), and walls. The Gaussian filter would be cheaper as one implicit
Helmholtz solve with nekRS's elliptic solver than as ~60 explicit steps per
field (7 fields).

### Verification (`settling_gpu`)

The `settling_twoway` case on the device path (same box, particles, filter
width and `qtl`), 200 steps, 1 and 4 ranks identical to all printed digits.
Against the ppiclF path at t = 0.2: fluid momentum -5.666e-2 vs -5.631e-2
(0.6 %), particle momentum -4.865e-2 vs -4.851e-2 (0.3 %), drag sum -0.3417
vs -0.3410 (0.2 %); imposed divergence RMS 1.98e-4 vs 2.01e-4 and the
solver's divergence error 5.9e-5 vs 6.0e-5. Projection conservation is
round-off (3e-16) instead of the Gaussian truncation (4.6e-4). The global
momentum drift is -1.2e-3 of the gravity impulse (ppiclF path: -1.0e-3) with
the RK-weighted force on the fluid; with the end-of-step force it was -5.6e-3.

### Where the time goes

With ppiclF on the CPU the particle side is 87–97 % of the step
(`udfExecuteStep` in `settling_twoway` and `el_ibm_collision`). The
`lpm_t` path moves everything but migration and the two index arrays onto the
device; with OCCA's serial backend it is already 8–10x faster per step than
ppiclF for the same case (0.24 s vs 2.0 s on 1 rank, 0.09 s vs 0.9 s on 4),
and nothing in it is host-bound on a GPU except the migration. Only the CPU
backend has been run so far; the kernels use nothing beyond `@atomic`, which
nekRS reports as supported on the device.
