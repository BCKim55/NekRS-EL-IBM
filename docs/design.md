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

## 2b'. Plane laminar jet from parcels (`cases/jet_bickley`)

Open slab (inflow with a 5 % coflow, outflow, free-slip sides, 40 x 32 x 1
elements, N = 5, `nu = 5e-4`), a slot stream of small droplets
(`tau_p = 0.016`) that stop within ~0.16 of the nozzle, so far downstream
the gas must be the Bickley plane laminar jet of the injected kinematic
momentum flux J = 0.02:
`u_c = 0.4543 (J^2/(nu x))^(1/3)`, `u = u_c sech^2(0.2752 (J/(nu^2 x^2))^(1/3) y)`.
The nekRS-only reference is the same J as a fixed Gaussian body force
(`source_only = 1`). Both run to t = 20 (steady).

Body force: `int u (u - u_co) dy` is constant along the jet to 0.3 % at
0.0173, 14 % below J because the free-slip sides confine the entrainment
(the coflow outside the jet is slowed, the excess velocity goes negative at
|eta| > 2 and that deficit is in the integral). Against Bickley with the
measured 0.0173: centreline decay within 3.7 %, half-width within 12 % and the
profile is sech^2 — the plane-jet similarity is reproduced; the 14 % is the
domain, not the solver.

First parcel attempt (u_jet = 1, mass loading 10): the gas carried only
0.0105 and the parcels the rest. Not an error: once the droplets have
decelerated to the gas speed (~0.5 u_jet) they still carry their mass flux
times that speed, half of J at this loading. The test therefore needs the
liquid *mass* flux to be small while J is kept: u_jet = 10, mass loading
0.01 (parcels keep 0.5 % of J). The stopping stream is a stronger, more
concentrated force than the Gaussian blob (J / (w u_jet tau_p) against
J / (2 pi sigma^2)), so the slot was widened to 0.04 and the gas still
peaked at 0.76 u_jet-independent units and crossed nekRS's CFL limit at
dt = 5e-3; the parcel runs use dt = 2.5e-3.

Parcel run (u_jet = 10, w = 0.04, phi_inj = 5e-6, 10 parcels per step, 8000
steps to t = 20, 11 700 parcels in the domain at the end), window
0.9 < x < 1.8:

| | parcels | body force |
|---|---|---|
| `int u (u - u_co) dy` over the window (set J = 0.02) | 0.0169, constant to 0.3 % | 0.0173, constant to 0.3 % |
| J from the u_c decay (`u_c^-3` linear in x) | 0.0168 | 0.0165 |
| virtual origin x0 | 0.572 (nozzle 0.5 + stopping length) | 0.325 |
| u_c vs Bickley at the measured J, max rel. err | 3.5 % | 3.7 % |
| half-width vs Bickley, max rel. err | 13 % | 12 % |
| profile at x = 1.2 | sech^2 (figure) | sech^2 |

The parcel jet carries 2.3 % less than the body-force jet (0.5 % stays with
the parcels, the rest is the slightly different near field); both reproduce
the x^-1/3 centreline decay, the sech^2 profile and a momentum flux that is
constant along the jet. The half-width error of both against Bickley is the
coflow/confinement (the excess velocity is negative outside the jet, which
narrows the half-width definition), not the particle coupling: the two
runs differ from each other by 6 % in b.

Step-2 sensitivity (same case to t = 15, window 0.9 < x < 1.8, J measured
at 0.0173 as above):

| run | J(x) mean | J from u_c decay | x0 | u_c vs Bickley | half-width vs Bickley |
|---|---|---|---|---|---|
| base: filter 0.04, 10 parcels/step | 0.0169 | 0.0168 | 0.572 | 3.5 % | 13.4 % |
| filter 0.06 | 0.0169 | 0.0166 | 0.511 | 4.2 % | 15.3 % |
| 40 parcels/step (47 000 parcels) | 0.0169 | 0.0169 | 0.563 | 2.5 % | 13.4 % |

The far-field jet does not depend on the filter width or the parcel count:
the momentum flux is identical to three digits and the decay constants agree
within 1–2 %; the filter only moves the virtual origin (a wider filter spreads
the source, x0 0.57 → 0.51). Filter width 0.03 at dt = 2.5e-3 crossed the
CFL limit (a narrower filter makes the stopping stream a stronger local
force); it needs dt = 1.25e-3 and was not run. A dt check was not run
either, beyond the two dt that bracket the CFL limit.

## 2b''. Droplet-driven jet vs a water jet (`cases/spray_jet_water`)

Question: how close does a gas jet driven by water-like droplets come to a
water jet? Quasi-2D slab 0 < x < 3, |y| < 0.8 (free-slip sides), z 0.1
periodic, 1920 elements N = 5, ν = 5e-4, coflow u_co through the inflow
plane, outflow at x = 3. Droplets: ρ_p/ρ_f = 833, d = 2e-3 (Stokes stopping
length ≈ 40 slot widths, as a 100 µm water drop at 10 m/s from a 2 mm slot),
injected at x = 0.5 from a slot w = 0.04 at u_jet = 10 with liquid volume
fraction φ in the stream, i.e. kinematic momentum flux J = (ρ_p/ρ_f) φ u_jet²
w. Sweep φ = 5e-6, 5e-5, 5e-4 (J 0.017, 0.17, 1.7; liquid/gas mass flux
0.004, 0.04, 0.42). References, same mesh: (i) `inlet_jet = 1`, the same J as
a top-hat velocity inlet of width w on x = 0, no droplets; (ii) Bickley with
x0 and J fitted to (i); (iii) the liquid jet, u = u_jet, all momentum in the
liquid, half-width w/2. Fields are time averages (`t_avg`); the droplet
velocity is the parcel mean in x bins and the liquid momentum flux J_l =
(ρ_p/ρ_f) φ u_jet w ⟨u_d⟩ (mass flux conserved). `sweep.sh`,
`scripts/spray_jet_water.py`.

| φ | u_co | u_d at x = 2 (u_jet 10) | liquid J at outflow | gas u_c / inlet jet u_c at x = 2 | min u (backflow) |
|---|---|---|---|---|---|
| 5e-6 | 0.05 | 1.9 | 4 % | 1.57 | 0.005 |
| 5e-5 | 0.11 | 2.5 | 8 % | 0.49 | −0.48 |
| 5e-4 | 0.25 | 4.1 | 23 % | 0.55 | −1.6 (inlet jet too) |

Findings. Dilute: the droplets give the gas 96 % of their momentum inside
the box, so the result is a spray-induced gas jet, not a water jet; because
the momentum is deposited along the path the induced jet is narrower than
the inlet jet (half-width 0.03–0.08 vs 0.05–0.18) and faster on the
centreline beyond x ≈ 1.3, while the inlet jet follows Bickley within 4 %.
Dense: the droplets keep 0.41 u_jet at x = 2 and 0.23 at the outflow, the
liquid holds 23 % of J at the outflow and the gas jet stays about the slot
width, i.e. "liquid core plus a thin gas sheath", but the gas carries more
than half of J from x ≈ 1.5 on. How water-jet-like the droplet jet is is set
by (1) whether the distance of interest is shorter than the droplet stopping
length and (2) the loading; with this d only x < 1 (12 slot widths) behaves
as a water jet, and the liquid core itself is outside the EL model.

Caveats. The droplet jet adds momentum without gas, so the coflow must
supply its entrainment: at u_co = 0.05 the φ = 5e-5 case reversed the outer
flow to −0.7 and its gas momentum flux read 111 % of J (backflow squared and
the pressure term), hence u_co ∝ J^(1/3). At J = 1.7 even the inlet jet
(Re_slot 516, unsteady) recirculates in the 1.6-wide box, Bickley is not
drawn, and the average (t 3–5) is short and asymmetric; a wider box and a
longer average need the GPU (0.15 s/step on 4 MI300A, 0.6 on 4 CPU ranks).
The dense case needs dt = 2e-4 with `advectionSubCyclingSteps = 2` (CFL
abort at 5e-4 without subcycling). The device path has no restart: the
parcel state is not checkpointed, so a run that overruns its allocation
leaves no fields.

## 2c. Round particle-laden jet (`cases/jet_modarress`, prepared)

The standard experiment for two-way EL codes with solid particles: Modarress,
Tan & Elghobashi (1984), air jet, D = 20 mm, Re_D = 13 300, 50 um glass
beads (rho_p/rho_g = 2490, St = tau_p U_0/D = 11.5), mass loading 0.32 and
0.85; centreline gas and particle velocity decay and radial profiles. The
case is set up in nondimensional form (D = U_0 = rho_g = 1): box 20D x 8D
x 8D, 40 x 24 x 24 elements at N = 5 (10^6 points, a coarse LES grid),
smoothed top-hat nozzle inflow with a 3 % coflow, outflow, symmetry sides,
parcels over the nozzle disc at the gas speed. It is turbulent, so it is a
GPU-node job (hours), and on this laptop the full mesh does not even link
(Nek5000's static arrays exceed the macOS arm64 linker's 2 GB range at any
rank count). A 16 x 10 x 10 copy (`jet_coarse.re2`) ran 20 steps on 4 CPU
ranks (0.7 s/step) as a smoke test: injection, coupling and the report are
exercised, nothing physical is claimed. The experimental curves still have to
be digitised from the paper (or from Apte, Mahesh & Lundgren 2003, who
simulated it) into `scripts/`.

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

### Second `lpm_t` finding: variable numbering

`addVariable` / `addProp` store `dofIds[name] = dofIds.size()` (the variable's
position in the list) while `getDOF` / `getProp` slice at `id * fieldOffset`
(a field index). The two agree only while every earlier variable has one
component, so a variable registered after a multi-component one is laid over
the wrong fields: a 2-component "tab" registered after the 3-component
"velocity" landed on velocity's y and z. `patches/nekRS/0002-lpm-dof-field-
index.patch` makes the id the first-field index. This repo works around it by
registering single-component variables first and the one multi-component
variable (velocity, fAcc) last.

### TAB droplet deformation and breakup (`cases/tab_breakup`)

`Config::tab` adds the TAB distortion y and its rate as particle DOFs
(O'Rourke & Amsden 1987: CF = 1/3, CB = 1/2, CK = 8, CD = 5), integrated
with the particle SSP-RK3, the dynamic drag `C_d (1 + 2.632 y)` of Liu,
Mather & Reitz, and in-place breakup when y > 1: the parcel becomes
`N_p (r/r_32)^3` droplets of the Sauter radius from the TAB energy balance,
y is reset, the liquid volume kept, no velocity kick. Verification on
droplets released into quiescent gas at We = 10–300 (`d = 0.002`,
`sigma = 1e-4`, `rho_l = 1000`, `mu_l = 1e-3`, one-way):

- the first three steps of y, y' and v match an independent Python SSP-RK3
  of the same ODE to 6 digits;
- the We = 10 droplet oscillates with period 0.224/0.222/0.224 against the
  analytic damped period 0.2230, about `y_eq = We/24`, first peak 1.61 y_eq
  (damped analytic 1.76; We decays through the drag);
- no breakup at We = 10, 12, 15 (TAB's critical We is 12 for an undamped
  droplet; liquid damping, Oh = 0.07, moves it above 15); breakup at We ≥ 18
  with `T = t_b U sqrt(rho_g/rho_l)/d` = 1.17, 1.06, 1.00, 1.00, 0.95, 0.93 for
  We = 18, 25, 35, 50, 80, 120 against the Pilch & Erdman (1987) *initiation*
  time `1.9 (We-12)^-1/4` = 1.21, 1.00, 0.87, 0.77, 0.66, 0.59: TAB's y = 1
  event is the onset of breakup and sits on that curve at low We, 1.6x above
  it at We = 120 (TAB saturates where the correlation keeps falling); the
  *total* breakup time is 4–5x longer, as expected;
- child size `d_32/d` = 0.40 (We = 18) to 0.10 (We = 120).
- The We = 200 and 300 droplets broke repeatedly within a few steps and
  were lost: each breakup leaves children whose own We is still above 12,
  and after three cascades `tau_p` of the children falls below dt, the
  explicit drag goes unstable and the parcel is deleted as unfound. After
  breakup dt must resolve the children's response time, or the drag needs
  an implicit (exponential) update — not done yet.

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

## 5. Validation against the original ppiclF and NekIBM

(`validation.md` collects every check in the repo, these included, in a
uniform set-up / reference / result / interpretation form.)

Three checks that the merged code reproduces the codes it was built from:
the ppiclF library itself before and after `patches/ppiclF`, the ppiclF
two-way coupling against ppiclF's own host (Nek5000), and the IBM against
NekIBM. Scripts and inputs: `scripts/ppiclf_regression.py`,
`cases/settling_twoway/nek5000_reference`, `cases/ibm_fixed_sphere`;
figure: `scripts/plot_validation.py` → `docs/img/validation.png`.

![validation against the originals](img/validation.png)

### 5a. ppiclF library, pristine vs patched (ppiclF's own examples)

Both libraries built the same way (gfortran, no `-DPPICLC`) and run on
ppiclF's `examples/stokes_2d` (250 particles per rank, Stokes drag,
tau_p = 1/9.8, from rest, dt = 1e-4, 10 000 steps) on 1, 2 and 3 ranks, and
on ppiclF's `tests/test1` (interpolation and projection diagnostics).
`stokes_2d` has an exact solution, v = -g tau (1 - e^{-t/tau}), so every
particle is checked (the VTU `TIME` tag is the start of the step; the state
in the file is one dt later).

| library | ranks | max error in v | max error in displacement | particles off |
|---|---|---|---|---|
| pristine | 1 | 3.7e-8 | 3.4e-8 | 0 |
| pristine | 2, 3 | 1.8e-3 | 1.7e-3 | 1 of 500 / 750 |
| patched | 1, 2, 3 | 3.7e-8 | 3.4e-8 | 0 |

The pristine library is exact on one rank and loses one particle's RK3
stage after a bin migration on 2 and 3 ranks (patch 0001, documented in
section 2); the patched library is exact on every rank count. `test1`
prints identical diagnostics for both libraries on 1–3 ranks (TimeI error
9.3e-9, Prjct error 1.3e-8), i.e. patches 0002–0005 change nothing on
inputs they are not triggered by.

### 5b. Two-way coupling: nekRS + ppiclF vs Nek5000 + original ppiclF

`cases/settling_twoway/nek5000_reference` runs the settling box on Nek5000
v19 with the unpatched ppiclF and its Nek5000 example wiring (momentum
source F / (1 - phi_p) / rho from `GetProFldIJKEF`, interpolation and
`IntegrateParticle` in `userchk`), BDF2, same mesh, dt, filter and
particle population (different random positions). Budget at the same times
(4 ranks each; P_f fluid momentum, P_p particle momentum, sum f total drag):

| t | code | P_f | P_p | sum f |
|---|---|---|---|---|
| 0.04 | Nek5000 + ppiclF | -5.682e-3 | -1.583e-2 | -0.2306 |
| 0.04 | nekRS + ppiclF (this repo) | -5.333e-3 | -1.552e-2 | -0.2283 |
| 0.12 | Nek5000 + ppiclF | -2.974e-2 | -3.389e-2 | -0.3297 |
| 0.12 | nekRS + ppiclF (this repo) | -2.919e-2 | -3.362e-2 | -0.3300 |
| 0.20 | Nek5000 + ppiclF | -5.696e-2 | -4.883e-2 | -0.3404 |
| 0.20 | nekRS + ppiclF (this repo) | -5.631e-2 | -4.851e-2 | -0.3410 |

At t = 0.2 the two codes agree to 1.2 % in fluid momentum, 0.7 % in
particle momentum and 0.2 % in total drag; the early-time differences are
larger (6 % in P_f at t = 0.04, when P_f is itself small) and come from the
different time integrators (BDF2 in Nek5000 against BDF3/EXT in nekRS,
with the source extrapolated as described in section 1) and the different
random positions. The device path (`settling_gpu`) sits within the same
band (section 4). Nek5000's own momentum budget closes to 1 % of the
gravity impulse (`total` column), nekRS's to 0.1 %.

### 5c. IBM: fixed sphere drag at Re = 10 (`cases/ibm_fixed_sphere`)

NekIBM's `IBM_test/simple` is a fixed sphere (d = 0.01, u_in = 1, nu = 1e-3,
Re = 10) in a 4d x 4d x 4d box, periodic laterally, inflow/outflow in y.
`cases/ibm_fixed_sphere` is the same set-up on this repo's IBM
(`elibm::Body` with the new `fixed` flag, 1963 markers at 4e-4 spacing,
Gaussian FWHM 8e-4, 20^3 elements, N = 5, dt = 2.5e-4, 200 steps). The
reported drag coefficient C_d = F_y / (rho u^2 pi d^2 / 8) from the total
marker force:

| step | t | C_d |
|---|---|---|
| 5 | 0.00125 | 11.89 |
| 20 | 0.005 | 6.11 |
| 50 | 0.0125 | 4.89 |
| 100 | 0.025 | 4.73 |
| 200 | 0.05 | 4.694 |

Reference values for an unbounded sphere at Re = 10: Schiller–Naumann
4.15, Clift–Grace–Weber standard drag 4.3. The 4d box is not unbounded: the
sphere sits 2d from the inflow plane and in a laterally periodic array of
spacing 4d (volume fraction 0.8 %), for which the Stokes-limit array
correction (Sangani–Acrivos) is a factor 1.5 and the finite-Re correction
is smaller. The 9–13 % excess over the unbounded correlations is therefore
the blockage of NekIBM's test box, not an IBM error; the marker slip
(mean |u - U_body| over the markers, 18 % of u_in) is the diffuse-interface
velocity left on a one-cell-wide Gaussian interface and is the same
quantity NekIBM reports as its interface error. Lateral forces are 1e-10,
as the symmetry demands.

Running NekIBM itself on this machine did not work: the author's build
targets Intel compilers on a cluster (`lelg` 72 000, `lpmin` 8). With
gfortran 16 on macOS arm64 it builds (`-std=legacy
-fallow-argument-mismatch`, `lpmin` 8, short `SESSION.NAME`) but the
author's `simple` mesh fails Nek's Jacobian check and a regenerated 10^3
mesh segfaults in Nek's dealiasing interpolation
(`set_dealias_rx` → `intp_rstd` → `mxm`) on the first step, after the IBM
initialisation has completed and printed its volume-fraction check; the
crash is in the Nek5000 v17 core, not in the IBM routines, and persists at
`-O0` and with a 512 MB main stack. The direct NekIBM-vs-this-repo drag
comparison therefore still needs a Linux/Intel node; the recipe above is
what to run there (`IBM_test/simple`, `QP*.dat` column `jf0`).

### 5d. IBM against NekIBM's own reference data

NekIBM's paper (Yang & Balachandar, J. Supercomput. 77, 2897, 2021) is
behind a paywall and its figures could not be read; its reference list and
`IBM_test` inputs show what it validated against: Johnson & Patel (1999) for
the fixed sphere, Mordant & Pinton (2000) for a settling sphere, and the
Richardson–Zaki / Yin–Koch correlations for clouds. The two single-sphere
references were run on this repo's IBM with NekIBM's own inputs
(`cases/ibm_sphere_re100`, `cases/ibm_settling_sphere`, refinements in
`*_fine`), at the resolution a laptop allows (NekIBM's Set1sphere mesh has
79 380 elements at N = 7); figure `scripts/settling_reference.py` →
`docs/img/nekibm_refs.png`.

![NekIBM reference cases](img/nekibm_refs.png)

Midplane velocity magnitude of the same runs at their last step
(`src/fieldDump.hpp` writes `plane_z0.csv`, `scripts/plot_slices.py` plots
it; red circle = nominal sphere). With the 0.08 d Gaussian the stagnant
region at Re 10 extends well outside the circle, with 0.04 d it hugs it;
at Re 100 the standing eddy behind the sphere is about one diameter long
(Johnson & Patel: 0.87 d); the settling sphere is shown in its own frame
with the wake behind it.

![midplane slices](img/ibm_slices.png)

| case | NekIBM input | reference | this repo | notes |
|---|---|---|---|---|
| settling glass sphere, Ga 49.1, rho_p/rho_f 2.56 | `Set1sphere` (nu 0.005416, d 1/6, g 9.81, 1.2 x 1.2 cross-section periodic) | Mordant & Pinton: Re_t = 41, u_t = 1.332 | u = 1.249 at t = 1 (N = 7, FWHM 0.12 d, laptop), 1.244 (N = 8, FWHM 0.10 d); on NekIBM's mesh (`hires`: 21 x 180 x 21, N = 7, markers d/20, 60 d column, Tuolumne 16 MI300A) u = 1.273 at t = 1 and 1.299 at t = 1.2 (Re 40.0), still rising slowly | −6 % → −2.5 % with resolution; filter-insensitive (0.4 %); the 7.2 d lateral period (same in `hires`) hinders settling by up to 14 % in the Stokes limit (Hasimoto, phi^(1/3) = 0.083), less at Re 41; standard-drag ODE gives 1.427 (Re 44), 7 % above the measurement |
| fixed sphere, Re 100 | `fixsphere/Re100` box (8d x 15d x 8d, inflow/outflow in y, periodic x, z) | Johnson & Patel: C_d = 1.09, standing eddy 0.87 d | C_d = 1.286 at U t/d = 10 (N = 5, FWHM 0.12 d, 16 x 30 x 16 elements, laptop), eddy 1.00 d; C_d = 1.201, eddy 0.95 d on NekIBM's mesh (`hires`: 32 x 60 x 32, N = 7, FWHM 0.04 d, Tuolumne, 8 MI300A) | +18 % → +10 % with the sharper interface (marker slip 3.9 % → 2.8 %); C_d still falling by 0.002 per d/U at U t/d = 10, so part of the rest is the developing wake |
| fixed sphere, Re 10 | `simple` box (4d)^3 | standard drag 4.15–4.3 | C_d = 4.69 with FWHM 0.08 d (two GLL spacings, marker slip 18 %), 4.09 with FWHM 0.04 d (one spacing, slip 31 %) | the two widths bracket the reference: a wide Gaussian enlarges the body, a one-cell one under-enforces no-slip; converging needs a finer grid, not a narrower filter |

What this shows: with NekIBM's own set-ups the IBM reproduces the reference
drag and settling velocity to 6–18 % on the laptop grids, with the sign and
size expected from the Gaussian interface width at d/10 resolution (fixed
spheres) and from the periodic-array hindrance of NekIBM's narrow column
(settling). On NekIBM's own meshes (section 6, Tuolumne) the Re 100 drag
error halves to +10 % and the standing eddy shortens from 1.00 d to 0.95 d
(reference 0.87 d), and the settling sphere comes within 2.5 % of the
measured terminal velocity (1.299 at t = 1.2, still rising, against
1.332), i.e. the interface width is the leading error and the scheme
converges towards the references as it is refined; the remaining
offset at U t/d = 10 includes a still-developing wake (Johnson & Patel's
value is the steady state). The paper's numbers were not available, so the
comparison is to the experiments and correlations NekIBM itself used.

Practical findings on the way: nekRS's mesh repair (`xyz repair`) diverges
and the Jacobian check aborts when a case's `.re2`/`.usr` are symlinks, and
on these periodic genbox meshes some 8-rank partitions failed the same way
where 4 ranks pass; a Re 100 run at N = 6 with dt = 2.5e-4 reached CFL > 1
on the element-edge GLL spacing and blew up at step 120 (dt = 2e-4 at N = 5
is stable); the nek interface data must stay under 2 GB per rank on macOS
arm64 (lelt x lx1^3 below about 6.5e5).

## 6. Running on Tuolumne (LLNL, MI300A, flux)

nekRS has no official recipe for Tuolumne; the Frontier `PrgEnv-gnu` one
(nekRS_HPCsupport) works with the gfx942 target. The ppiclF build script uses
gfortran flags and the UDF links `libgfortran`, so everything is built with
the GNU programming environment:

```bash
# no cpe/25.03: mixing it with the system default cray-mpich links libmpifort 9.x against libmpi 8.x
module reset && module load PrgEnv-gnu gcc-native craype-accel-amd-gfx942 \
  rocm/6.3.1 cmake python/3.12.2 && module unload cray-libsci
export NEKRS_HOME=/p/lustre5/<user>/nekrs-install
cd external/nekRS && CC=cc CXX=CC FC=ftn ./build.sh -DCMAKE_INSTALL_PREFIX=$NEKRS_HOME \
  -DENABLE_HYPRE_GPU=off -DENABLE_ADIOS=off -DOCCA_ENABLE_CUDA=off -DOCCA_ENABLE_DPCPP=off  # ADIOS links the python module's old libzmq
# cases (ppiclF with the Cray wrappers, one rank per GPU)
export FC=ftn CC=cc LAUNCHER="flux run -g 1 -n"   # this flux takes -g, not --gpus-per-task
flux alloc -q pdebug -N 1 -B <bank> -t 60m
scripts/run_case.sh cases/ibm_fixed_sphere sphere 4 --backend HIP --device-id 0  # flux shows one GPU per task
# batch (from the repo root; flux copies the script to /var/tmp, so it finds the repo through $PWD)
flux batch -q pbatch -N 1 -B <bank> -t 3h -n 4 -g 1 --output=build/x.out scripts/tuolumne_batch.sh <case dir> <casename> 4
```

Measured on one node (4 MI300A): `ibm_fixed_sphere` 0.27 s/step (11x the
laptop, C_d bit-identical to the laptop run), `ibm_sphere_re100/hires`
1.9 s/step on 4 GPUs and 1.1 s/step on 8 (two nodes; the host-staged ppiclF
side dominates at 21M points). `pbatch` carried a backlog of several
thousand jobs, so the Re 100 hires run went through `pdebug` on two nodes
(2000 steps in 41 min, `cases/ibm_sphere_re100/hires/tuolumne_8gpu.log`, the
`sphere:` history plus every 50th step time); the settling hires case waits
in `pbatch` via `scripts/tuolumne_batch.sh`; it was finally run in `pdebug` on
4 nodes (16 GPUs, 2.5 s/step: 1.7x faster than 8 GPUs; 32 ranks segfaulted
in the set-up, not pursued) for 1207 steps to t = 1.2
(`cases/ibm_settling_sphere/hires/tuolumne_16gpu.log`). ppiclF is pre-built on a login
node (gslib download) into `build/<case path with / as _>`; `| head` on a
running job's output file closes the pipe and kills the job.

NekIBM-resolution inputs for the section 5d comparison are in
`cases/ibm_sphere_re100/hires` (32 x 60 x 32 elements, N = 7, 21M points,
markers d/50, FWHM 0.04 d) and `cases/ibm_settling_sphere/hires`
(21 x 180 x 21 elements on the 60 d column, N = 7, 28M points, markers d/20
as NekIBM's `n_dh = 20`, FWHM 0.1 d). The `.re2` meshes are committed
(13 and 17 MB) because genbox is not available there.

## 7. Pilobolus sporangium launch (preparation)

Target: the discharge of a *Pilobolus kleinii* sporangium (Yafetto et al.
2008, PLoS ONE 3:e3237): a 350 µm sporangium is pushed off the vesicle by a
sap jet to ~9 m/s within 0.2 ms, then decelerates at ~1300 m/s² during
0.3–1.2 ms of flight, 7–9x what an equivalent sphere's drag gives, while
the sap spreads into a disc around it and a tail that breaks into ~100 µm
drops. Re ≈ 200. The five steps of the hand-off, in order; each has its
check.

### 7.1 Separate EL and IBM filter widths

ppiclF has one Gaussian filter, and its width sets the ghost-exchange and
truncation distances. The EL projection (force, φ_p) wants a width of a
few particle spacings; the IBM spreading wants one to two grid spacings.
Patch `0006-per-particle-gaussian-width` gives every ppiclF particle its
own FWHM in an rprop slot (`PPICLF_R_JFILT`, registered through the new
`ppiclf_solve_InitParticleFilter(jfilt)`): the projection loop uses that
particle's σ for both the exponent and the normalisation (real and ghost
particles, grid and sub-bin projections), and the global filter is set to
the larger of the two widths so every neighbour is still reached. A
per-particle width wider than the global one aborts. `Config::filterWidthIbm`
(0: same as `filterWidth`) is what the markers get; EL particles get
`filterWidth`, including injected parcels; the body-stability check and the
marker thickness use the IBM width.

Checks: with equal widths `el_ibm_collision` reproduces its previous result
to all printed digits but the last of the slip (4.843e-3 → 4.844e-3) and
`ibm_sphere_re100` bit for bit (C_d 1.2860). With EL 0.06 and IBM 0.03
(`filter_width_ibm` in the collision par) the run is stable (the body
threshold ρ_b/ρ_f > 1.64 instead of 2.28), the marker slip halves
(4.8e-3 → 6.7e-3 of the body speed at the end, 2.0e-2 → 1.2e-2 during the
impact), the sphere ends at u = 0.2718 instead of 0.2710, and the global
momentum closes to −4.5e-3 of P_0 against −3.2e-3 with equal widths: the
same explicit-lag loss, slightly larger because the narrower shell reacts
faster at the contact.

### 7.2 Prescribed-motion bodies

`Body::prescribed`: the UDF sets x, v (and w) for the new time before
`step()`; the coupling pushes that state to ppiclF before the particle
sub-steps so the markers force the fluid to the new velocity, skips the
dynamics for that body (periodic wrap only) and still reports the marker
force sum in `hydroForceOnBodies()`. For an accelerating body the force on
the body proper is that sum plus ρ_f V dU/dt, the inertia of the fluid
inside the shell (Uhlmann); the Pilobolus UDF adds it.

Check (`cases/ibm_sphere_moving`): the `ibm_sphere_re100` box made fully
periodic, fluid at rest, the sphere driven at −u_in in y from y = +5d. In
the body frame that is the fixed sphere in a stream. Over U t/d = 8–10 the
moving sphere gives C_d = 1.257 ± 0.029 (mean ± std over 11 reports) against
1.289 for the fixed sphere on the same grid, i.e. 2.4 % lower with a
±2 % step-to-step oscillation the fixed case does not have: the markers of
a moving body cross GLL cells, and the explicit direct forcing samples the
grid velocity at new relative positions every step. The oscillation
averages out; a finer grid or a wider IBM Gaussian reduces it. The
periodic 15d column also lets the sphere approach its own start-up
disturbance after 10d of travel (the last report drops 5 %), so the
comparison window is U t/d < 10. The run itself was slow (25 s/step on the
laptop against 2.5 s/step fixed, with one step of 2000 s): the machine
slept during the run; the moving-body cost on Tuolumne (settling sphere,
4.2 s/step on 8 GPUs against 1.1 s/step fixed) is a separate, real item.

### 7.3 Fixed sphere at Re 300 (prepared, needs Tuolumne)

`cases/ibm_sphere_re300`: the Re 100 hires mesh and interface (32 x 60 x 32,
N = 7, markers d/50, FWHM 0.04 d), ν = 1/300. Johnson & Patel (1999):
C_d = 0.656 (time mean), St = f d/U = 0.137 (hairpin shedding, period
U t/d = 7.3). The wake needs U t/d ≈ 50 to settle and ~5 periods to
measure St, so 20 000 steps at dt = 5e-5: `numSteps = 5000` chunks with
`startFrom`, 1.1 s/step on 8 MI300A → 1.5 h per chunk. Acceptance: C_d within
the Re 100 hires error (+10 %) and St within 5 %.

```bash
cases/ibm_sphere_re300/setup.sh
flux batch -q pbatch -N 2 -B <bank> -t 3h -n 8 -g 1 --output=build/re300_1.out \
    scripts/tuolumne_batch.sh cases/ibm_sphere_re300 sphere 8
# then startFrom = "sphere0.f00001" in sphere.par and resubmit for chunks 2-4
```

On pdebug (2 nodes, 60 min) use `cases/ibm_sphere_re300/chunk.sh <k> [steps]`
inside each allocation instead: 2000 steps for chunk 1 (JIT), 2500 after, 8
chunks; it sets `startFrom`, keeps each chunk's last checkpoint as
`chunk_<k>.fld` and writes the force history to `build/re300_chunk_<k>.out`.

### 7.4 Single-phase launch (`cases/pilobolus_launch`)

SI units. Air at rest (ρ 1.2, μ 1.8e-5), a box x ∈ [−2, 10] mm
(hires: [−3, 13]), y, z ∈ [−5, 5] mm, far-field u = 0 on five faces and an
outflow at x+; a fine band |y|, |z| < 0.6 mm (hires 0.44) along the flight
path with 0.2 mm elements at N = 5 (40 µm ≈ d/9; hires 0.08 mm at N = 7,
11 µm ≈ d/30), graded outside; 11 760 elements (hires 105 800, 54M points).
The sporangium is an equivalent sphere d = 0.35 mm, ρ_b = 1150, driven from
x = 0.25 mm by the measured history: a half-cosine rise to v_max = 8.6 m/s
in t_rise = 0.22 ms (peak acceleration 6e4 m/s²), then v = v_max −
1300 (t − t_rise); or, with `trajectory = file.csv` (t_us, speed), the
smoothed measurement. The UDF reports x, v, a, the marker force sum F_m, the
force on the body F = F_m + ρ_f V a, the drag coefficient and the
deceleration F/(ρ_b V) that force alone would give; `free_flight = 1`
releases the body at t_rise and reports the deceleration of a dry sphere.
The difference between the measured 1300 m/s² and the dry deceleration is
the sap's contribution, the target of step 5. dt from the CFL on the band's
edge GLL spacing: 1 µs (hires 0.3 µs). The laptop grid needs 8 ranks on
macOS (nek interface < 2 GB per rank). Added mass is irrelevant at
ρ_b/ρ_f ≈ 960; the wake history effect is what the simulation resolves.

Laptop result (prescribed, 1200 steps, 8 ranks, 3.4–6.8 s/step): during
the coasting phase the fluid force on the sporangium is 2.4–3.3 µN,
C_d = 0.785 against Schiller–Naumann 0.832 at Re 185 (6 % low with a marker
slip of 20 % of the body speed on this d/9 grid), and the deceleration that
force alone gives the sporangium is 111 m/s² (93–141). The measured
1300 m/s² is 12x that: **the air accounts for 9 % of the observed
deceleration, the sap must supply the remaining 1190 m/s²**. During the
acceleration phase the marker sum reaches −5.3 µN, of which 2.1 µN is the
inertia of the fluid inside the shell (ρ_f V a at 6e4 m/s²); the force on
the body proper peaks at −4.2 µN. The force carries a sawtooth of ±0.2 µN
with a period of 27 µs = the time to cross one 0.2 mm element at 7.5 m/s:
the grid-crossing noise of 7.2 on a coarse grid (figure
`scripts/pilobolus_launch.py` → `docs/img/pilobolus_launch.png`).

Free flight (`free_flight = 1`, released at t_rise with 8.60 m/s): the dry
sporangium loses only 0.12 m/s in the next 0.96 ms, a mean deceleration of
125 m/s² (C_d 0.70–0.75, 10 % below Schiller–Naumann on this grid), against
the measured 1300 m/s². Both runs agree on the answer to the hand-off's
question: without the sap the sporangium would coast almost undecelerated.

![Pilobolus launch](img/pilobolus_launch.png)

### 7.5 Sap droplets and adhesion (`cases/pilobolus_launch/sap`)

The sap leaves the nozzle as a jet that breaks into 100–150 µm drops with
15–20 µm satellites (We ≈ 0.1 in air: no aerodynamic breakup, so TAB and
KH-RT stay off; τ_p ≈ 40 ms ≫ dt, so the explicit drag is fine). The case
injects that population directly: for t < t_inj rank 0 emits, every step,
the main drops and satellites whose volume adds up to sap_volume · dt/t_inj
(fractional counts carried over), uniformly over the nozzle disc at the
origin, at v_jet along x with a Gaussian lateral spread; `Coupling::inject`
bins them. EL projection width filter_width (2.5e-4, a few drop spacings),
IBM width filter_width_ibm (5e-5) — the split of 7.1.

Adhesion (`Config::attachOnContact`): a droplet sticks to the first body it
touches. In ppiclF's RHS the first stage that finds a contact force calls
`elibm_Attach`: the parcel stores the body index and its offset in the body
frame (`PPICLF_R_JATT`, `JAX..JAZ`), takes the body's velocity at once, and
the momentum it loses, m_p (u_b − v_p), goes to the body through the same
impulse channel as the contact impulse (`JIX..JIZ`). From then on the
parcel's RHS is the rigid velocity at its attachment point (no drag, no
contact; it is inside the no-slip shell), `PlaceMarkers` snaps it to the
rigid position after each step like a marker, its volume still projects to
φ_p, and its mass is added to the body's inertia and weight in the body
update (`attachedMass()`), so a free body decelerates as the loaded
sporangium would. A tether (ligament tension between tail drops and the
body) is the next model if the stuck mass alone does not explain the
measured deceleration; it is not implemented.

Check (`cases/el_ibm_adhesion`): the collision demo with `attach = 1`.
107 of the 200 particles stick (the rest miss the sphere), none bounces or
enters it (closest approach −0.02 d_p, the overlap at which contact is
first detected), the sphere ends at u = 0.276 against 0.271 with elastic
contact (the stuck mass, 2.6x the sphere's excess mass, rides along), and
the total momentum closes to −3.2e-3 of P_0, the same as without
attachment. A first version lost 18 %: free particles pushing stuck ones
(EL–EL contact) got their force but the stuck ones, moving rigidly, passed
no reaction to the body; the stuck particle now feeds that reaction into
the body's impulse channel and uses its rigid velocity in the dashpot.
With `attach = 0` the collision case is unchanged to all digits.

Laptop result (prescribed v(t), 1200 steps, 8 ranks, 3–23 s/step with the
EL work): of the 1011 droplets injected in 0.3 ms (30 main drops, ~980
satellites) 775 are stuck by 0.86 ms and 1005 are still in the box at the
end; the stuck mass is 23.2 ng = 77 % of the sap and **90 % of the
sporangium's own mass** (the 30 nL of sap outweighs the 26 ng sporangium),
reached in a saturating curve: the main drops, which carry the mass, hit
within 0.4 ms, the satellites trickle in until 0.9 ms. The fluid force on
the sporangium during coasting is 1.8 µN against 2.9 µN dry (−36 %,
C_d 0.4–0.6 against 0.79): the sap jet at 12 m/s drags the air behind the
sporangium along and lowers its relative velocity. With the stuck mass
counted as inertia the air alone would decelerate the loaded projectile by
~40 m/s², 3 % of the measured 1300 m/s², less than the dry sphere's 9 %.
Sticking droplets therefore cannot explain the measured deceleration
through their mass or their drag; the momentum they bring is forward
(impacts at v_jet > v_b push the sporangium). What remains is the tension
of the sap ligament/tail between the launched mass and the sporangiophore
(the tether of the paragraph above), or a larger effective cross-section of
the sap disc than the equivalent sphere; both are open. Figure: bottom row
of `docs/img/pilobolus_launch.png` (droplet counts, stuck mass, sap
momentum; the sap curves in the upper panels are the loaded sporangium).

### 7.6 Particle stream at a held sphere (`cases/el_ibm_stream`)

The continuous-injection counterpart of `el_ibm_collision`, as a test bed
for EL-IBM and EL-EL contact before the liquid-bridge force: box 0 < x < 2,
y and z periodic (0..1), inflow u = 0.5 at x = 0, outflow at x = 2,
ν = 0.005 (Re_D 30), IBM sphere D = 0.3 held at (0.8, 0.5, 0.5), one EL
particle per step (d_p = 0.01 = D/30, ρ_p 2000, St 3.7) born at x = 0.1 with
the stream velocity in the mid-plane sheet, |y − 0.5| < 0.12 (birth
positions 1.5 d_p apart). Soft-sphere contact k_n 500, e 0.9, 10 sub-steps,
`attach = 0`. ppiclF removes the particles that leave through the outflow.
Reported every step on rank 0 (`gatherElRows`, `src/elSnapshot.hpp`):
particles in contact, distinct particles ever in contact and ever moving
upstream (the robust collision count, since only a collision sends one
upstream), sphere force and C_d. `scripts/plot_stream.py` (summary, gif),
`scripts/snapshots_to_vtk.py` + `scripts/paraview_el_movie.py --view xy
--color-by contact --state` for ParaView.

Result (2400 steps, dt 5e-3, 0.49 s/step on 4 CPU ranks): 2400 injected,
1445 collided (60 %), 1265 left, ~1130 in the box at steady state from
t ≈ 7; no particle inside the sphere at any snapshot; C_d 4.7 → 2.9 while
the particle-laden wake and the momentum the particles give the gas
develop. A contact lasts about one step, so the end-of-step contact flag
misses some; the upstream count is the one to use. After the collisions the
sheet spreads over the whole periodic z (median |z − z_c| 0.18 at t = 12):
EL-EL contacts and the 3D flow round the sphere scatter the particles, so
the x-y view is a projection. A planar option (zero z force on the
particles) would be the next step for a truly 2D picture.

### 7.7 Cylinder bodies and the quasi-2D droplet stream (`cases/el_ibm_cylinder`)

`Body::length > 0` makes a body a cylinder with its axis along z spanning
that length (the periodic slab thickness): volume πR²L, markers in rings on
the lateral surface (no end caps, `nt x nz`, rings staggered by half a
spacing), contact on the in-plane distance (`elibm_ContactForce` gets an
`icyl` flag, passed through `ppiclc_elibm_SetBodies`), translation in x, y
and spin about z only, I = ½ m R², and the IBM feedback bound uses 2 instead
of 3 (shell volume 2πRLt over πR²L). Everything else (markers, impulses,
attachment) is shared with the sphere.

Case, SI units: slab 30 x 10 x 2 mm (24 x 8 x 1 elements, N 5), air
(ρ 1.2, μ 1.8e-5), inflow 0.5 m/s, outflow, y and z periodic; cylinder
D = 3 mm, ρ 1374, free, at (8, 5) mm; water droplets 100 µm every third
step (3333/s) at x = 2 mm, |y − y_c| < 1.2 mm, mid-plane; k_n 1 N/m, e 0.9,
20 sub-steps; dt 1e-4 with 2 advection sub-cycles (CFL 0.67 at the
cylinder), 0.3 s. The slab must be thicker than the ppiclF ghost cutoff
(~ filter width 0.6 mm; 0.5 mm was refused), hence 2 mm. No gravity (the
plane is horizontal).

Result: 1000 droplets injected, 619 sent upstream by a collision, 756 left,
~250 in the slab; none inside the cylinder; droplets stay within 0.3 mm of
the mid-plane (median), unlike the 3D sphere case. The cylinder accelerates
downstream at ~0.15 m/s² (droplet impulses + drag; mass 1.9e-5 kg), reaching
0.043 m/s and x = 14.7 mm at 0.3 s, and spins weakly (|ω| < 0.02 rad/s).
C_d (on u_in, D, L) 2.2 falling to 2.0 as the cylinder moves with the
stream; the 30 % blockage of the periodic y images raises it above the
unbounded value (~1.4 at Re 100). `scripts/plot_stream.py --scale 1000
--xlim 0 30 --ylim 0 10 --body ...`.

