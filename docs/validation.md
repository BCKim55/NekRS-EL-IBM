# Verification and validation

Every check that has been run on this code, in one place. Each entry says
what is being tested and why, how the case is set up, what it is compared
against, how to run it, what came out, and what the result does and does not
show. The design notes (`design.md`) hold the derivations and the findings
along the way; this file is the evidence.

Terminology used below:

- **verification**: the code solves the equations it claims to solve
  (exact solutions, conservation budgets, agreement between two independent
  implementations of the same model);
- **validation**: the model reproduces experiments or established
  correlations (Johnson & Patel, Mordant & Pinton, Bickley, Pilch & Erdman).

All laptop runs: macOS arm64, nekRS v26.0 (`96b3cf9e`) built with gcc 16,
OCCA serial (CPU) backend, Open MPI, 1–8 ranks. GPU runs: Tuolumne (LLNL,
MI300A, HIP backend), see `design.md` §6. Figures are in `img/`, the
scripts that make them in `../scripts/`.

## Contents

| # | check | kind | case / script | status |
|---|---|---|---|---|
| V1 | one-way coupling: tracers in solid-body rotation | verification | `rotation_tracer` | round-off (2e-14) |
| V2 | two-way coupling budgets (projection, momentum, QTL divergence) | verification | `settling_twoway` | 1e-3 of the gravity impulse |
| V3 | device particle path vs ppiclF path | verification | `settling_gpu` | 0.2–0.6 % |
| V4 | spray parcels: single-droplet ODE, budgets, both paths | verification | `spray_jet2d`, `spray_jet2d_gpu` | 7e-8; 0.8 % |
| V5 | EL–IBM contact and momentum exchange | verification | `el_ibm_collision` | momentum closed to 0.3 % |
| V6 | ppiclF library before and after the patches | verification | `scripts/ppiclf_regression.py` | patched exact on 1–3 ranks |
| V7 | nekRS + ppiclF vs Nek5000 + original ppiclF | verification (cross-code) | `settling_twoway/nek5000_reference` | 0.2–1.2 % at t = 0.2 |
| V8 | TAB oscillator and breakup onset | verification + validation | `tab_breakup` | period 0.5 %; on Pilch & Erdman |
| V9 | plane jet from parcels vs Bickley | validation | `jet_bickley` | u_c within 3.5 % |
| V10 | fixed sphere Re 10 (NekIBM `simple` box), filter-width bracket | validation | `ibm_fixed_sphere`, `_fine` | 4.09 / 4.69 vs 4.15–4.3 |
| V11 | fixed sphere Re 100 vs Johnson & Patel, two resolutions | validation | `ibm_sphere_re100`, `hires` | C_d +18 % → +10 % |
| V12 | settling glass sphere vs Mordant & Pinton | validation | `ibm_settling_sphere`, `_fine`, `hires` | u_t −6 % (laptop); hires running |
| V13 | separate EL and IBM filter widths (ppiclF patch 0006) | verification | `el_ibm_collision`, `ibm_sphere_re100` | equal widths: previous results reproduced; split: budget −4.5e-3 |
| V14 | prescribed-motion body vs fixed body in a stream | verification | `ibm_sphere_moving` | C_d 1.257 ± 0.029 vs 1.289 |

Not done, and why, is listed at the end.

---

## V1. One-way coupling: tracers in a frozen solid-body rotation (`cases/rotation_tracer`)

**What is tested.** The plumbing between nekRS and ppiclF with nothing else
in the way: element layout handed to ppiclF as its overlap mesh,
interpolation of the fluid velocity at particle positions, particle
redistribution between ranks as particles cross element and rank
boundaries.

**Set-up.** nekRS's `ethier.re2` (a unit sphere, 32 elements), N = 7, the
velocity frozen at a solid-body rotation about z through the centre
(`u = -w y, v = w x`, w = 1), 1000 tracers at random positions, dt = 2π/600,
600 steps = one revolution, 1, 2, 4 and 8 ranks.

**Reference.** The velocity is linear, so GLL interpolation is exact. ppiclF
samples the fluid only at the first RK stage, so with this right-hand side
its RK3 is exactly forward Euler: after n steps a tracer is at its initial
position rotated by n·atan(w dt) and scaled by (1 + (w dt)²)^(n/2). This
discrete solution is what the coupling must reproduce; the difference from
the exact rotation (`errExact`) is ppiclF's own O(dt) time error and is only
reported.

**Run.**
```bash
scripts/run_case.sh cases/rotation_tracer rotation 8 --backend CPU
```

**Result.** `errDiscrete` (max over the 1000 tracers of the distance to the
discrete solution) is 1.9e-14 at step 550 and 2.1e-14 at step 600 on 8
ranks, the same on 1–4 ranks; `errExact` 1.7e-2 (= the Euler error of one
revolution at this dt). No tracer is lost.

**What it shows.** Interpolation, element ownership and migration are exact
to round-off, independently of the rank count. It does not test the force
projection or any two-way term.

## V2. Two-way coupling budgets (`cases/settling_twoway`)

**What is tested.** The three transfers from particles to fluid: force
projection (conservation), the momentum source (the fluid receives what it
is given), and the volume-fraction (QTL) divergence; plus global momentum
and rank-independence.

**Set-up.** Fully periodic unit box, 6 × 6 × 6 elements, N = 5, fluid
ρ = 1, ν = 0.01, 10 000 point particles of d_p = 0.01, ρ_p = 100 (particle
volume fraction 0.52 %, mass loading 0.52) at rest at random positions,
gravity g = −1 in z, Gaussian projection FWHM 0.15, `qtl = 1`, dt = 1e-3,
200 steps, 1 and 4 ranks.

**Reference.** Budgets normalised by the gravity impulse Σm_p |g| t:

- `proj` = ∫F_proj dV − Σf_p (normalised by |Σf_p|): the Gaussian is
  truncated at 1e-4 of its peak and integrated with GLL quadrature, so this
  is the projection's conservation error;
- `fluid` = P_f − P_f0 − ∫∫FU dV dt: the fluid momentum change equals the
  integrated source;
- `total` = P_f + P_p − P_0 − Σm_p g t: with periodic walls gravity on the
  particles is the only external force, so this is zero up to the
  1/(1 − φ_p) factor and the explicit lag between the particle and fluid
  updates;
- `div`, `divErr`: RMS of the imposed ∇·u = −(1/φ_f) Dφ_f/Dt and of the
  difference between it and the solved divergence.

**Run.**
```bash
scripts/run_case.sh cases/settling_twoway twoway 4 --backend CPU
```

**Result** at t = 0.2 (step 200):

| quantity | value |
|---|---|
| proj (projection conservation) | 3.6e-4 (Gaussian truncation + quadrature) |
| fluid (source received) | 1.6e-3 of the gravity impulse |
| total (global momentum) | −1.0e-3 of the gravity impulse |
| div, divErr | 2.0e-4, 6.0e-5 (nekRS's own divergence error on this grid) |
| 1 vs 4 ranks | identical to all printed digits |

**What it shows.** The projection conserves to the level set by its cutoff,
the fluid solver applies the source exactly, and the only momentum loss is
the explicit coupling lag at 1e-3 of the driving impulse. The QTL divergence
is imposed to the solver's own tolerance. It does not show that the drag law
or the filter width are physically right (see V7, V9).

## V3. Device particle path against the ppiclF path (`cases/settling_gpu`)

**What is tested.** `src/elgpu.hpp` (nekRS `lpm_t`, own deposit/filter/QTL
kernels, SSP-RK3) solves the same model as the ppiclF path. Two independent
implementations of binning, ghosting, projection (ppiclF's Gaussian sum vs an
L2 projection plus diffusion filter) and time integration.

**Set-up.** V2 on the device path, same box, particles, filter width,
`qtl = 1`, 200 steps, 1 and 4 ranks.

**Run.**
```bash
scripts/run_case.sh cases/settling_gpu twoway 4 --backend CPU
```

**Result** at t = 0.2:

| quantity | ppiclF path (V2) | device path | difference |
|---|---|---|---|
| fluid momentum P_f | −5.631e-2 | −5.666e-2 | 0.6 % |
| particle momentum P_p | −4.851e-2 | −4.865e-2 | 0.3 % |
| drag sum Σf | −0.3410 | −0.3417 | 0.2 % |
| imposed divergence RMS | 2.01e-4 | 1.98e-4 | |
| solver divergence error | 6.0e-5 | 5.9e-5 | |
| projection conservation | 4.6e-4 (truncation) | 3e-16 (round-off) | |
| global momentum drift | −1.0e-3 | −1.2e-3 | |
| 1 vs 4 ranks | identical | identical | |
| step time, 1 / 4 ranks | 2.0 s / 0.9 s | 0.24 s / 0.09 s | 8–10x |

**What it shows.** The two code paths agree to well under 1 % with no shared
code in the particle side, and both close their budgets. The remaining
difference is the filter (Gaussian sum vs diffusion) and the RK stage at
which the force is taken. The device path has only been run on the CPU
(OCCA serial) so far; its kernels use nothing beyond `@atomic`.

## V4. Spray parcels and injection, both paths (`cases/spray_jet2d`, `spray_jet2d_gpu`)

**What is tested.** Parcels of N_p droplets (per-droplet equation of motion,
per-parcel force and volume), run-time injection from a deterministic
generator shared by both paths, and the momentum budget of an open injection.

**Set-up.** Quasi-2D slot jet: box 2 × 1 × 0.1 (40 × 20 × 1 elements, N = 5),
periodic x and z, free-slip y (so the x momentum of gas plus parcels equals
the injected momentum exactly), gas ρ = 1, ν = 5e-4, slot width 0.02 at
x = 0.1, liquid ρ = 1000, d_p = 0.002 (τ_p = 0.44, Re_p0 = 4), u_jet = 1,
liquid volume fraction 0.01 in the slot (mass loading 10 at the nozzle),
10 parcels of 2.4 droplets per step, dt = 2e-3, 1000 steps, no gravity,
4 ranks.

**Reference.** (i) A single droplet released into quiescent gas against the
drag ODE integrated independently in Python; (ii) the injected momentum
(known exactly from the generator); (iii) the other path.

**Run.**
```bash
scripts/run_case.sh cases/spray_jet2d jet 4 --backend CPU
scripts/run_case.sh cases/spray_jet2d_gpu jet 4 --backend CPU
python3 scripts/jet_reference.py single <single-droplet logs>   # droplet vs drag ODE
python3 scripts/jet_reference.py jet cases/spray_jet2d cases/spray_jet2d_gpu   # budgets, figure
```

**Result.**

| | ppiclF path | device path |
|---|---|---|
| single droplet vs drag ODE, max error / U0 | 7e-8 | 7e-8 |
| momentum budget at t = 2 (max over the run) | +8.6e-4 (2.8e-3) | +9.4e-4 (2.8e-3) |
| fluid budget P_gas − ∫FU | 1.6e-4 | 1.6e-4 |
| P_gas, P_parcels at t = 2 | 1.0720e-2, 2.9314e-2 | 1.0806e-2, 2.9232e-2 |
| centreline gas u_x, max | 0.8333 | 0.8336 |
| centreline difference | — | max 0.0055, rms 0.0013 |
| step time | 1.25 s (10 k parcels) | 0.19 s |

![spray jet, both paths](img/spray_jet2d.png)

**What it shows.** Both paths integrate the droplet equation to RK3
accuracy, conserve momentum to 1e-3 of the injected amount, and agree with
each other to 0.8 % in momenta and 0.7 % in centreline velocity. A dt = 5e-3
run blew up at t = 0.5 when the induced gas velocity (0.83 u_jet) crossed
nekRS's advection CFL: that is the fluid solver's limit, not the coupling.

## V5. EL–IBM contact and momentum exchange (`cases/el_ibm_collision`)

**What is tested.** Soft-sphere contact between point particles and a
resolved (IBM) body, the body's rigid-body update including the contact
impulse, and that the three momentum reservoirs (fluid, EL, body) exchange
without loss.

**Set-up.** Fully periodic unit box, 8³ elements, N = 5, a free IBM sphere
of radius 0.15 (ρ_b/ρ_f = 4; the explicit body update needs
ρ_b/ρ_f > 1 + 3√(2π)σ/R = 2.28, see `design.md` §1) at rest at the centre;
a cylindrical cloud of 200 EL particles (d_p = 0.01, ρ_p = 2000) fired at it
at 0.5 along +x; contact stiffness k_n = 500, restitution 0.9, 10 particle
sub-steps per fluid step (contact time ≈ 4.5e-3 against dt = 5e-3);
Gaussian FWHM 0.06, marker spacing 0.025; 240 steps to t = 1.2; 4 ranks.

**Reference.** No external force acts, so P_fluid + P_EL + P_body must stay
at the initial EL momentum P_0. Geometric checks: `gap/dp` (closest approach
of any particle to the surface in particle diameters; a small negative value
is the soft-contact overlap, a large one would be a particle inside the
sphere), `contacts`, `bounced` (particles with u_x < 0).

**Run.**
```bash
scripts/run_case.sh cases/el_ibm_collision collision 4 --backend CPU
python3 scripts/plot_collision.py cases/el_ibm_collision          # figures
python3 scripts/snapshots_to_vtk.py cases/el_ibm_collision        # ParaView
```

**Result** at t = 1.2:

| quantity | value |
|---|---|
| body velocity | 0 → 0.271 (the cloud's momentum, minus what the fluid keeps) |
| momentum split P_f / P_EL / P_b | 3.74e-2 / 5.55e-2 / 1.15e-2 |
| (P_f + P_EL + P_b − P_0)/P_0 | −3.2e-3 |
| min gap/d_p | 0.109 (no particle inside the sphere) |
| bounced particles | 4 of 200 |
| IBM marker slip, RMS | 4.8e-3 ≈ 2 % of the body speed |

![collision summary](img/collision_summary.png)

**What it shows.** The contact keeps every particle outside the body, the
impulse reaches the body, and the total momentum closes to 0.3 %. The loss
is the explicit lag of the coupling at a stiff contact; it does not change
with the marker spacing. The sphere does not deform the far field much at
this loading, so this is a verification of the exchange, not of the drag.

## V6. ppiclF library, pristine vs patched (`scripts/ppiclf_regression.py`)

**What is tested.** `patches/ppiclF/0001–0005` (RK3 bound after bin
migration, overlap-mesh node test, projection exit, periodic ghost
de-duplication, filter-vs-period abort) do not change results on inputs that
do not trigger them, and fix the ones that do.

**Set-up.** Both libraries built the same way (gfortran, no `-DPPICLC`),
ppiclF's own `examples/stokes_2d` (250 particles per rank, Stokes drag
τ_p = 1/9.8, gravity 9.8, from rest, dt = 1e-4, 10 000 steps) on 1, 2 and 3
ranks, and ppiclF's `tests/test1` (interpolation and projection diagnostics).

**Reference.** The exact solution v(t) = −gτ(1 − e^(−t/τ)) and its integral,
checked for every particle (matched between VTU outputs by ppiclF tag; the
VTU `TIME` stamp is the start of the step, the state is one dt later).

**Result.**

| library | ranks | max error in v | max error in displacement | particles off by > 1e-5 |
|---|---|---|---|---|
| pristine | 1 | 3.7e-8 | 3.4e-8 | 0 |
| pristine | 2, 3 | 1.8e-3 | 1.7e-3 | 1 of 500 / 750 |
| patched | 1, 2, 3 | 3.7e-8 | 3.4e-8 | 0 |

`test1` prints identical diagnostics for both libraries on 1–3 ranks (TimeI
error 9.3e-9, Prjct error 1.3e-8).

**What it shows.** The pristine library loses one particle's RK stage after
a migration on 2 or more ranks (patch 0001); the patched one is exact on
every rank count, and the other patches are inert where not triggered.

## V7. nekRS + ppiclF against Nek5000 + original ppiclF (`cases/settling_twoway/nek5000_reference`)

**What is tested.** The whole two-way coupling against the code ppiclF was
written for: the same particle library (unpatched, one rank count where it
is exact) in its own Nek5000 example wiring.

**Set-up.** V2's box on Nek5000 v19 with the unpatched ppiclF, the Nek5000
example wiring (momentum source F/(1 − φ_p)/ρ from `GetProFldIJKEF`,
interpolation and `IntegrateParticle` in `userchk`), BDF2, same mesh, dt,
filter and particle population (different random positions), 4 ranks each.
Inputs, README and log are in the case directory.

**Result** (P_f fluid momentum, P_p particle momentum, Σf total drag):

| t | code | P_f | P_p | Σf |
|---|---|---|---|---|
| 0.04 | Nek5000 + ppiclF | −5.682e-3 | −1.583e-2 | −0.2306 |
| 0.04 | nekRS + ppiclF | −5.333e-3 | −1.552e-2 | −0.2283 |
| 0.12 | Nek5000 + ppiclF | −2.974e-2 | −3.389e-2 | −0.3297 |
| 0.12 | nekRS + ppiclF | −2.919e-2 | −3.362e-2 | −0.3300 |
| 0.20 | Nek5000 + ppiclF | −5.696e-2 | −4.883e-2 | −0.3404 |
| 0.20 | nekRS + ppiclF | −5.631e-2 | −4.851e-2 | −0.3410 |

![validation against the originals](img/validation.png)

**What it shows.** At t = 0.2 the two codes agree to 1.2 % (P_f), 0.7 %
(P_p) and 0.2 % (Σf). The larger early difference (6 % in P_f at t = 0.04,
when P_f is small) comes from the time integrators (BDF2 vs BDF3/EXT with
the source handling of `design.md` §1) and the different random positions.
Nek5000's own momentum budget closes to 1 % of the gravity impulse, nekRS's
to 0.1 %. The device path (V3) sits inside the same band.

## V8. TAB deformation and breakup (`cases/tab_breakup`)

**What is tested.** The TAB droplet model on the device path: the distortion
ODE and its SSP-RK3 integration, the dynamic drag, the breakup event, and
the breakup time against experiment.

**Set-up.** Single droplets released into quiescent gas at We = 10, 12, 15,
18, 25, 35, 50, 80, 120, 200, 300 (d = 0.002, σ = 1e-4, ρ_l = 1000,
μ_l = 1e-3, Oh = 0.07, gas ρ = 1, ν = 5e-5), one-way, dt = 2e-3, 600 steps.

**Reference.** (i) The first three steps of y, y' and v against an
independent Python SSP-RK3 of the same ODE; (ii) the damped-oscillator
period 2π/√(ω² − ζ²) with ω² = 8σ/(ρ_l r³), ζ = 5μ_l/(2ρ_l r²); (iii) TAB's
critical We = 12; (iv) Pilch & Erdman (1987) breakup-initiation time
T_ini = 1.9 (We − 12)^(−1/4) in T = t_b U √(ρ_g/ρ_l)/d.

**Run.**
```bash
scripts/run_case.sh cases/tab_breakup jet 1 --backend CPU | tee build/tab.log
python3 scripts/tab_reference.py build/tab.log --out docs/img/tab_breakup.png
```

**Result.**

| check | result |
|---|---|
| first three steps vs Python SSP-RK3 | 6 digits |
| We = 10 oscillation period | 0.224 / 0.222 / 0.224 vs analytic 0.2230 (0.5 %) |
| first peak / y_eq (y_eq = We/24) | 1.61 (damped analytic 1.76; We decays through the drag) |
| breakup at We = 10, 12, 15 | none (damping moves the threshold above 12) |
| T at We = 18, 25, 35, 50, 80, 120 | 1.17, 1.06, 1.00, 1.00, 0.95, 0.93 |
| Pilch & Erdman T_ini, same We | 1.21, 1.00, 0.87, 0.77, 0.66, 0.59 |
| child size d_32/d | 0.40 (We 18) → 0.10 (We 120) |

![TAB validation](img/tab_breakup.png)

**What it shows.** The oscillator is integrated correctly and TAB's y = 1
event sits on the measured initiation curve at low We, 1.6x above it at
We = 120 (TAB saturates where the correlation keeps falling); the total
breakup time is 4–5x longer, as expected for an initiation criterion. **Open
problem:** at We = 200 and 300 the children break again within a few steps,
their τ_p falls below dt and the explicit drag goes unstable; breakup needs
either a dt that resolves the children or an implicit drag update.

## V9. Plane laminar jet from parcels against Bickley (`cases/jet_bickley`)

**What is tested.** That a dilute stream of small droplets hands its momentum
to the gas correctly: far downstream the gas must be the Bickley plane
laminar jet of the injected kinematic momentum flux J, whatever the near
field looks like. The nekRS-only control is the same J as a fixed Gaussian
body force.

**Set-up.** Open slab, inflow with a 5 % coflow, outflow, free-slip sides,
40 × 32 × 1 elements, N = 5, ν = 5e-4; slot width 0.04 at x = 0.5,
u_jet = 10, d_p = 0.0005 (τ_p = 0.016, stopping length ≈ 0.16), liquid
volume fraction 5e-6 in the slot (mass loading 0.01, so the parcels keep only
0.5 % of J after they slow down), 10 parcels per step, dt = 2.5e-3, 8000
steps to t = 20 (steady), 11 700 parcels in the domain at the end.

**Reference.** Bickley: u_c = 0.4543 (J²/(νx))^(1/3),
u = u_c sech²(0.2752 (J/(ν²x²))^(1/3) y), evaluated with the momentum flux
measured in the run (0.0173 for the body force, 14 % below the nominal 0.02
because the free-slip sides confine the entrainment: the coflow outside the
jet is slowed and that deficit is in the integral).

**Run.**
```bash
scripts/run_case.sh cases/jet_bickley jet 4 --backend CPU                 # parcels
# source_only = 1 in jet.par for the body-force reference
python3 scripts/jet_bickley.py build/bick/parcels build/bick/source --xmin 0.9 --xmax 1.8
```

**Result** in the window 0.9 < x < 1.8:

| | parcels | body force |
|---|---|---|
| ∫u(u − u_co)dy over the window | 0.0169, constant to 0.3 % | 0.0173, constant to 0.3 % |
| J from the u_c decay (u_c⁻³ linear in x) | 0.0168 | 0.0165 |
| virtual origin x0 | 0.572 (nozzle 0.5 + stopping length) | 0.325 |
| u_c vs Bickley at the measured J, max error | 3.5 % | 3.7 % |
| half-width vs Bickley, max error | 13 % | 12 % |
| profile at x = 1.2 | sech² | sech² |

Sensitivity (to t = 15): filter width 0.04 → 0.06 and 10 → 40 parcels per
step leave J identical to three digits and the decay constants within 1–2 %;
the filter only moves the virtual origin (0.57 → 0.51).

![plane jet vs Bickley](img/bickley.png)

**What it shows.** The parcel jet carries 2.3 % less than the body-force
jet (0.5 % stays with the parcels), and both reproduce the x^(−1/3)
centreline decay, the sech² profile and a momentum flux constant along the
jet. The half-width error is the confinement (both runs share it and differ
from each other by 6 %), not the coupling. A first attempt at mass loading
10 put half of J in the parcels at the end of the domain: not an error, the
liquid keeps its mass flux times the gas speed, which is why the test runs
dilute.

## V10. Fixed sphere at Re 10 in NekIBM's `simple` box (`cases/ibm_fixed_sphere`, `_fine`)

**What is tested.** The IBM direct forcing on a resolved body: total force,
symmetry, and the effect of the Gaussian interface width.

**Set-up.** NekIBM's `IBM_test/simple`: sphere d = 0.01 at the centre of a
(4d)³ box, inflow u = 1 and outflow in y, periodic x and z, ν = 1e-3
(Re = 10), 20³ elements, N = 5 (GLL spacing ≈ d/25), 1963 markers at
d/25, dt = 2.5e-4, 200 steps (U t/d = 5), 4 ranks. `_fine` halves the marker
spacing and the FWHM (0.04 d, one GLL spacing; base 0.08 d, two spacings).
C_d = F_y/(ρ u² π d²/8) from the summed marker force.

**Reference.** Unbounded sphere at Re 10: Schiller–Naumann 4.15,
Clift–Grace–Weber standard drag 4.30. The box is not unbounded: the sphere
is 2d from the inflow and in a laterally periodic array of spacing 4d
(volume fraction 0.8 %), for which the Stokes-limit array correction
(Sangani–Acrivos) is a factor 1.5 and the finite-Re correction smaller.

**Run.**
```bash
scripts/run_case.sh cases/ibm_fixed_sphere sphere 4 --backend CPU | tee build/c_re10.log
cases/ibm_fixed_sphere_fine/setup.sh && scripts/run_case.sh cases/ibm_fixed_sphere_fine sphere 4 --backend CPU | tee build/c_re10f.log
```

**Result** at U t/d = 5:

| FWHM | C_d | marker slip / u_in | lateral force |
|---|---|---|---|
| 0.08 d (two GLL spacings) | 4.694 | 18 % | 1e-10 |
| 0.04 d (one GLL spacing) | 4.090 | 31 % | 1e-10 |
| Tuolumne, 4 MI300A, FWHM 0.08 d | 4.6937 (bit-identical) | | |

![midplane slices](img/ibm_slices.png)

**What it shows.** The two widths bracket the unbounded references: a wide
Gaussian enlarges the effective body (the stagnant region in the slice
extends outside the circle), a one-cell one under-enforces no-slip (slip
31 %). The excess of the wide case over the correlations is of the size of
the box blockage. Converging the drag needs a finer grid, not a narrower
filter (V11 shows that). Lateral forces vanish as the symmetry demands.
The direct comparison with NekIBM's own run of this box was not possible
(NekIBM does not run on macOS, see the end of this file).

## V11. Fixed sphere at Re 100 vs Johnson & Patel (`cases/ibm_sphere_re100`, `hires`)

**What is tested.** The IBM drag and wake against a classic reference, and
their convergence with resolution.

**Set-up.** NekIBM's `fixsphere/Re100` box: 8d × 15d × 8d, d = 0.01, inflow
u = 1 and outflow in y, periodic x and z, ν = 1e-4 (Re = 100).

| | laptop | hires (Tuolumne) |
|---|---|---|
| elements | 16 × 30 × 16 | 32 × 60 × 32 (NekIBM's mesh) |
| N (mean GLL spacing) | 5 (d/10) | 7 (d/28, edge ≈ d/60) |
| markers | d/17 | d/50 (7854 markers) |
| Gaussian FWHM | 0.12 d | 0.04 d |
| dt, steps | 2e-4, 500 | 5e-5, 2000 |
| ranks | 4 (CPU) | 8 MI300A, 1.1 s/step, 41 min |

Both run to U t/d = 10. `src/fieldDump.hpp` writes the z = 0 plane at the
last step for the slices.

**Reference.** Johnson & Patel (1999): C_d = 1.09, standing eddy length
0.87 d (steady state).

**Run.**
```bash
scripts/run_case.sh cases/ibm_sphere_re100 sphere 4 --backend CPU | tee build/c_re100.log
# Tuolumne, 2 nodes: REPO_ROOT=$PWD scripts/tuolumne_batch.sh cases/ibm_sphere_re100/hires sphere 8
python3 scripts/settling_reference.py --settle build/c_settle.log --re100 build/c_re100.log \
    --re100-hires cases/ibm_sphere_re100/hires/tuolumne_8gpu.log --re10 build/c_re10.log --out docs/img/nekibm_refs.png
```

**Result** at U t/d = 10:

| | laptop | hires | Johnson & Patel |
|---|---|---|---|
| C_d | 1.286 (+18 %) | 1.201 (+10 %) | 1.09 |
| standing eddy length (centreline u < 0) | 1.00 d | 0.95 d | 0.87 d |
| marker slip / u_in | 3.9 % | 2.8 % | |
| dC_d/d(U t/d) at the end | | −0.002 | 0 |

![NekIBM reference cases](img/nekibm_refs.png)

**What it shows.** Halving the interface width and doubling the resolution
halves the drag error and moves the eddy length toward the reference: the
Gaussian interface width is the leading error and the scheme converges
toward Johnson & Patel as it is refined. Part of the remaining 10 % is the
wake still developing at U t/d = 10 (C_d still falling slowly; the
reference is the steady state), part the interface. The lateral period of 8d
(volume fraction 0.5 %) adds a small blockage.

## V12. Settling glass sphere vs Mordant & Pinton (`cases/ibm_settling_sphere`, `_fine`, `hires`)

**What is tested.** A free IBM body under gravity: hydrodynamic force on a
moving body, the rigid-body update with buoyancy, and the terminal velocity
against experiment.

**Set-up.** NekIBM's `Set1sphere` input: Mordant & Pinton's 0.5 mm glass
bead in water nondimensionalised with d = 1/6, ν = 0.005416 (Ga = 49.1),
ρ_p/ρ_f = 2.56, g = 9.81, released from rest.

| | laptop | `_fine` | hires (Tuolumne) |
|---|---|---|---|
| domain | 1.2 × 3 × 1.2 (7.2d × 18d × 7.2d), periodic | same | 1.2 × 10 × 1.2 (60d column) |
| elements, N | 10 × 25 × 10, 7 (d/10) | N = 8 | 21 × 180 × 21, 7 (d/20; NekIBM's mesh) |
| markers, FWHM | d/17, 0.12 d | d/17, 0.10 d | d/20 (NekIBM n_dh = 20), 0.10 d |
| dt, steps | 2e-3, 500 (t = 1) | same | 1e-3, 2000 (t = 2) |

**Reference.** Mordant & Pinton (2000): Re_t = 41, u_t = 1.332 in these
units. Also the standard-drag ODE (Schiller–Naumann + added mass,
`scripts/settling_reference.py`), which gives 1.427 (Re 44), 7 % above the
measurement.

**Run.**
```bash
scripts/run_case.sh cases/ibm_settling_sphere settle 4 --backend CPU | tee build/c_settle.log
# Tuolumne: flux batch -q pbatch -N 2 -B <bank> -t 3h -n 8 -g 1 --output=build/settle_hires.out \
#     scripts/tuolumne_batch.sh cases/ibm_settling_sphere/hires settle 8
```

**Result** at t = 1:

| | v | Re | v / u_t (M&P) | force / net weight |
|---|---|---|---|---|
| laptop, N = 7 | 1.249 | 38.4 | 0.937 | 98 % |
| `_fine`, N = 8, FWHM 0.10 d | 1.244 | 38.3 | 0.933 | |
| hires, 60d column | pdebug run reached t = 0.79 (4.2 s/step on 8 GPUs); full run queued in pbatch | | | |
| Mordant & Pinton | 1.332 | 41 | 1 | |

**What it shows.** The sphere settles 6 % slower than measured and the
result is insensitive to the filter width (0.4 %), so the deficit is not the
interface. The 7.2d lateral period of NekIBM's box is a periodic array of
spheres (φ^(1/3) = 0.083), which hinders settling by up to 14 % in the
Stokes limit (Hasimoto) and less at Re 41, i.e. the right size and sign. The
hires case has the same cross-section, so it tests the marker and grid
resolution, not the column width; a wider periodic box would be the next
check.

## V13. Separate EL and IBM filter widths (`el_ibm_collision`, `ibm_sphere_re100`)

**What is tested.** ppiclF patch 0006 (per-particle Gaussian width in the
projection, `PPICLF_R_JFILT`) and `Config::filterWidthIbm` change nothing
when both widths are equal, and keep the budgets when they differ.

**Set-up.** V5 with `filter_width_ibm = 0` (equal) and `0.03` (IBM half the
EL width 0.06); V11's laptop case (markers only).

**Result.**

| run | body u at t = 1.2 | (P_f + P_EL + P_b − P_0)/P_0 | marker slip (end) | note |
|---|---|---|---|---|
| equal widths, before patch | 0.27100 | −3.155e-3 | 4.843e-3 | |
| equal widths, after patch | 0.27100 | −3.155e-3 | 4.844e-3 | all other digits identical |
| EL 0.06, IBM 0.03 | 0.27176 | −4.542e-3 | 6.670e-3 | slip during impact 1.2e-2 vs 2.0e-2 |
| `ibm_sphere_re100`, C_d at U t/d = 10 | 1.2860 before and after | | | bit-identical |

**What it shows.** The patch is inert where the widths agree. With the
narrower IBM shell the momentum budget closes to 0.45 % instead of 0.32 %
(same explicit-lag mechanism at the contact), the body threshold
ρ_b/ρ_f drops from 2.28 to 1.64, and the slip halves. See `design.md` §7.1.

## V14. Prescribed-motion body against the fixed body in a stream (`ibm_sphere_moving`)

**What is tested.** `Body::prescribed` (design doc §7.2): the markers
enforce a given x(t), v(t) and the reported fluid force is the same as for
the equivalent fixed body in a stream.

**Set-up.** The `ibm_sphere_re100` box (16 × 30 × 16, N = 5) made fully
periodic, fluid at rest, sphere driven at −u_in in y from y = +5d, 500 steps
to U t/d = 10, same markers and FWHM as the fixed case, 4 ranks.

**Reference.** `ibm_sphere_re100` on the same grid (Galilean transform; the
only differences are the inflow/outflow boundaries and the start-up).

**Run.**
```bash
scripts/run_case.sh cases/ibm_sphere_moving moving 4 --backend CPU | tee build/moving_re100.log
```

**Result** over U t/d = 8–10 (11 reports):

| | C_d mean | C_d std | marker slip |
|---|---|---|---|
| moving sphere, fluid at rest | 1.257 | 0.029 | 3.7–4.0 % |
| fixed sphere in a stream | 1.289 | 0.001 | 3.9 % |

**What it shows.** The prescribed body gives the same drag within 2.4 %,
with a ±2 % step-to-step oscillation the fixed case lacks: the markers of a
moving body cross GLL cells and the explicit direct forcing samples new
grid velocities every step. It averages out and shrinks with resolution.
The last report (U t/d = 10) drops 5 % because the sphere reaches its own
start-up disturbance through the 15d periodic column, so the window stops
at 9.8.

---

## Not done, and why

- **NekIBM itself on the same inputs.** Its Nek5000 v17 core builds with
  gfortran 16 on macOS arm64 but segfaults in the dealiasing interpolation
  (`set_dealias_rx` → `intp_rstd` → `mxm`) on the first step, at `-O0` and
  with a 512 MB stack; it targets Intel compilers on Linux. The IBM routines
  themselves initialise and print their checks. Recipe for a Linux/Intel
  node: `IBM_test/simple`, compare `QP*.dat` column `jf0` with V10.
- **NekIBM's paper values** (Yang & Balachandar 2021) are behind a paywall;
  the comparison is to the experiments and correlations NekIBM itself used.
- **GPU correctness of the device particle path (V3)**: run only on the OCCA
  serial backend. The IBM/ppiclF path has run on MI300A (V10, V11) and
  reproduces the CPU result bit for bit for V10.
- **Round particle-laden jet** (`jet_modarress`, Modarress 1984): the full
  mesh does not link on macOS arm64 (nek interface data > 2 GB per rank);
  coarse smoke test only. Needs the Tuolumne install.
- **TAB cascade at We ≥ 200** (V8): explicit drag unstable once the
  children's τ_p < dt.
- **A dt-refinement study** for the coupling itself: only the two dt that
  bracket nekRS's CFL limit were run (V4, V9).
- **Wider settling column** to isolate the hindrance in V12.
