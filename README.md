# NekRS-EL-IBM

Coupled Euler–Lagrange (EL) point particles and immersed-boundary (IBM) rigid
particles on top of [nekRS](https://github.com/Nek5000/nekRS), using
[ppiclF](https://github.com/dpzwick/ppiclF) for the Lagrangian side and
[NekIBM](https://github.com/YunchaoYang/NekIBM) as the reference for the IBM
direct forcing.

![EL particles fired at a free IBM sphere](docs/img/collision.gif)

*Left: x–y projection, where particles above and below the sphere overlap its
outline. Right: x against the distance from the axis through the sphere
centre, where the sphere is an exact half disc; no particle enters it (closest
centre 0.47 d_p from the surface, touching at 0.5 d_p).*

- EL particles: drag, gravity, two-way momentum coupling (Gaussian projection,
  `F / (1 - phi_p)`), volume-filtered continuity
  `div u = -(1/phi_f) D phi_f/Dt`, soft-sphere EL–EL contact.
- IBM spheres: surface markers in ppiclF with explicit direct forcing, rigid
  6-DOF bodies in C++, soft-sphere EL–IBM contact.
- Two particle back ends. `src/elibm.hpp` runs ppiclF on the CPU (EL + IBM +
  contact, host-staged; nekRS's GPU backends should work unchanged but only
  the CPU backend has been tested). `src/elgpu.hpp` is a prototype that keeps
  EL particles on the device with nekRS's `lpm_t` and adds the projection and
  filter as OCCA kernels (EL + two-way + QTL, no contact or IBM yet); see
  [docs/design.md](docs/design.md) §4.

## Layout

| Path | What |
|---|---|
| `external/nekRS` | nekRS v26.0 (`96b3cf9e`), submodule |
| `external/ppiclF` | ppiclF `72483db` (last upstream commit), submodule |
| `external/NekIBM` | NekIBM `9631f97`, submodule; reference only, not built |
| `patches/ppiclF/` | fixes applied to ppiclF at build time (the submodule stays pristine) |
| `patches/nekRS/` | upstream fix for `lpm_t`'s RK3/RK4 stage state (for reference; not needed by this repo) |
| `src/elibm.hpp` | `elibm::Coupling`: EL + IBM coupling for a nekRS UDF (ppiclF, CPU) |
| `src/elgpu.hpp`, `src/elgpu.okl` | `elgpu::Coupling`: EL particles on the device with `lpm_t` (prototype) |
| `src/ppiclf_user_elibm.f`, `src/PPICLF_ELIBM.h`, `src/ELIBM.inc` | ppiclF user routines and particle layout for EL + IBM |
| `src/ppiclfCoupler.hpp`, `src/ppiclf_cbind.f` | low-level nekRS ↔ ppiclF helpers and C bindings |
| `scripts/build_ppiclf.sh` | builds a case-specific `libppiclF.a` and `udf.env` |
| `scripts/run_case.sh` | builds ppiclF for a case and runs it with nekRS |
| `scripts/plot_collision.py` | figures for the collision demo |
| `cases/` | verification cases and the collision demo |
| `docs/design.md` | model, numerics, ppiclF findings, verification, GPU roadmap |

## Build and run

nekRS must be installed from `external/nekRS` (or the same commit) with double
`dfloat`. The UDF must be compiled with the same C++ compiler and standard
library as nekRS; with Open MPI wrappers set `OMPI_CXX` to the compiler nekRS
was built with (a libc++ UDF loaded into a libstdc++ nekRS segfaults in
`dlopen`).

```bash
git clone --recursive https://github.com/BCKim55/NekRS-EL-IBM.git
cd NekRS-EL-IBM
export NEKRS_HOME=$HOME/.local/nekrs
export OMPI_CXX=g++-16            # whatever built nekRS
scripts/run_case.sh cases/el_ibm_collision collision 4 --backend CPU
python3 scripts/plot_collision.py cases/el_ibm_collision
```

`run_case.sh` rebuilds `libppiclF.a` only when the case's `user_routines`,
the shared model in `src/`, or the patches change, and always relinks the UDF.

## Cases

| Case | What it checks | Result |
|---|---|---|
| `rotation_tracer` | one-way: tracers in a frozen solid-body rotation | positions match ppiclF's discrete solution to 2e-14 on 1–8 ranks |
| `settling_twoway` | two-way: 10k particles settling in a periodic box, with `qtl` | projection error = Gaussian truncation (3.6e-4); momentum to 1e-3 of the gravity impulse; QTL imposed to nekRS's own divergence error; 1 and 4 ranks identical |
| `settling_gpu` | `settling_twoway` on the device path (`src/elgpu.hpp`) | matches the ppiclF path to 0.2–0.6 % in momenta and drag, same divergence error; 1 and 4 ranks identical; 8–10x faster per step even on CPU |
| `spray_jet2d`, `spray_jet2d_gpu` | spray step 1: parcels injected from a slot into quiescent gas, both paths on identical input | single droplet matches the drag ODE to 7e-8 of U0 on both; momentum budget closed to 1e-3; the two paths agree to 0.8 % in momenta, 0.7 % in centreline velocity |
| `el_ibm_collision` | EL particles fired at a free IBM sphere | no particle enters the sphere, the sphere picks up the cloud's momentum (u 0 → 0.27), total momentum within 0.3 %, IBM slip ≈ 2 % of the body speed |

![Collision demo summary](docs/img/collision_summary.png)

![Spray jet, both particle paths](docs/img/spray_jet2d.png)

## Status

- [x] ppiclF linked into a nekRS UDF (host-staged; ppiclF runs on CPU)
- [x] One-way and two-way momentum coupling
- [x] Volume-fraction (QTL) divergence
- [x] IBM rigid spheres (NekIBM-style direct forcing) on ppiclF markers
- [x] EL–EL and EL–IBM soft-sphere contact
- [x] Device EL path (`lpm_t` + projection/filter kernels), two-way + QTL
- [x] Spray step 1: parcels (`N_p` droplets each) and run-time injection, both paths
- [ ] Spray steps 2–3: droplet deformation drag (TAB), breakup (TAB / KH-RT), dense-regime corrections
- [ ] IBM–IBM contact, particle–wall contact (ppiclF walls are not wired yet)
- [ ] Device path: contact (cell list), IBM markers, implicit Helmholtz filter

## License

BSD 3-Clause, see [LICENSE](LICENSE). The submodules keep their own licenses
(nekRS and NekIBM: UChicago Argonne BSD-style; ppiclF: MIT). The patches in
`patches/ppiclF/` modify ppiclF source and remain under ppiclF's MIT license,
[patches/ppiclF/LICENSE.ppiclF](patches/ppiclF/LICENSE.ppiclF).
