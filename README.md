# NekRS-EL-IBM

Coupled Euler–Lagrange point particles and immersed-boundary (IBM) particles on
top of nekRS, using ppiclF for the Lagrangian side and NekIBM as the reference
for the IBM algorithms.

## Layout

| Path | What |
|---|---|
| `external/nekRS` | nekRS, pinned at v26.0 (`96b3cf9e`) |
| `external/ppiclF` | ppiclF, pinned at `72483db` (last upstream commit, 2020) |
| `external/NekIBM` | NekIBM, pinned at `9631f97`; reference only, not built |
| `patches/ppiclF/` | fixes applied to ppiclF at build time (submodule stays pristine) |
| `src/ppiclf_cbind.f` | C bindings for ppiclF, appended to each case's user file |
| `src/ppiclfCoupler.hpp` | nekRS-side coupling helpers, included from a `.udf` |
| `scripts/build_ppiclf.sh` | builds a case-specific `libppiclF.a` and `udf.env` |
| `scripts/run_case.sh` | builds ppiclF for a case and runs it with nekRS |
| `cases/` | verification and application cases |
| `docs/design.md` | coupling design, ppiclF findings, GPU roadmap |

## Build and run

nekRS must be installed from `external/nekRS` (or the same commit) with double
`dfloat`. The UDF must be compiled with the same C++ compiler and standard
library as nekRS; with Open MPI wrappers that means setting `OMPI_CXX` to the
compiler nekRS was built with (a libc++ UDF loaded into a libstdc++ nekRS
segfaults in `dlopen`).

```bash
git clone --recursive <this repo>
export NEKRS_HOME=$HOME/.local/nekrs
export OMPI_CXX=g++-16            # whatever built nekRS
scripts/run_case.sh cases/rotation_tracer rotation 4 --backend CPU
```

`run_case.sh` rebuilds `libppiclF.a` only when the case's `user_routines`,
`src/ppiclf_cbind.f` or the patches change, and always relinks the UDF.

## Status

- [x] ppiclF linked into a nekRS UDF (host-staged; ppiclF runs on CPU)
- [x] One-way tracer verification (`cases/rotation_tracer`), 1–8 ranks
- [x] Two-way momentum coupling (`cases/settling_twoway`), 1 and 4 ranks
- [ ] Volume-fraction (QTL) divergence via `nrs->userDivergence`
- [ ] IBM particles ported from NekIBM onto ppiclF
- [ ] EL–IBM collisions in one neighbour search
- [ ] GPU particle path (see `docs/design.md`)
