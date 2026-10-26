# `settling_twoway` on Nek5000 v19 + the original (unpatched) ppiclF

Reference run for the validation in `docs/design.md`, section 5: the same
10 000 particles (rho_p 100, d_p 0.01, from rest, g = -1 in z) in the same
periodic unit box (6^3 elements, N = 5, nu = 0.01, dt = 1e-3, 200 steps),
two-way coupled through ppiclF's own Gaussian projection (FWHM 0.15), with
ppiclF's Nek5000 example wiring (`userf` = F / (1 - phi_p) / rho,
`userchk` = interpolate + `IntegrateParticle`).

```bash
# ppiclF built from external/ppiclF WITHOUT patches/ppiclF (plain `make`),
# Nek5000 v19 in $NEK_SOURCE_ROOT
./makenek settle && mpirun -np 4 ./nek5000   # SESSION.NAME: "settle" and this directory
grep twoway: logfile
```

Prints `twoway: step t Pf Pp sumf total` every 20 steps, the same budget
as `cases/settling_twoway/settling.udf`. Positions come from a different
random stream than the nekRS case, so the comparison is of the cloud
statistics, not particle by particle. Nek5000 limits here: the box path in
`SESSION.NAME` must be short (< ~130 characters) and `include "PPICLF"` in
the `.usr` duplicates DATA-initialised commons at link time on macOS, hence
`ppiclf_user_budget` lives in `user_routines/ppiclf_user.f`.
