#!/usr/bin/env python3
"""Checks for cases/tab_breakup.

    python3 scripts/tab_reference.py <log> [--out fig.png]

1. Oscillator: the distortion history of the We < 12 droplet (lines "tab:")
   against the analytic damped oscillator about its equilibrium
   y_eq = (CF/CB) rho_g u^2 r / (CK sigma) = We/24 (so the undamped overshoot
   2 y_eq reaches 1 at We = 12, TAB's critical Weber number): frequency
   sqrt(omega^2 - zeta^2), omega^2 = 8 sigma/(rho_l r^3), zeta = 5 mu_l/(2 rho_l r^2).
   (We decays slowly through the drag; the check uses the first periods.)
2. Breakup time of every droplet (lines "tbreak:") against the Pilch & Erdman
   (1987) correlation for the total breakup time
   T = t_b U sqrt(rho_g/rho_l) / d:
     6 (We-12)^-0.25            12 < We <= 18
     2.45 (We-12)^0.25          18 < We <= 45
     14.1 (We-12)^-0.25         45 < We <= 351
     0.766 (We-12)^0.25         351 < We <= 2670
     5.5                        We > 2670
   TAB's own breakup time is the time for y to reach 1, i.e. the *initiation*
   of breakup, for which Pilch & Erdman give T_ini = 1.9 (We-12)^-0.25; the
   total breakup time above is 4-5 times longer. What must hold: no breakup
   below We = 12, T(TAB) close to T_ini, and the trend with We.
"""
import argparse

import numpy as np

DP, RHOP, SIGMA, MUL, RHOG = 0.002, 1000.0, 1e-4, 1e-3, 1.0


def pilch(we):
    we = np.asarray(we, dtype=float)
    T = np.full_like(we, np.nan)
    m = (we > 12) & (we <= 18)
    T[m] = 6 * (we[m] - 12) ** -0.25
    m = (we > 18) & (we <= 45)
    T[m] = 2.45 * (we[m] - 12) ** 0.25
    m = (we > 45) & (we <= 351)
    T[m] = 14.1 * (we[m] - 12) ** -0.25
    m = (we > 351) & (we <= 2670)
    T[m] = 0.766 * (we[m] - 12) ** 0.25
    T[we > 2670] = 5.5
    return T


def pilch_init(we):
    we = np.asarray(we, dtype=float)
    return np.where(we > 12, 1.9 * np.abs(we - 12) ** -0.25, np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log")
    ap.add_argument("--out", default="tab_breakup.png")
    a = ap.parse_args()
    tab, tb = [], []
    for l in open(a.log, errors="ignore"):
        if l.startswith("tab: ") and not l.startswith("tab: omega"):
            tab.append([float(v) for v in l.split()[1:]])
        if l.startswith("tbreak:"):
            tb.append([float(v) for v in l.split()[1:]])
    tab, tb = np.array(tab), np.array(tb)
    tb = tb[np.argsort(tb[:, 0])]

    r = 0.5 * DP
    omega = np.sqrt(8 * SIGMA / (RHOP * r**3))
    zeta = 5 * MUL / (2 * RHOP * r**2)
    omega_d = np.sqrt(omega**2 - zeta**2)
    # measured period from successive maxima of y (first droplet, We = 10)
    t, y = tab[:, 1], tab[:, 2]
    peaks = [i for i in range(1, len(y) - 1) if y[i] > y[i - 1] and y[i] >= y[i + 1]]
    periods = np.diff(t[peaks][:4])
    print(f"oscillator: analytic damped period {2 * np.pi / omega_d:.4f} (omega {omega:.3f}, zeta {zeta:.3f}), "
          f"measured first periods {np.array2string(periods, precision=4)}, y max {y.max():.3f}, "
          f"y_eq (We/24 at release) {tb[0, 1] / 24:.3f}, first peak / y_eq {y.max() / (tb[0, 1] / 24):.3f} "
          f"(undamped 2, damped {1 + np.exp(-zeta * np.pi / omega_d):.3f})")

    print("breakup times, T = t_b U sqrt(rho_g/rho_l)/d:")
    print("  We    t_b(TAB)   T(TAB)   T_ini(P&E)   T_total(P&E)   d_after/d   N_p")
    for gid, we, tbr, vx, d, np_ in tb:
        U = np.sqrt(we * SIGMA / (RHOG * DP))
        T = tbr * U * np.sqrt(RHOG / RHOP) / DP if tbr >= 0 else np.nan
        print(f"  {we:5.0f}  {tbr if tbr >= 0 else float('nan'):8.4f}  {T:7.3f}  {pilch_init([we])[0]:9.3f}  {pilch([we])[0]:12.3f}      {d / DP:6.3f}  {np_:7.1f}")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    axs[0].plot(t, y, label="TAB y(t), We = %.0f" % tb[0, 1])
    yeq = tb[0, 1] / 24
    axs[0].plot(t, yeq * (1 - np.exp(-zeta * t) * (np.cos(omega_d * t) + zeta / omega_d * np.sin(omega_d * t))), "k--", lw=0.8,
                label="damped oscillator about We/24 (fixed We)")
    axs[0].set_xlabel("t")
    axs[0].set_ylabel("distortion y")
    axs[0].legend(fontsize=8)
    axs[0].grid(alpha=0.3)
    we = np.geomspace(12.01, 400, 200)
    axs[1].plot(we, pilch(we), "k-", label="Pilch & Erdman: total breakup time")
    axs[1].plot(we, pilch_init(we), "k--", label="Pilch & Erdman: initiation time")
    m = tb[:, 2] >= 0
    U = np.sqrt(tb[:, 1] * SIGMA / (RHOG * DP))
    T = tb[:, 2] * U * np.sqrt(RHOG / RHOP) / DP
    axs[1].plot(tb[m, 1], T[m], "o", label="TAB first breakup (this code)")
    axs[1].axvline(12, color="C3", ls=":", label="We = 12")
    axs[1].set_xscale("log")
    axs[1].set_xlabel("We")
    axs[1].set_ylabel(r"$T = t_b U \sqrt{\rho_g/\rho_l}/d$")
    axs[1].legend(fontsize=8)
    axs[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(a.out, dpi=120)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
