#!/usr/bin/env python3
"""NekIBM's reference cases on this repo's IBM, against NekIBM's own references.

    python3 scripts/settling_reference.py --settle <log> [--settle-fine <log>] --re100 <log>
        --re10 <log> --re10-fine <log> [--out docs/img/nekibm_refs.png]

(a) cases/ibm_settling_sphere[_fine] (NekIBM Set1sphere = Mordant & Pinton
    2000, 0.5 mm glass bead in water, Ga 49.1, rho_p/rho_f 2.56): body
    velocity against the measured terminal Re_t = 41 (u_t = 0.074 m/s) and
    the standard-drag ODE (Schiller-Naumann + added mass, no history force).
(b) cases/ibm_sphere_re100 (NekIBM fixsphere/Re100 box): C_d history against
    Johnson & Patel (1999) 1.09 and Schiller-Naumann 1.09.
(c) cases/ibm_fixed_sphere[_fine] (NekIBM IBM_test/simple, Re 10): C_d with
    the Gaussian FWHM 0.08 d (two GLL spacings) and 0.04 d (one), against
    Schiller-Naumann 4.15 / standard drag 4.3.
"""
import argparse
import re

import numpy as np

D, RHO, G, NU = 1.0 / 6.0, 2.56, 9.81, 0.005416368
RE_T_MP = 41.0  # Mordant & Pinton (2000), 0.5 mm glass bead


def ode(t_end, dt=1e-4):
    """Standard-drag settling ODE: (rho + 1/2) dv/dt = (rho - 1) g - 3/4 Cd |v| v / d."""
    t, v = [0.0], [0.0]
    while t[-1] < t_end:
        re = abs(v[-1]) * D / NU
        drag = 0.75 * (24.0 / re * (1 + 0.15 * re**0.687)) * abs(v[-1]) * v[-1] / D if re > 0 else 0.0
        v.append(v[-1] + dt * ((RHO - 1) * G - drag) / (RHO + 0.5))
        t.append(t[-1] + dt)
    return np.array(t), np.array(v)


def parse(path, pat):
    rows = []
    for l in open(path, errors="ignore"):
        m = re.search(pat, l)
        if m:
            rows.append([float(x) for x in m.groups()])
    return np.array(rows)


SETTLE = r"^settle: step \d+\s+t\s+([\d.eE+-]+)\s+y\s+([\d.eE+-]+)\s+v\s+([\d.eE+-]+)\s+Re\s+([\d.eE+-]+)"
CD = r"^sphere: step \d+\s+t\s+([\d.eE+-]+).*Cd\s+([\d.eE+-]+)\s+slip\s+([\d.eE+-]+)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--settle", required=True)
    ap.add_argument("--settle-fine")
    ap.add_argument("--re100", required=True)
    ap.add_argument("--re10", required=True)
    ap.add_argument("--re10-fine")
    ap.add_argument("--out", default="nekibm_refs.png")
    a = ap.parse_args()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    vt_mp = RE_T_MP * NU / D
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.4))

    ax = axs[0]
    for path, lab, st in ((a.settle, "N = 7, FWHM 0.12 d", "C0-"), (a.settle_fine, "N = 8, FWHM 0.10 d", "C1-")):
        if not path:
            continue
        S = parse(path, SETTLE)
        t_o, v_o = ode(S[-1, 0])
        print(f"settling {lab}: v {-S[-1, 2]:.4f} (Re {S[-1, 3]:.2f}) at t {S[-1, 0]:.3f}; ODE {v_o[-1]:.4f} (Re {v_o[-1] * D / NU:.2f}); "
              f"Mordant & Pinton u_t {vt_mp:.4f} (Re 41); ours/MP {-S[-1, 2] / vt_mp:.3f}, ours/ODE {-S[-1, 2] / v_o[-1]:.3f}")
        ax.plot(S[:, 0], -S[:, 2], st, label=f"this repo's IBM, {lab} (NekIBM Set1sphere box)")
    t_o, v_o = ode(1.0)
    ax.plot(t_o, v_o, "k--", lw=0.9, label="standard-drag ODE (Schiller-Naumann + added mass)")
    ax.axhline(vt_mp, color="C3", ls=":", label="Mordant & Pinton (2000) measured $u_t$, Re$_t$ = 41")
    ax.set_xlabel("t  (x 0.054 s)")
    ax.set_ylabel("settling velocity  (x 0.0556 m/s)")
    ax.set_title("(a) glass sphere settling, Ga = 49, $\\rho_p/\\rho_f$ = 2.56")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    ax = axs[1]
    C = parse(a.re100, CD)
    print(f"Re 100: Cd {C[-1, 1]:.4f} at U t / d = {C[-1, 0] / 0.01:.1f}, slip {C[-1, 2]:.3f}; Johnson & Patel 1.09, ratio {C[-1, 1] / 1.09:.3f}")
    ax.plot(C[:, 0] / 0.01, C[:, 1], "C0-", label="this repo's IBM, N = 5, FWHM 0.12 d (NekIBM fixsphere/Re100 box)")
    ax.axhline(1.09, color="k", ls="--", lw=0.9, label="Johnson & Patel (1999): 1.09")
    ax.set_ylim(0.8, 2.0)
    ax.set_xlabel("U t / d")
    ax.set_ylabel("$C_d$")
    ax.set_title("(b) fixed sphere, Re = 100 (8d x 15d x 8d box)")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    ax = axs[2]
    for path, lab, st in ((a.re10, "FWHM 0.08 d = 2 GLL spacings", "C0-"), (a.re10_fine, "FWHM 0.04 d = 1 GLL spacing", "C1-")):
        if not path:
            continue
        C = parse(path, CD)
        print(f"Re 10 {lab}: Cd {C[-1, 1]:.4f} at U t / d = {C[-1, 0] / 0.01:.1f}, marker slip {C[-1, 2]:.3f}")
        ax.plot(C[:, 0] / 0.01, C[:, 1], st, label=f"this repo's IBM, N = 5, {lab}")
    ax.axhline(4.15, color="k", ls="--", lw=0.9, label="Schiller-Naumann, unbounded: 4.15")
    ax.axhline(4.30, color="k", ls=":", lw=0.9, label="standard drag curve, unbounded: 4.30")
    ax.set_ylim(3, 7)
    ax.set_xlabel("U t / d")
    ax.set_ylabel("$C_d$")
    ax.set_title("(c) fixed sphere, Re = 10 (4d box): Gaussian width")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
