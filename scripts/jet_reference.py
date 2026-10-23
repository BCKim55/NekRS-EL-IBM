#!/usr/bin/env python3
"""Checks for cases/spray_jet2d (ppiclF path) and spray_jet2d_gpu (device path).

    python3 scripts/jet_reference.py single <log> [<log> ...]
        single-droplet runs: droplet velocity vs the drag ODE
        m dv/dt = 3 pi mu d (1 + 0.15 Re^0.687) (0 - v), integrated with RK4
        at 1e-3 of the fluid step
    python3 scripts/jet_reference.py jet <case dir> [<case dir> ...]
        jet runs: momentum budget from the log and figures from field_final.csv
"""
import csv
import glob
import os
import re
import sys

import numpy as np

MU, RHOF, DP, RHOP, U0 = 5e-4, 1.0, 0.002, 1000.0, 1.0


def drag_ode(t_end, dt):
    m = RHOP * np.pi / 6 * DP**3

    def f(v):
        re = abs(v) * DP * RHOF / MU
        return -3 * np.pi * MU * DP * (1 + 0.15 * re**0.687) * v / m

    v, t, out = U0, 0.0, []
    h = dt / 1000
    n = int(round(t_end / h))
    for i in range(1, n + 1):
        k1 = f(v)
        k2 = f(v + 0.5 * h * k1)
        k3 = f(v + 0.5 * h * k2)
        k4 = f(v + h * k3)
        v += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if i % 1000 == 0:
            out.append(v)
    return np.array(out)


def single(logs):
    for log in logs:
        rows = [l.split() for l in open(log, errors="ignore") if l.startswith("drop:") and int(l.split()[1]) > 0]
        step = np.array([int(r[1]) for r in rows])
        t = np.array([float(r[2]) for r in rows])
        v = np.array([float(r[3]) for r in rows])
        dt = t[1] - t[0]
        ref = drag_ode(t[-1], dt)[step - 1]
        err = np.abs(v - ref) / U0
        tau = RHOP * DP**2 / (18 * MU)
        print(f"{os.path.basename(log)}: {len(v)} steps, dt {dt:g}, tau_p(Stokes) {tau:.4f}, "
              f"v(t_end) {v[-1]:.6f} ref {ref[-1]:.6f}, max |dv|/U0 {err.max():.2e}, "
              f"Re_p0 {U0 * DP * RHOF / MU:.1f}")


def budget(logdir):
    log = sorted(glob.glob(os.path.join(logdir, "*.log")))
    return log


def read_log(log):
    """jet runs: (t, Pgas, Pdrop, injected, budget, fluid) from the report lines."""
    rows = []
    for l in open(log, errors="ignore"):
        if l.startswith("jet: step"):
            f = l.split()
            rows.append([float(f[4]), float(f[8]), float(f[10]), float(f[12]), float(f[14]), float(f[16])])
    return np.array(rows)


def jet(cases):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(len(cases) + 2, 1, figsize=(9, 3.0 * len(cases) + 6.4))
    prof = []
    for case in cases:
        logs = sorted(glob.glob(os.path.join(case, "*.log")))
        if logs:
            r = read_log(logs[-1])
            ax = axs[-2]
            lab = os.path.basename(case)
            ax.plot(r[:, 0], r[:, 1], "-", label=f"{lab}: gas")
            ax.plot(r[:, 0], r[:, 2], "--", label=f"{lab}: parcels")
            ax.plot(r[:, 0], r[:, 1] + r[:, 2], ":", label=f"{lab}: sum")
            ax.plot(r[:, 0], r[:, 3], "k-", lw=0.8, label="injected" if case == cases[0] else None)
            print(f"{lab}: t {r[-1, 0]:.2f}  P_gas {r[-1, 1]:.4e}  P_parcels {r[-1, 2]:.4e}  injected {r[-1, 3]:.4e}  "
                  f"budget {r[-1, 4]:+.2e}  max|budget| {np.abs(r[:, 4]).max():.2e}  fluid {r[-1, 5]:+.2e}")
    axs[-2].set_xlabel("t")
    axs[-2].set_ylabel("x momentum")
    axs[-2].legend(fontsize=7, ncol=2)
    axs[-2].grid(alpha=0.3)
    for k, case in enumerate(cases):
        d = np.genfromtxt(os.path.join(case, "field_final.csv"), delimiter=",", names=True)
        x, y, u = d["x"], d["y"], d["u"]
        xs, ys = np.unique(np.round(x, 6)), np.unique(np.round(y, 6))
        ix = np.searchsorted(xs, np.round(x, 6))
        iy = np.searchsorted(ys, np.round(y, 6))
        acc = np.zeros((len(ys), len(xs)))
        cnt = np.zeros_like(acc)
        np.add.at(acc, (iy, ix), u)
        np.add.at(cnt, (iy, ix), 1)
        U = acc / np.maximum(cnt, 1)
        ax = axs[k]
        im = ax.pcolormesh(xs, ys, U, cmap="viridis", shading="nearest", vmin=0, vmax=0.5)
        pf = os.path.join(case, "parcels_final.csv")
        if os.path.exists(pf):
            pp = np.genfromtxt(pf, delimiter=",", names=True)
            ax.scatter(pp["x"], pp["y"], s=1, c="w", alpha=0.6)
        ax.set_aspect("equal")
        ax.set_ylim(0.2, 0.8)
        ax.set_title(f"{os.path.basename(case)}: gas $u_x$ (z-averaged), parcels in white")
        fig.colorbar(im, ax=ax)
        iyc = np.argmin(np.abs(ys - 0.5))
        axs[-1].plot(xs, U[iyc], label=os.path.basename(case))
        prof.append((xs, U[iyc]))
    axs[-1].set_xlabel("x")
    axs[-1].set_ylabel("centreline gas $u_x$")
    axs[-1].legend()
    axs[-1].grid(alpha=0.3)
    out = os.path.join(cases[0], "jet_fields.png")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    print("wrote", out)
    if len(prof) == 2:
        a, b = prof
        n = min(len(a[1]), len(b[1]))
        print(f"centreline gas u_x: max {a[1].max():.4f} vs {b[1].max():.4f}, "
              f"max |diff| {np.abs(a[1][:n] - b[1][:n]).max():.4f}, rms diff {np.sqrt(np.mean((a[1][:n] - b[1][:n])**2)):.4f}")


if __name__ == "__main__":
    if sys.argv[1] == "single":
        single(sys.argv[2:])
    else:
        jet(sys.argv[2:])
