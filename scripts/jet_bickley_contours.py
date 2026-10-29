#!/usr/bin/env python3
"""Contours of the excess velocity u - u_co of two cases/jet_bickley runs and
their difference, plus the Bickley solution evaluated with each run's fitted
virtual origin.

    python3 scripts/jet_bickley_contours.py build/bick/parcels build/bick/source [--out contours.png]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from jet_bickley import load  # noqa: E402


def fit_x0(xs, uc, w):
    A = np.vstack([xs[w], np.ones(w.sum())]).T
    slope, icpt = np.linalg.lstsq(A, uc[w] ** -3, rcond=None)[0]
    return -icpt / slope


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs=2)
    ap.add_argument("--nu", type=float, default=5e-4)
    ap.add_argument("--uco", type=float, default=0.05)
    ap.add_argument("--J", type=float, default=0.02)
    ap.add_argument("--xmin", type=float, default=0.9)
    ap.add_argument("--xmax", type=float, default=1.8)
    ap.add_argument("--out", default="bickley_contours.png")
    a = ap.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    runs = []
    for c in a.cases:
        xs, ys, U = load(c)
        E = U - a.uco
        uc = E[np.argmin(np.abs(ys))]
        w = (xs >= a.xmin) & (xs <= a.xmax)
        x0 = fit_x0(xs, uc, w)
        runs.append((os.path.basename(c.rstrip("/")), xs, ys, E, x0))

    (n1, xs, ys, E1, x01), (n2, _, _, E2, x02) = runs
    X, Y = np.meshgrid(xs, ys)
    vmax = max(E1.max(), E2.max())
    lev = np.linspace(-0.05 * vmax, vmax, 22)

    fig, axs = plt.subplots(2, 2, figsize=(15, 7.5), sharex=True, sharey=True)
    axs = list(axs.flat)
    for ax, (name, E, x0) in zip(list(axs)[:2], [(n1, E1, x01), (n2, E2, x02)]):
        cf = ax.contourf(X, Y, E, lev, cmap="RdBu_r", extend="both")
        ax.contour(X, Y, E, [0.0], colors="k", linewidths=0.6)
        fig.colorbar(cf, ax=ax, label="$u - u_{co}$")
        ax.set_title(f"{name}: excess velocity (black: $u - u_{{co}} = 0$), virtual origin x0 = {x0:.3f}")
    # Bickley with the parcels run's x0, inside the comparison window
    xb = xs[(xs >= a.xmin) & (xs <= a.xmax)]
    Xb, Yb = np.meshgrid(xb, ys)
    ucb = 0.4543 * (a.J**2 / (a.nu * (Xb - x01))) ** (1 / 3)
    eta = 0.2752 * (a.J / (a.nu**2 * (Xb - x01) ** 2)) ** (1 / 3) * Yb
    Eb = ucb / np.cosh(eta) ** 2
    ax = list(axs)[2]
    cf = ax.contourf(Xb, Yb, Eb, lev, cmap="RdBu_r", extend="both")
    fig.colorbar(cf, ax=ax, label="$u - u_{co}$")
    ax.set_title(f"Bickley with J = {a.J}, x0 = {x01:.3f} ({n1}), only inside the window {a.xmin} < x < {a.xmax}")
    ax = list(axs)[3]
    D = E1 - E2
    dm = np.abs(D[:, (xs >= a.xmin)]).max()
    cf = ax.contourf(X, Y, D, np.linspace(-dm, dm, 21), cmap="PuOr", extend="both")
    fig.colorbar(cf, ax=ax, label=f"({n1}) - ({n2})")
    ax.set_title(f"difference {n1} - {n2} (colour range from x > {a.xmin}; the near field exceeds it)")
    for ax in axs:
        ax.axvline(a.xmin, color="k", ls=":", lw=0.6)
        ax.axvline(a.xmax, color="k", ls=":", lw=0.6)
        ax.set_ylabel("y")
        ax.set_aspect("equal")
    for ax in axs[2:]:
        ax.set_xlabel("x")
    fig.tight_layout()
    fig.savefig(a.out, dpi=120)
    print("wrote", a.out)
    j = np.argmin(np.abs(xs - 1.2))
    print(f"at x = 1.2: max |diff| {np.abs(D[:, j]).max():.4f} vs u_c {E1[np.argmin(np.abs(ys)), j]:.4f}; "
          f"min excess velocity {n1} {E1[:, j].min():.4f}, {n2} {E2[:, j].min():.4f}")


if __name__ == "__main__":
    main()
