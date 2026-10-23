#!/usr/bin/env python3
"""Compare cases/jet_bickley runs with the Bickley plane laminar jet.

    python3 scripts/jet_bickley.py <run dir> [<run dir> ...] [--nu 5e-4] [--uco 0.05] [--J 0.02]

Each run dir holds field_final.csv (x,y,z,u,v,w on the GLL points). The
excess velocity u - u_co is z-averaged, the centreline value u_c(x) and the
half-width b(x) (where u - u_co = u_c/2) are extracted, and compared with
    u_c = 0.4543 (J^2 / (nu (x - x0)))^(1/3),   b = 1.2465 / 0.2752 ... see below
    u   = u_c sech^2(0.2752 (J / (nu^2 (x-x0)^2))^(1/3) y)
where x0 is the virtual origin fitted on u_c^-3 (linear in x). Also checks
that the excess momentum flux  J(x) = int u (u - u_co) dy  is constant.
"""
import argparse
import os

import numpy as np


def load(case):
    d = np.genfromtxt(os.path.join(case, "field_final.csv"), delimiter=",", names=True)
    xs, ys = np.unique(np.round(d["x"], 6)), np.unique(np.round(d["y"], 6))
    ix = np.searchsorted(xs, np.round(d["x"], 6))
    iy = np.searchsorted(ys, np.round(d["y"], 6))
    acc = np.zeros((len(ys), len(xs)))
    cnt = np.zeros_like(acc)
    np.add.at(acc, (iy, ix), d["u"])
    np.add.at(cnt, (iy, ix), 1)
    return xs, ys, acc / np.maximum(cnt, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--nu", type=float, default=5e-4)
    ap.add_argument("--uco", type=float, default=0.05)
    ap.add_argument("--J", type=float, default=0.02)
    ap.add_argument("--xmin", type=float, default=0.5, help="start of the fit/compare window")
    ap.add_argument("--xmax", type=float, default=1.7)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    for case in a.cases:
        xs, ys, U = load(case)
        E = U - a.uco  # excess velocity
        iyc = np.argmin(np.abs(ys))
        uc = E[iyc]
        # half-width: first y > 0 where E drops below uc/2 (linear interpolation)
        b = np.full_like(uc, np.nan)
        for j in range(len(xs)):
            prof = E[iyc:, j]
            k = np.where(prof < 0.5 * uc[j])[0]
            if len(k) and k[0] > 0 and uc[j] > 1e-6:
                k0 = k[0]
                y0, y1, p0, p1 = ys[iyc + k0 - 1], ys[iyc + k0], prof[k0 - 1], prof[k0]
                b[j] = y0 + (0.5 * uc[j] - p0) / (p1 - p0) * (y1 - y0)
        # excess momentum flux J(x) = int u (u - uco) dy (trapezoid on the GLL rows)
        Jx = np.trapezoid(U * E, ys, axis=0)

        w = (xs >= a.xmin) & (xs <= a.xmax)
        # Bickley: uc^-3 = (nu / (0.4543^3 J^2)) (x - x0): fit slope and x0
        A = np.vstack([xs[w], np.ones(w.sum())]).T
        slope, icpt = np.linalg.lstsq(A, uc[w] ** -3, rcond=None)[0]
        x0 = -icpt / slope
        J_fit = np.sqrt(a.nu / (0.4543**3 * slope))
        uc_th = 0.4543 * (a.J**2 / (a.nu * (xs - x0))) ** (1.0 / 3.0)
        # sech^2(0.2752 (J/(nu^2 x^2))^(1/3) y) = 1/2 -> y_half = 0.8814 / (0.2752 (J/(nu^2 x^2))^(1/3))
        b_th = 0.8814 / (0.2752 * (a.J / (a.nu**2 * (xs - x0) ** 2)) ** (1.0 / 3.0))

        lab = os.path.basename(case.rstrip("/"))
        err_uc = np.abs(uc[w] / uc_th[w] - 1).max()
        err_b = np.nanmax(np.abs(b[w] / b_th[w] - 1))
        print(f"{lab}: virtual origin x0 {x0:.3f}, J from u_c decay {J_fit:.4f} (set {a.J}), "
              f"u_c vs Bickley max rel err {err_uc:.3f}, half-width max rel err {err_b:.3f}, "
              f"J(x) = int u(u-u_co)dy over window: mean {Jx[w].mean():.4f}, spread {Jx[w].std() / Jx[w].mean():.3f}")

        axs[0, 0].plot(xs, uc, label=f"{lab}")
        axs[0, 0].plot(xs[w], uc_th[w], "k--", lw=0.8, label="Bickley" if case == a.cases[0] else None)
        axs[0, 1].plot(xs, b, label=lab)
        axs[0, 1].plot(xs[w], b_th[w], "k--", lw=0.8, label="Bickley" if case == a.cases[0] else None)
        axs[1, 0].plot(xs, Jx, label=lab)
        # profile at x = 1.2 in similarity variables
        j = np.argmin(np.abs(xs - 1.2))
        eta = 0.2752 * (a.J / (a.nu**2 * (xs[j] - x0) ** 2)) ** (1.0 / 3.0) * ys
        axs[1, 1].plot(eta, E[:, j] / uc[j], label=f"{lab}, x = {xs[j]:.2f}")
    axs[1, 0].axhline(a.J, color="k", ls="--", lw=0.8, label="injected J")
    e = np.linspace(-3, 3, 200)
    axs[1, 1].plot(e, 1 / np.cosh(e) ** 2, "k--", lw=0.8, label="sech$^2$")
    axs[0, 0].set_ylabel("centreline excess $u_c$")
    axs[0, 1].set_ylabel("half-width $b$")
    axs[1, 0].set_ylabel(r"$\int u(u-u_{co})\,dy$")
    axs[1, 1].set_xlabel(r"$\eta = 0.2752\,(J/\nu^2 x^2)^{1/3} y$")
    axs[1, 1].set_ylabel("$(u-u_{co})/u_c$")
    axs[1, 1].set_xlim(-3, 3)
    for ax in axs.flat[:3]:
        ax.set_xlabel("x")
    for ax in axs.flat:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("Plane laminar jet from EL parcels vs a fixed body force vs the Bickley similarity solution")
    fig.tight_layout()
    out = a.out or os.path.join(a.cases[0], "bickley.png")
    fig.savefig(out, dpi=120)
    print("wrote", out)


if __name__ == "__main__":
    main()
