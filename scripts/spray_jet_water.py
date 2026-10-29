#!/usr/bin/env python3
"""Compare cases/spray_jet_water runs: droplet-driven plane jets (p*) against
the gas-only velocity-inlet jet of the same momentum flux (i*), the Bickley
similarity solution and a liquid jet (which keeps its velocity and all of
its momentum).

    python3 scripts/spray_jet_water.py build/spray/p1:build/spray/i1 [build/spray/p2:build/spray/i2 ...] [--out fig.png]

Per run pair: gas centreline excess velocity u_c(x), mean droplet velocity
u_d(x) (parcels_final.csv binned in x), gas momentum flux J_g(x) = int u(u -
u_co) dy, liquid momentum flux J_l(x) = (rho_p/rho_f) phi u_jet w u_d(x)
(the liquid mass flux is conserved), and the half-width b(x). Fields are
the time average (field_avg.csv) when present.
"""
import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from jet_bickley import load  # noqa: E402


def par(path):
    p = {}
    for l in open(os.path.join(path, "jet.par")):
        m = re.match(r"^\s*([A-Za-z_]+)\s*=\s*([-+0-9.eE]+)\s*$", l)
        if m:
            p[m.group(1)] = float(m.group(2))
    return p


def field(path):
    f = os.path.join(path, "field_avg.csv")
    return load_named(f) if os.path.exists(f) else load(path)


def load_named(f):
    d = np.genfromtxt(f, delimiter=",", names=True)
    xs, ys = np.unique(np.round(d["x"], 6)), np.unique(np.round(d["y"], 6))
    ix = np.searchsorted(xs, np.round(d["x"], 6))
    iy = np.searchsorted(ys, np.round(d["y"], 6))
    acc = np.zeros((len(ys), len(xs)))
    cnt = np.zeros_like(acc)
    np.add.at(acc, (iy, ix), d["u"])
    np.add.at(cnt, (iy, ix), 1)
    return xs, ys, acc / np.maximum(cnt, 1)


def analyse(path, uco):
    xs, ys, U = field(path)
    E = U - uco
    iyc = np.argmin(np.abs(ys))
    uc = E[iyc]
    b = np.full_like(uc, np.nan)
    for j in range(len(xs)):
        prof = E[iyc:, j]
        k = np.where(prof < 0.5 * uc[j])[0]
        if len(k) and k[0] > 0 and uc[j] > 0.05 * uc.max():
            k0 = k[0]
            y0, y1, p0, p1 = ys[iyc + k0 - 1], ys[iyc + k0], prof[k0 - 1], prof[k0]
            b[j] = y0 + (0.5 * uc[j] - p0) / (p1 - p0) * (y1 - y0)
    Jx = np.trapezoid(U * E, ys, axis=0)
    return xs, uc, b, Jx, U.min(axis=0)


def droplets(path, xs):
    f = os.path.join(path, "parcels_final.csv")
    d = np.genfromtxt(f, delimiter=",", names=True)
    edges = np.linspace(xs[0], xs[-1], 61)
    xc = 0.5 * (edges[1:] + edges[:-1])
    ud = np.full(len(xc), np.nan)
    for i in range(len(xc)):
        m = (d["x"] >= edges[i]) & (d["x"] < edges[i + 1])
        if m.sum() >= 5:
            ud[i] = d["u"][m].mean()
    return xc, ud, len(d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pairs", nargs="+", help="parcel_run:inlet_run")
    ap.add_argument("--nu", type=float, default=5e-4)
    ap.add_argument("--xmin", type=float, default=1.0)
    ap.add_argument("--xmax", type=float, default=2.7)
    ap.add_argument("--out", default="spray_jet_water.png")
    a = ap.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(a.pairs)
    fig, axs = plt.subplots(4, n, figsize=(5.2 * n, 13), squeeze=False)
    for col, pair in enumerate(a.pairs):
        pdir, idir = pair.split(":")
        P = par(pdir)
        uco, uj, w, phi, rr = P["u_co"], P["u_jet"], P["w_slot"], P["phi_inj"], P["rhop"] / 1.0
        J = rr * phi * uj * uj * w
        xs, ucp, bp, Jp, ump = analyse(pdir, uco)
        _, uci, bi, Ji, umi = analyse(idir, uco)
        xc, ud, nd = droplets(pdir, xs)
        Jl = rr * phi * uj * w * ud
        wdw = (xs >= a.xmin) & (xs <= a.xmax)
        # Bickley with x0 fitted to the gas-only inlet jet (u_c^-3 linear in x) and its measured J;
        # the droplet-driven jet is still being accelerated in the window, so it has no similarity origin
        A = np.vstack([xs[wdw], np.ones(wdw.sum())]).T
        slope, icpt = np.linalg.lstsq(A, uci[wdw] ** -3, rcond=None)[0]
        x0 = -icpt / slope
        Jm = Ji[wdw].mean()
        uc_th = 0.4543 * (Jm**2 / (a.nu * np.maximum(xs - x0, 1e-9))) ** (1 / 3)
        b_th = 0.8814 / (0.2752 * (Jm / (a.nu**2 * np.maximum(xs - x0, 1e-9) ** 2)) ** (1 / 3))
        err = np.abs(uci[wdw] / uc_th[wdw] - 1).max()
        print(f"{os.path.basename(pdir)} phi {phi:g} J_inj {J:.4f} u_co {uco}: inlet jet J_gas window mean {Jm:.4f} "
              f"({Jm / J:.1%} of injected), x0 {x0:.3f}, u_c vs Bickley max err {err:.3f}, min u {umi[wdw].min():.3f}")
        print(f"   droplet-driven: J_gas window mean {Jp[wdw].mean():.4f} ({Jp[wdw].mean() / J:.1%} of injected), "
              f"u_c/u_c(inlet) at x = 2 {np.interp(2.0, xs, ucp) / np.interp(2.0, xs, uci):.3f}, "
              f"min u {ump[wdw].min():.3f} (backflow if < 0)")
        m = ~np.isnan(ud)
        print(f"   droplets: {nd} in the domain, u_d at nozzle+0.1 {np.interp(P['x_nozzle'] + 0.1, xc[m], ud[m]):.3f}, "
              f"at x = 2 {np.interp(2.0, xc[m], ud[m]):.3f}, at outflow {ud[m][-1]:.3f} (u_jet {uj}); "
              f"liquid momentum at outflow {Jl[m][-1] / J:.1%} of injected, gas {Jp[-1] / J:.1%}")

        ax = axs[0, col]
        ax.plot(xs, ucp, "C0", label="gas, droplet-driven")
        ax.plot(xs, uci, "C1", label="gas, velocity-inlet jet")
        ax.plot(xs[wdw], uc_th[wdw], "k--", lw=0.8, label="Bickley (x0, J from the inlet jet)")
        ax.set_ylabel("centreline excess $u_c$")
        ax.set_title(f"$\\phi_{{inj}}$ = {phi:g}, J = {J:.3g}, $u_{{co}}$ = {uco:g}, liquid/gas mass flux {rr * phi:.2g}", fontsize=10)
        ax = axs[1, col]
        ax.plot(xc, ud, "C2.-", label="droplets $u_d$ (mean in x bins)")
        ax.plot(xs, ucp + uco, "C0", label="gas centreline u")
        ax.axhline(uj, color="k", ls="--", lw=0.8, label="liquid jet: u = $u_{jet}$")
        ax.set_ylabel("velocity")
        ax.set_yscale("log")
        ax = axs[2, col]
        ax.plot(xs, Jp / J, "C0", label="gas $J_g$ / J, droplet-driven")
        ax.plot(xc, Jl / J, "C2.-", label="liquid $J_l$ / J")
        ax.plot(xc, (np.interp(xc, xs, Jp) + Jl) / J, "C3", label="gas + liquid")
        ax.plot(xs, Ji / J, "C1", label="gas $J_g$ / J, velocity-inlet jet")
        ax.axhline(1, color="k", ls="--", lw=0.8, label="liquid jet: all in the liquid")
        ax.set_ylabel("momentum flux / injected J")
        ax.set_ylim(-0.1, 1.3)
        ax = axs[3, col]
        ax.plot(xs, bp, "C0", label="droplet-driven")
        ax.plot(xs, bi, "C1", label="velocity-inlet jet")
        ax.plot(xs[wdw], b_th[wdw], "k--", lw=0.8, label="Bickley")
        ax.axhline(0.5 * w, color="k", ls=":", lw=0.8, label="liquid jet: w/2")
        ax.set_ylabel("half-width b")
        ax.set_xlabel("x")
        for ax in axs[:, col]:
            ax.axvline(P["x_nozzle"], color="gray", lw=0.5)
            ax.grid(alpha=0.3)
            ax.legend(fontsize=7)
    fig.suptitle("Droplet-driven plane jet (water-like droplets) vs a gas jet of the same momentum flux vs a liquid jet")
    fig.tight_layout()
    fig.savefig(a.out, dpi=120)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
