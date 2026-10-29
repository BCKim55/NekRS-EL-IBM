#!/usr/bin/env python3
"""Time-averaged u contours of the spray_jet_water sweep: droplet-driven jet
(top, with the final parcel positions) over the gas velocity-inlet jet of the
same momentum flux (bottom), one column per loading, shared colour scale per
column.

    python3 scripts/spray_jet_contours.py build/spray/p1:build/spray/i1 ... [--out fig.png]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from spray_jet_water import load_named, par  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pairs", nargs="+")
    ap.add_argument("--out", default="spray_jet_contours.png")
    a = ap.parse_args()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(a.pairs)
    fig, axs = plt.subplots(2, n, figsize=(6.5 * n, 6.2), squeeze=False)
    for col, pair in enumerate(a.pairs):
        pdir, idir = pair.split(":")
        P = par(pdir)
        J = P["rhop"] * P["phi_inj"] * P["u_jet"] ** 2 * P["w_slot"]
        fields = [load_named(os.path.join(d, "field_avg.csv")) for d in (pdir, idir)]
        vmax = max(U.max() for _, _, U in fields)
        vmin = min(min(U.min() for _, _, U in fields), -0.05 * vmax)
        levels = np.linspace(vmin, vmax, 41)
        for row, (d, (xs, ys, U)) in enumerate(zip((pdir, idir), fields)):
            ax = axs[row, col]
            c = ax.contourf(xs, ys, U, levels=levels, cmap="RdBu_r", extend="both")
            ax.contour(xs, ys, U, levels=[0], colors="k", linewidths=0.6)
            if row == 0:
                pf = os.path.join(d, "parcels_final.csv")
                if os.path.exists(pf):
                    p = np.genfromtxt(pf, delimiter=",", names=True)
                    ax.plot(p["x"][::4], p["y"][::4], ".", color="k", ms=0.4, alpha=0.15)
            ax.set_aspect("equal")
            ax.set_xlim(xs[0], xs[-1])
            ax.set_ylim(ys[0], ys[-1])
            kind = "droplet-driven (dots: parcels at the end)" if row == 0 else "gas velocity inlet, same J"
            ax.set_title(f"{os.path.basename(d)}: {kind}\nJ = {J:.3g}, u_co = {P['u_co']:g}, time-averaged u (black line: u = 0)", fontsize=9)
            if row == 1:
                ax.set_xlabel("x")
            ax.set_ylabel("y")
            fig.colorbar(c, ax=ax, fraction=0.03, pad=0.02)
    fig.tight_layout()
    fig.savefig(a.out, dpi=110)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
