#!/usr/bin/env python3
"""Midplane (z = 0) contours of the IBM reference cases from plane_z0.csv.

    python3 scripts/plot_slices.py --re10 cases/ibm_fixed_sphere --re10-fine cases/ibm_fixed_sphere_fine
        --re100 cases/ibm_sphere_re100 [--re100-hires build/tuo] --settle cases/ibm_settling_sphere --settle-log build/c_settle.log
        [--out docs/img/ibm_slices.png]

Each panel: velocity magnitude on the GLL points of the z = 0 plane (tricontourf),
a few streamlines of the in-plane velocity, and the nominal sphere outline.
"""
import argparse
import os
import re

import numpy as np


def load(case):
    d = np.genfromtxt(os.path.join(case, "plane_z0.csv"), delimiter=",", names=True)
    # GLL points are duplicated on element faces; keep one copy per location
    key = np.round(d["x"], 9) + 1e3 * np.round(d["y"], 9)
    _, idx = np.unique(key, return_index=True)
    return {k: d[k][idx] for k in d.dtype.names}


def panel(ax, d, xc, yc, r, title, zoom=None, uref=1.0):
    import matplotlib.tri as mtri
    from matplotlib.patches import Circle

    x, y, u, v = d["x"], d["y"], d["u"], d["v"]
    tri = mtri.Triangulation(x, y)
    # the solid: triangles whose centroid lies inside the sphere are masked
    # out of the contours, the streamline field is undefined there, and the
    # disc is drawn filled on top (the IBM has fluid inside the shell; it is
    # rigid, not flow)
    cx = x[tri.triangles].mean(axis=1) - xc
    cy = y[tri.triangles].mean(axis=1) - yc
    tri.set_mask(cx * cx + cy * cy < r * r)
    mag = np.hypot(u, v) / uref
    cf = ax.tricontourf(tri, mag, levels=np.linspace(0, 1.25, 26), cmap="viridis", extend="max")
    # streamlines on a regular grid interpolated from the triangulation
    if zoom:
        xlo, xhi, ylo, yhi = zoom
    else:
        xlo, xhi, ylo, yhi = x.min(), x.max(), y.min(), y.max()
    gx, gy = np.meshgrid(np.linspace(xlo, xhi, 120), np.linspace(ylo, yhi, 160))
    iu = np.ma.filled(mtri.LinearTriInterpolator(tri, u)(gx, gy), np.nan)
    iv = np.ma.filled(mtri.LinearTriInterpolator(tri, v)(gx, gy), np.nan)
    inside = (gx - xc) ** 2 + (gy - yc) ** 2 < r * r
    iu[inside] = np.nan
    iv[inside] = np.nan
    ax.streamplot(gx, gy, np.nan_to_num(iu), np.nan_to_num(iv), color="w", linewidth=0.5, density=1.1, arrowsize=0.6)
    ax.add_patch(Circle((xc, yc), r, fc="0.35", ec="r", lw=1.2, zorder=5))
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    return cf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--re10", required=True)
    ap.add_argument("--re10-fine")
    ap.add_argument("--re100", required=True)
    ap.add_argument("--re100-hires", help="dir with the Tuolumne plane_z0.csv of ibm_sphere_re100/hires")
    ap.add_argument("--settle", required=True)
    ap.add_argument("--settle-log", required=True)
    ap.add_argument("--out", default="ibm_slices.png")
    a = ap.parse_args()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # final body position of the settling sphere from the log
    yb = None
    for l in open(a.settle_log, errors="ignore"):
        m = re.search(r"^settle: step \d+\s+t\s+[\d.eE+-]+\s+y\s+([\d.eE+-]+)\s+v\s+([\d.eE+-]+)", l)
        if m:
            yb, vb = float(m.group(1)), float(m.group(2))

    widths = [1, 1, 1.4] + ([1.4] if a.re100_hires else []) + [1.1]
    fig, axs = plt.subplots(1, len(widths), figsize=(4 * len(widths) + 1, 5.2), gridspec_kw={"width_ratios": widths})
    cf = panel(axs[0], load(a.re10), 0, 0, 0.005, "Re 10, FWHM 0.08 d (2 GLL spacings)\n$C_d$ 4.69, |u|/U")
    if a.re10_fine:
        panel(axs[1], load(a.re10_fine), 0, 0, 0.005, "Re 10, FWHM 0.04 d (1 GLL spacing)\n$C_d$ 4.09, |u|/U")
    panel(axs[2], load(a.re100), 0, 0, 0.005, "Re 100, N = 5, FWHM 0.12 d, $C_d$ 1.29, |u|/U", zoom=(-0.03, 0.03, -0.03, 0.06))
    if a.re100_hires:
        panel(axs[3], load(a.re100_hires), 0, 0, 0.005, "Re 100 hires (Tuolumne), N = 7, FWHM 0.04 d\n$C_d$ 1.20, |u|/U",
              zoom=(-0.03, 0.03, -0.03, 0.06))
    d = load(a.settle)
    # frame moving with the sphere: subtract the body velocity so the wake shows
    d = dict(d)
    d["v"] = d["v"] - vb
    panel(axs[-1], d, 0, yb, 1 / 12, f"settling sphere, t = 1, Ga 49\nsphere frame, |u - v_b|/|v_b|, v_b = {vb:.3f}",
          zoom=(-0.6, 0.6, yb - 0.7, yb + 0.5), uref=abs(vb))
    for ax in axs:
        ax.set_xlabel("x")
    axs[0].set_ylabel("y (flow direction)")
    fig.colorbar(cf, ax=axs, shrink=0.8, label="velocity magnitude / reference")
    fig.savefig(a.out, dpi=130, bbox_inches="tight")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
