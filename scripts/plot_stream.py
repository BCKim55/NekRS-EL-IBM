#!/usr/bin/env python3
"""Figures for cases/el_ibm_stream from its snapshots/ directory and run log.

    python3 scripts/plot_stream.py <run dir> [--log run.log] [--gif]

Writes stream_summary.png (four x-y snapshots + histories) and, with --gif,
stream.gif next to the snapshots. The log supplies the cumulative collision
counts and the drag history ("stream: step ..." lines).
"""
import argparse
import csv
import glob
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_csv(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]} if rows else {}


def load(case):
    snap = os.path.join(case, "snapshots")
    body = read_csv(os.path.join(snap, "body.csv"))
    frames = []
    for path in sorted(glob.glob(os.path.join(snap, "el_*.csv"))):
        step = int(os.path.basename(path)[3:8])
        i = np.searchsorted(body["step"], step)
        if i < len(body["step"]) and body["step"][i] == step:
            el = read_csv(path)
            if el:
                frames.append((step, body["t"][i], el, {k: v[i] for k, v in body.items()}))
    return frames


def read_log(path):
    pat = re.compile(r"stream: step (\d+)\s+t ([\d.]+)\s+nEl (\d+)\s+injected (\d+)\s+exited (\d+)\s+contacts (\d+)\s+"
                     r"u<0 now (\d+)\s+touched (\d+)\s+bounced (\d+)\s+F (\S+) (\S+) (\S+)\s+Cd (\S+)")
    keys = ["step", "t", "nEl", "injected", "exited", "contacts", "backwards", "touched", "bounced", "Fx", "Fy", "Fz", "Cd"]
    rows = [[float(v) for v in m.groups()] for m in (pat.search(l) for l in open(path, errors="replace")) if m]
    return {k: np.array([r[i] for r in rows]) for i, k in enumerate(keys)} if rows else {}


def draw(ax, el, b, title, uin, xlim, ylim):
    ax.clear()
    ax.add_patch(plt.Circle((b["x"], b["y"]), b["radius"], color="0.75", ec="0.3", lw=1.0, zorder=1))
    touching = el["contact"] > 0
    sc = ax.scatter(el["x"], el["y"], c=el["u"], cmap="coolwarm", vmin=-uin, vmax=uin, s=5, zorder=2)
    ax.scatter(el["x"][touching], el["y"][touching], facecolors="none", edgecolors="k", s=28, lw=0.8, zorder=3)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--log", default=None)
    ap.add_argument("--uin", type=float, default=0.5)
    ap.add_argument("--xlim", type=float, nargs=2, default=(0.0, 2.0))
    ap.add_argument("--ylim", type=float, nargs=2, default=(0.1, 0.9))
    ap.add_argument("--gif", action="store_true")
    a = ap.parse_args()
    frames = load(a.case)
    log = read_log(a.log) if a.log else {}
    out = os.path.join(a.case, "snapshots")

    fig = plt.figure(figsize=(13, 8))
    picks = [frames[int(round(i * (len(frames) - 1) / 3))] for i in range(4)]
    sc = None
    for k, (step, t, el, b) in enumerate(picks):
        ax = fig.add_subplot(3, 2, k + 1)
        sc = draw(ax, el, b, f"t = {t:.2f} (step {step}): {len(el['x'])} particles, {int((el['contact'] > 0).sum())} in contact",
                  a.uin, a.xlim, a.ylim)
    cax = fig.add_axes([0.92, 0.42, 0.012, 0.5])
    fig.colorbar(sc, cax=cax, label="particle $u_x$")

    t = np.array([f[1] for f in frames])
    ax = fig.add_subplot(3, 2, 5)
    ax.plot(t, [len(f[2]["x"]) for f in frames], "k-", label="particles in the domain")
    ax.plot(t, [(f[2]["u"] < 0).sum() for f in frames], "C0-", label="moving upstream ($u_x<0$)")
    ax.plot(t, [(f[2]["contact"] > 0).sum() for f in frames], "C3-", label="in contact (end of step)")
    if log:
        ax.plot(log["t"], log["bounced"], "C2--", label="collided so far (distinct particles sent upstream)")
        ax.plot(log["t"], log["exited"], "C1:", label="left through the outflow")
    ax.set_xlabel("t")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(3, 2, 6)
    if log:
        ax.plot(log["t"], log["Cd"], "k-", label="sphere $C_d$ (hydrodynamic force on $u_{in}$, D)")
        ax.set_ylabel("$C_d$")
        ax.legend(fontsize=7, loc="upper left")
    ax2 = ax.twinx()
    ax2.plot(t, [f[2]["u"].mean() / a.uin for f in frames], "C3-", alpha=0.7)
    ax2.set_ylabel("mean particle $u_x / u_{in}$", color="C3")
    ax.set_xlabel("t")
    ax.grid(alpha=0.3)
    fig.suptitle("EL particles ($d_p$ = D/30) streaming at a held IBM sphere: soft-sphere collisions, no adhesion\n"
                 "(x-y view of the mid-plane sheet; circled: in contact at the end of the step)")
    fig.savefig(os.path.join(out, "stream_summary.png"), dpi=130, bbox_inches="tight")
    print("wrote", os.path.join(out, "stream_summary.png"))

    if a.gif:
        from matplotlib.animation import FuncAnimation, PillowWriter

        fig, ax = plt.subplots(figsize=(9, 4))

        def update(i):
            step, tt, el, b = frames[i]
            draw(ax, el, b, f"t = {tt:.3f}  particles {len(el['x'])}  in contact {int((el['contact'] > 0).sum())}",
                 a.uin, a.xlim, a.ylim)

        anim = FuncAnimation(fig, update, frames=len(frames))
        anim.save(os.path.join(out, "stream.gif"), writer=PillowWriter(fps=10), dpi=80)
        print("wrote", os.path.join(out, "stream.gif"))


if __name__ == "__main__":
    main()
