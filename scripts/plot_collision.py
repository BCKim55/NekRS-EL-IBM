#!/usr/bin/env python3
"""Figures for cases/el_ibm_collision from its snapshots/ directory.

    python3 scripts/plot_collision.py cases/el_ibm_collision

Writes collision_summary.png (snapshots + time histories) and
collision.gif (x-y view over time) next to the snapshots.
"""
import csv
import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter


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
            frames.append((step, body["t"][i], read_csv(path), {k: v[i] for k, v in body.items()}))
    return frames, body


def draw(ax, el, b, title):
    ax.clear()
    ax.add_patch(plt.Circle((b["x"], b["y"]), b["radius"], color="0.75", ec="0.3", lw=1.0, zorder=1))
    touching = el["contact"] > 0
    sc = ax.scatter(el["x"], el["y"], c=el["u"], cmap="coolwarm", vmin=-0.5, vmax=0.5, s=6, zorder=2)
    ax.scatter(el["x"][touching], el["y"][touching], facecolors="none", edgecolors="k", s=30, lw=0.8, zorder=3)
    ax.annotate("", xy=(b["x"] + 0.4 * b["u"] + 1e-9, b["y"]), xytext=(b["x"], b["y"]),
                arrowprops=dict(arrowstyle="->", color="k", lw=1.5), zorder=4)
    ax.set_xlim(0, 1)
    ax.set_ylim(0.2, 0.8)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return sc


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else "cases/el_ibm_collision"
    frames, body = load(case)
    out = os.path.join(case, "snapshots")

    # summary: four snapshots + histories
    fig = plt.figure(figsize=(12, 7))
    picks = [frames[int(round(i * (len(frames) - 1) / 3))] for i in range(4)]
    sc = None
    for k, (step, t, el, b) in enumerate(picks):
        ax = fig.add_subplot(2, 4, k + 1)
        sc = draw(ax, el, b, f"t = {t:.2f}  (step {step})")
    cax = fig.add_axes([0.92, 0.55, 0.012, 0.33])
    fig.colorbar(sc, cax=cax, label="EL particle $u_x$")

    t = np.array([f[1] for f in frames])
    ax = fig.add_subplot(2, 2, 3)
    ax.plot(body["t"], body["u"], "k-", label="IBM sphere $u_x$")
    ax.plot(t, [f[2]["u"].mean() for f in frames], "C3-", label="EL mean $u_x$")
    ax.plot(t, [(f[2]["u"] < 0).sum() / len(f[2]["u"]) for f in frames], "C0--", label="EL fraction with $u_x<0$")
    ax.set_xlabel("t")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(2, 2, 4)
    gap = []
    for _, _, el, b in frames:
        d = np.sqrt((el["x"] - b["x"]) ** 2 + (el["y"] - b["y"]) ** 2 + (el["z"] - b["z"]) ** 2)
        gap.append(d - b["radius"])
    gap = np.array(gap)
    ax.plot(t, gap.min(axis=1) / 0.01, "k-", label="min (EL centre to sphere surface) / $d_p$")
    ax.axhline(0.5, color="C3", ls=":", label="touching: $d_p/2$")
    ax.set_ylim(0, 2)
    ax.set_xlabel("t")
    ax.set_ylabel("gap / $d_p$")
    ax2 = ax.twinx()
    ax2.step(t, [(f[2]["contact"] > 0).sum() for f in frames], "C0-", where="mid", alpha=0.6)
    ax2.set_ylabel("EL in contact", color="C0")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)
    fig.suptitle("EL point particles fired at a free IBM sphere (view along z; circled: in contact with the sphere; arrow: sphere velocity)")
    fig.savefig(os.path.join(out, "collision_summary.png"), dpi=130, bbox_inches="tight")

    # animation
    fig, ax = plt.subplots(figsize=(7, 4.5))

    def update(i):
        step, t, el, b = frames[i]
        draw(ax, el, b, f"EL-IBM collision, t = {t:.3f}")

    anim = FuncAnimation(fig, update, frames=len(frames))
    anim.save(os.path.join(out, "collision.gif"), writer=PillowWriter(fps=8), dpi=90)
    print("wrote", os.path.join(out, "collision_summary.png"), "and collision.gif")


if __name__ == "__main__":
    main()
