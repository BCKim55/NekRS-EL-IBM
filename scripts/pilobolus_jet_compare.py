#!/usr/bin/env python3
"""Compare pilobolus sap-jet runs (e.g. liquid bridges on/off).

Each run is a directory holding the nekRS log (*.out) and snapshots/ written by
elibm::writeSnapshot.  Produces a 2x3 summary figure and a side/rear-view gif of
the parcels and the body.

  python3 scripts/pilobolus_jet_compare.py --run label=dir [--run label=dir ...]
      --png docs/img/pilobolus_jet_compare.png --gif docs/img/pilobolus_jet.gif
"""
import argparse
import glob
import io
import os
import re

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

LINE = re.compile(r"^launch: step (\d+)\s+t (\S+)\s+x (\S+)\s+v (\S+)\s+a (\S+)\s+Fm (\S+) (\S+) (\S+)\s+F (\S+)\s+Cd (\S+)"
                  r"\s+a_f (\S+)\s+slip (\S+)\s+drops (\d+)\s+stuck (\d+)\s+m_stuck (\S+)\s+P_sap (\S+)")
COLS = ["step", "t", "x", "v", "a", "fmx", "fmy", "fmz", "F", "Cd", "a_f", "slip", "drops", "stuck", "m_stuck", "P_sap"]


def read_log(d):
    rows = []
    for f in glob.glob(os.path.join(d, "*.out")):
        with open(f, errors="replace") as fh:
            for line in fh:
                m = LINE.match(line)
                if m:
                    rows.append([float(g) for g in m.groups()])
    a = np.array(sorted(rows))
    return {c: a[:, i] for i, c in enumerate(COLS)}


def read_body(d):
    b = np.loadtxt(os.path.join(d, "snapshots", "body.csv"), delimiter=",", skiprows=1)
    return {"step": b[:, 0].astype(int), "t": b[:, 1], "x": b[:, 2], "r": b[:, 8]}


def snapshots(d):
    out = {}
    for f in sorted(glob.glob(os.path.join(d, "snapshots", "el_*.csv"))):
        step = int(os.path.basename(f)[3:8])
        a = np.loadtxt(f, delimiter=",", skiprows=1, ndmin=2)
        if a.shape[1] >= 7:  # skip empty snapshots (before injection)
            out[step] = a
    return out


def summary(runs, png):
    fig, ax = plt.subplots(2, 3, figsize=(14, 7.5))
    ax = ax.ravel()
    for label, r in runs.items():
        L, ms = r["log"], r["log"]["t"] * 1e3
        ax[0].plot(ms, L["Cd"], label=label)
        ax[1].plot(ms, -L["a_f"], label=label)
        ax[2].plot(ms, L["slip"], label=label)
        ax[3].plot(ms, L["P_sap"] * 1e9, label=label)
        # parcel x relative to the body, last common snapshot
        snaps, body = r["snap"], r["body"]
        for s, style in zip(r["pick"], ["-", "--", ":"]):
            if s in snaps:
                xb = np.interp(s, body["step"], body["x"])
                dx = (snaps[s][:, 1] - xb) * 1e3
                ax[4].hist(dx, bins=np.linspace(-6, 0.5, 66), histtype="step", linestyle=style,
                           label=f"{label}, t = {np.interp(s, body['step'], body['t'])*1e3:.2f} ms")
                u = snaps[s][:, 4]
                ax[5].scatter(dx, u, s=1, alpha=0.4, label=f"{label}, t = {np.interp(s, body['step'], body['t'])*1e3:.2f} ms")
    ax[0].set_ylabel("C_d (frontal area, body only)")
    ax[1].set_ylabel("-a_f = -F_x / m_body  [m/s$^2$]")
    ax[2].set_ylabel("slip rms  [m/s]")
    ax[3].set_ylabel("sap x-momentum  [µg m/s]")
    for a in ax[:4]:
        a.set_xlabel("t  [ms]")
        a.grid(alpha=0.3)
        a.legend(fontsize=8)
    ax[4].set_xlabel("x - x_body  [mm]")
    ax[4].set_ylabel("parcels per 0.1 mm")
    ax[5].set_xlabel("x - x_body  [mm]")
    ax[5].set_ylabel("parcel u  [m/s]")
    ax[5].axhline(0, color="k", lw=0.5)
    for a in ax[4:]:
        a.grid(alpha=0.3)
        a.legend(fontsize=7, markerscale=5)
    fig.suptitle("Pilobolus sap jet: liquid bridges on vs off (prescribed body motion, d 282 µm, 40 nL sap as 25 µm parcels)")
    fig.tight_layout()
    fig.savefig(png, dpi=130)
    print("wrote", png)


def gif(runs, path, every=2, fps=12):
    from PIL import Image
    steps = sorted(set.intersection(*[set(r["snap"]) for r in runs.values()]))[::every]
    n = len(runs)
    frames = []
    for s in steps:
        fig, ax = plt.subplots(n, 2, figsize=(12, 3.4 * n), gridspec_kw={"width_ratios": [7, 3]}, squeeze=False,
                               constrained_layout=True)
        for i, (label, r) in enumerate(runs.items()):
            p, body = r["snap"][s], r["body"]
            xb = np.interp(s, body["step"], body["x"]) * 1e3
            R = np.interp(s, body["step"], body["r"]) * 1e3
            t = np.interp(s, body["step"], body["t"]) * 1e3
            u = p[:, 4]
            for a, (c1, c2) in zip(ax[i], [(1, 2), (2, 3)]):
                sc = a.scatter(p[:, c1] * 1e3 - (xb if c1 == 1 else 0), p[:, c2] * 1e3, c=u, s=2, cmap="viridis", vmin=0, vmax=15)
                if c1 == 1:  # side view: hemisphere, flat face aft
                    th = np.linspace(-np.pi / 2, np.pi / 2, 60)
                    a.fill(np.r_[R * np.cos(th), 0], np.r_[R * np.sin(th), -R], color="0.3", zorder=3)
                    a.set_xlim(-6, 1)
                    a.set_ylim(-1.5, 1.5)
                    if i == n - 1:
                        a.set_xlabel("x - x_body  [mm]")
                    a.set_ylabel("y  [mm]")
                    a.set_title(f"{label}:  t = {t:.3f} ms,  x_body = {xb:.2f} mm,  {len(p)} parcels", loc="left", fontsize=10)
                else:
                    a.add_patch(plt.Circle((0, 0), R, color="0.3", fill=False, lw=1.5, zorder=3))
                    a.set_xlim(-1.5, 1.5)
                    a.set_ylim(-1.5, 1.5)
                    a.set_xlabel("y  [mm]")
                    a.set_ylabel("z  [mm]")
                a.set_aspect("equal")
                a.grid(alpha=0.2)
        fig.colorbar(sc, ax=ax[:, 1].tolist(), label="parcel u  [m/s]", shrink=0.6)
        buf = io.BytesIO()
        fig.savefig(buf, dpi=72, format="png")
        plt.close(fig)
        buf.seek(0)
        frames.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE))
    frames[0].save(path, save_all=True, append_images=frames[1:], duration=int(1000 / fps), loop=0)
    print("wrote", path, len(frames), "frames")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True, help="label=dir")
    ap.add_argument("--png", default="docs/img/pilobolus_jet_compare.png")
    ap.add_argument("--gif", default="")
    ap.add_argument("--pick", default="600,1200", help="snapshot steps for the histogram panels")
    ap.add_argument("--every", type=int, default=2)
    a = ap.parse_args()
    runs = {}
    for spec in a.run:
        label, d = spec.split("=", 1)
        runs[label] = {"log": read_log(d), "body": read_body(d), "snap": snapshots(d), "pick": [int(s) for s in a.pick.split(",")]}
    summary(runs, a.png)
    if a.gif:
        gif(runs, a.gif, every=a.every)


if __name__ == "__main__":
    main()
