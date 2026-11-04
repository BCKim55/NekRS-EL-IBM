#!/usr/bin/env python3
"""Sporangium position and speed: simulation body history vs the Video S3 track.

Measured file: Spore/results/trajectory.csv (frame, t_us, x_mm, y_mm, vx, vy,
speed, accel_g) tracked at 15.3 px/mm, t = 0 at frame 123, the first frame in
flight; the cap was at rest in frame 122 and the release is estimated at frame
122.25 (force_analysis/REPORT.md), 2.3 px behind the frame-123 centroid.
Lengths and speeds scale with 15.3/ppm for another calibration (--ppm 20 is
the morphology-based scale of the hand-off).

Simulation: `launch: step` lines of a pilobolus_launch run (x, v, a prescribed;
a_f = F/(rho_b V), the acceleration the fluid force alone would give). The
air-only curve integrates a_f from t_rise as a free flight.
"""
import argparse
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PPM_FILE = 15.3      # calibration of trajectory.csv
FRAME_US = 20.0
FRAME0 = 123         # first tracked frame
RELEASE = 122.25     # release frame estimate
OFFSET_PX = 2.3      # frame-123 centroid displacement from rest
G = 9.81


def measured(path, ppm):
    d = np.genfromtxt(path, delimiter=",", names=True)
    s = PPM_FILE / ppm
    t = (d["t_us"] + (FRAME0 - RELEASE) * FRAME_US) * 1e-6
    x = (d["x_mm"] + OFFSET_PX / PPM_FILE) * s           # mm from rest
    v = d["speed_mps"] * s
    a = d["accel_g"] * G * s
    return t, x, v, a


def sim(path):
    pat = re.compile(r"^launch: step \d+\s+t ([\d.eE+-]+)\s+x ([\d.eE+-]+)\s+v ([\d.eE+-]+)\s+a ([\d.eE+-]+)"
                     r".*?\s+a_f ([\d.eE+-]+)")
    rows = [[float(g) for g in m.groups()] for l in open(path, errors="ignore") if (m := pat.match(l))]
    r = np.array(rows)
    t, x, v, a, af = r.T
    hdr = next((l for l in open(path, errors="ignore") if l.startswith("launch: d ")), "")
    m = re.search(r"t_rise ([\d.eE+-]+)", hdr)
    tr = float(m.group(1)) if m else 2.2e-4
    x = (x - x[0] + v[0] * t[0]) * 1e3   # mm from the start position
    # air-only free flight after t_rise: v' = a_f, x' = v
    i0 = np.searchsorted(t, tr)
    va, xa = v.copy(), x.copy()
    for i in range(i0 + 1, len(t)):
        dt = t[i] - t[i - 1]
        va[i] = va[i - 1] + 0.5 * (af[i] + af[i - 1]) * dt
        xa[i] = xa[i - 1] + 0.5 * (va[i] + va[i - 1]) * dt * 1e3
    return t, x, v, a, af, va, xa, tr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--measured", required=True)
    ap.add_argument("--sim", action="append", default=[], help="log[:label], repeatable")
    ap.add_argument("--ppm", type=float, nargs="*", default=[15.3, 20.0], help="calibrations to show")
    ap.add_argument("--out", default="pilobolus_trajectory.png")
    a = ap.parse_args()

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    mk = ["o", "s"]
    for k, ppm in enumerate(a.ppm):
        t, x, v, acc = measured(a.measured, ppm)
        kw = dict(color="k", ms=3.5, lw=0, marker=mk[k % 2], mfc="k" if k == 0 else "none",
                  label=f"video S3, {ppm:g} px/mm")
        ax[0].plot(t * 1e3, x, **kw)
        ax[1].plot(t * 1e3, v, **kw)
        ax[2].plot(t * 1e3, acc * 1e-3, **kw)
        fl = t > 0.3e-3
        p = np.polyfit(t[fl], v[fl], 1)
        print(f"video {ppm:g} px/mm: launch v {v[0]:.2f} at {t[0]*1e6:.0f} us, peak {v.max():.2f}, "
              f"x(1.2 ms) {np.interp(1.2e-3, t, x):.2f} mm, flight fit {p[0]:.0f} m/s2, v(1.2 ms) {np.interp(1.2e-3, t, v):.2f}")
    for k, spec in enumerate(a.sim):
        path, _, label = spec.partition(":")
        label = label or path
        t, x, v, acc, af, va, xa, tr = sim(path)
        c = f"C{k}"
        ax[0].plot(t * 1e3, x, c, lw=1.8, label=f"{label}: prescribed")
        ax[0].plot(t * 1e3, xa, c, lw=1.2, ls="--", label=f"{label}: air only after t_rise")
        ax[1].plot(t * 1e3, v, c, lw=1.8)
        ax[1].plot(t * 1e3, va, c, lw=1.2, ls="--")
        ax[2].plot(t * 1e3, acc * 1e-3, c, lw=1.8)
        ax[2].plot(t * 1e3, af * 1e-3, c, lw=1.2, ls="--")
        fl = t > tr
        print(f"sim {label}: v_max {v.max():.2f} at {t[np.argmax(v)]*1e6:.0f} us, x(1.2 ms) {np.interp(1.2e-3, t, x):.2f} mm, "
              f"v(1.2 ms) {np.interp(1.2e-3, t, v):.2f}; air-only mean a_f {af[fl].mean():.0f} m/s2 -> v(1.2 ms) {va[-1]:.2f}, "
              f"x(1.2 ms) {xa[-1]:.2f} mm")
    ax[0].set_ylabel("position from rest [mm]")
    ax[1].set_ylabel("speed [m/s]")
    ax[2].set_ylabel("acceleration [10³ m/s²]")
    ax[2].set_ylim(-5, 60)
    ax[2].axhline(0, color="0.6", lw=0.6)
    for x_ in ax:
        x_.set_xlabel("t since release [ms]")
        x_.grid(alpha=0.3)
    ax[0].legend(fontsize=8)
    ax[2].text(0.98, 0.95, "measured launch ≤ 40 µs (≥ 2.5e5 m/s², off scale)\nflight decel 1592 ± 103 m/s² (15.3 px/mm)",
               transform=ax[2].transAxes, ha="right", va="top", fontsize=8)
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print("wrote", a.out)


if __name__ == "__main__":
    sys.exit(main())
