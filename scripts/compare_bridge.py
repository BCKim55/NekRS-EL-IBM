#!/usr/bin/env python3
"""Liquid bridge on/off comparison for cases/el_ibm_bridge.

    python3 scripts/compare_bridge.py build/el_ibm_bridge build/el_ibm_bridge_ref --out fig.png

Left: the bridge force law F(S) of the run (Willett et al. 2000) for the
droplet-droplet and droplet-body pairs. Right: histories from the two run
logs (contacts, collided, body velocity, C_d) and, from the last snapshot,
the droplet distribution round the body.
"""
import argparse
import glob
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from plot_stream import read_csv, read_log  # noqa: E402


def par(case, key, default=None):
    for l in open(os.path.join(case, "cyl.par")):
        m = re.match(rf"\s*{key}\s*=\s*(\S+)", l)
        if m:
            return float(m.group(1))
    return default


def body_u(path):
    pat = re.compile(r"stream: step (\d+)\s+t ([\d.]+).* body x (\S+) (\S+) u (\S+) (\S+)")
    rows = [[float(v) for v in m.groups()] for m in (pat.search(l) for l in open(path, errors="replace")) if m]
    return np.array(rows).T  # step, t, x, y, u, v


def bridge_force(S, r1, r2, sigma, theta, vb):
    R = 2 * r1 * r2 / (r1 + r2)
    srup = (1 + 0.5 * theta) * vb ** (1 / 3)
    sh = 0.5 * S * np.sqrt(R / vb)
    F = 2 * np.pi * sigma * R * np.cos(theta) / (1 + 1.05 * sh + 2.5 * sh**2)
    return np.where(S < srup, F, 0.0), srup


a = argparse.ArgumentParser()
a.add_argument("bridge")
a.add_argument("ref")
a.add_argument("--out", default=None)
a.add_argument("--scale", type=float, default=1000.0, help="length scale factor for the snapshot panel (m -> mm)")
a = a.parse_args()

c = a.bridge
dp, R, sigma = par(c, "dp"), par(c, "body_radius"), par(c, "sigma")
theta = np.radians(par(c, "contact_angle_deg", 0))
vb = par(c, "bridge_volume_fraction", 0.1) * np.pi / 6 * dp**3
rhop = par(c, "rhop")
m = rhop * np.pi / 6 * dp**3

fig, axs = plt.subplots(2, 3, figsize=(15, 8.5))
fig.suptitle("Liquid bridge (Willett et al. 2000) on / off: water droplets $d_p$ = D/30 at a free IBM cylinder, air 0.5 m/s",
             fontsize=11)

ax = axs[0, 0]
S = np.linspace(0, 0.6 * dp, 400)
for r2, lab in [(0.5 * dp, "droplet-droplet"), (R, "droplet-body")]:
    F, srup = bridge_force(S, 0.5 * dp, r2, sigma, theta, vb)
    ax.plot(S / dp, F * 1e6, label=f"{lab}: $F_{{max}}$ {F[0]*1e6:.1f} µN, rupture {srup/dp:.2f} $d_p$")
ax.axhline(m * 9.81 * 1e6, color="0.5", ls=":", label=f"droplet weight {m*9.81*1e6:.4f} µN")
ax.set_xlabel("gap $S / d_p$")
ax.set_ylabel("bridge force [µN]")
ax.set_title(f"force law: σ {sigma} N/m, θ {np.degrees(theta):.0f}°, $V_b/V_p$ {par(c, 'bridge_volume_fraction', 0.1)}", fontsize=9)
ax.legend(fontsize=7)
ax.grid(alpha=0.3)

runs = [(a.bridge, "bridge on", "C3"), (a.ref, "bridge off", "C0")]
logs = {n: read_log(n + ".log") for n, _, _ in runs}
bods = {n: body_u(n + ".log") for n, _, _ in runs}

for ax, key, lab in [(axs[0, 1], "contacts", "droplets in contact / bridged (end of step)"),
                     (axs[0, 2], "bounced", "distinct droplets sent upstream (collided)")]:
    for n, name, col in runs:
        ax.plot(logs[n]["t"], logs[n][key], color=col, label=name)
    ax.set_title(lab, fontsize=9)
    ax.set_xlabel("t [s]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

ax = axs[1, 0]
for n, name, col in runs:
    ax.plot(bods[n][1], bods[n][4], color=col, label=name)
ax.set_title("cylinder velocity $u_x$ [m/s] (droplet impulses + drag + bridge pull)", fontsize=9)
ax.set_xlabel("t [s]")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

ax = axs[1, 1]
for n, name, col in runs:
    t, cd = logs[n]["t"], logs[n]["Cd"]
    k = t > 0.02
    ax.plot(t[k], cd[k], color=col, label=name)
ax.set_title("$C_d$ (hydrodynamic force on $u_{in}$, D, L)", fontsize=9)
ax.set_xlabel("t [s]")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

ax = axs[1, 2]
for n, name, col in runs:
    last = sorted(glob.glob(os.path.join(n, "snapshots", "el_*.csv")))[-1]
    el = read_csv(last)
    b = read_csv(os.path.join(n, "snapshots", "body.csv"))
    bx, by = b["x"][-1], b["y"][-1]
    ax.scatter((el["x"] - bx) * a.scale, (el["y"] - by) * a.scale, s=6, color=col, alpha=0.6, label=f"{name}: {len(el['x'])} droplets")
ax.add_patch(plt.Circle((0, 0), R * a.scale, fill=False, ec="k", lw=1))
ax.set_xlim(-3 * R * a.scale, 5 * R * a.scale)
ax.set_ylim(-3 * R * a.scale, 3 * R * a.scale)
ax.set_aspect("equal")
ax.set_title("droplets round the cylinder at the end (body frame, mm)", fontsize=9)
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

fig.tight_layout()
out = a.out or os.path.join(a.bridge, "snapshots", "bridge_compare.png")
fig.savefig(out, dpi=130)
print("wrote", out)
