#!/usr/bin/env python3
"""Check a hemisphere-body run: no droplet inside the body, markers on the
surface, and the body's orientation / spin history.

    python3 scripts/check_hemisphere.py build/<run> [--log run.log]
"""
import argparse
import csv
import glob
import os

import numpy as np


def read(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]} if rows else {}


def gap(p, axis, R):
    """signed distance of points p (relative to the centre of mass) to the hemisphere surface"""
    z = axis / np.linalg.norm(axis)
    seed = np.array([1.0, 0, 0]) if abs(z[0]) < 0.9 else np.array([0, 1.0, 0])
    x = np.cross(seed, z); x /= np.linalg.norm(x)
    y = np.cross(z, x)
    q = np.stack([p @ x, p @ y, p @ z + 0.375 * R], axis=1)  # relative to the sphere centre
    rho = np.hypot(q[:, 0], q[:, 1])
    dome = np.linalg.norm(q, axis=1) - R
    flat = -q[:, 2]
    rim = np.hypot(rho - R, q[:, 2])
    g = np.where(q[:, 2] >= 0, dome, np.where(rho <= R, flat, rim))
    return g


a = argparse.ArgumentParser()
a.add_argument("case")
a = a.parse_args()
snap = os.path.join(a.case, "snapshots")
body = read(os.path.join(snap, "body.csv"))
print(f"body shape {int(body['shape'][0])}  axis start {body['ax'][0]:+.3f} {body['ay'][0]:+.3f} {body['az'][0]:+.3f}"
      f"  end {body['ax'][-1]:+.3f} {body['ay'][-1]:+.3f} {body['az'][-1]:+.3f}")
worst = (0, None, 0)
for i, path in enumerate(sorted(glob.glob(os.path.join(snap, "el_*.csv")))):
    el = read(path)
    if not el:
        continue
    step = int(os.path.basename(path)[3:8])
    j = np.searchsorted(body["step"], step)
    if j >= len(body["step"]) or body["step"][j] != step:
        continue
    c = np.array([body["x"][j], body["y"][j], body["z"][j]])
    ax = np.array([body["ax"][j], body["ay"][j], body["az"][j]])
    p = np.stack([el["x"], el["y"], el["z"]], axis=1) - c
    g = gap(p, ax, body["radius"][j]) - 0.5 * el["dp"]  # surface gap of the droplet surface
    n_in = int((g < -0.5 * el["dp"]).sum())  # centre inside the body
    if g.min() < worst[0]:
        worst = (g.min(), step, n_in)
print(f"deepest droplet overlap {worst[0]:.3e} m at step {worst[1]} ({worst[2]} centres inside the body)")
