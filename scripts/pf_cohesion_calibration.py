#!/usr/bin/env python3
"""Static lattice sums for the pairwise-force cohesion bond (cohesion_model 2).

F = s H f(r/H), f = 1 + a t^2 - 2(2+a) t^3 + (3+a) t^4, a = -23.5 (Whebell et
al. 2024), attractive part only (the spring is the core), H = d/t0 so the
zero crossing sits at contact. Prints, per unit s and with d = 1:
  - the tension of a single-file chain at contact spacing (target pi d sigma/2),
  - the edge tension of a hex monolayer at contact spacing (target 2 sigma),
and the s/sigma each target implies. Whebell's 3D Laplace calibration
(sigma = 30.96 s, FCC at 1.1 dx) is quoted for comparison.
"""
import numpy as np
from scipy.optimize import brentq

A = -23.5
f = lambda t: 1 + A * t**2 - 2 * (2 + A) * t**3 + (3 + A) * t**4
t0 = brentq(f, 0.1, 0.5)
H = 1 / t0
ts = np.linspace(t0, 1, 4001)


def F(r):  # attractive branch only, per unit s, d = 1; negative = attraction
    t = r / H
    return np.where((t > t0) & (t < 1), H * f(t), 0.0)


k = np.arange(1, 20)
chain = -(k * F(k * 1.0)).sum()  # bonds of order k cross a cut k times

pts = np.array([(i + 0.5 * j, j * np.sqrt(3) / 2) for i in range(-12, 13) for j in range(-12, 13)])
c = 0.25 * np.sqrt(3) / 2
below = pts[(pts[:, 1] < c) & (pts[:, 0] >= -0.5) & (pts[:, 0] < 0.5)]
above = pts[pts[:, 1] > c]
mono = 0.0
for p in below:
    d = above - p
    r = np.hypot(d[:, 0], d[:, 1])
    mono += (-F(r) * d[:, 1] / r).sum()

print(f"zero crossing t0 = {t0:.5f}  H = {H:.3f} d  max |f| = {-f(ts).min():.3f} at t = {ts[f(ts).argmin()]:.3f}")
print(f"max pair attraction        = {-f(ts).min() * H:.3f} s d")
print(f"chain tension at contact   = {chain:.3f} s d   -> pi/2 sigma d needs s = {np.pi / 2 / chain:.3f} sigma")
print(f"hex monolayer edge tension = {mono:.3f} s     -> 2 sigma needs s = {2 / mono:.4f} sigma")
print(f"3D bulk (Whebell 2024)     : sigma = 30.96 s   -> s = {1 / 30.96:.4f} sigma")
