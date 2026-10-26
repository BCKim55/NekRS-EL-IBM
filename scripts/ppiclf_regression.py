#!/usr/bin/env python3
"""ppiclF regression: pristine vs patched library on ppiclF's stokes_2d example.

    python3 scripts/ppiclf_regression.py <dir with out_{pristine,patched}_stokes_2d_np{1,2,3}>

stokes_2d: 250 particles per rank, Stokes drag with tau_p = 1/9.8 and gravity
9.8, from rest, dt = 1e-4, VTU every 100 steps. The exact solution is
    v(t) = -g tau (1 - exp(-t/tau)),   y(t2) - y(t1) = int v dt,
independent of the (random, rank-seeded) initial position, so every particle
of every run can be checked: y04 (v) at the last output against v(t), and
the displacement between the first and last output against the integral.
Particles are matched between outputs by their ppiclF tags.
"""
import glob
import os
import re
import struct
import sys

import numpy as np

G, TAU, DT = 9.8, 1.0 / 9.8, 1e-4


def read_vtu(path):
    raw = open(path, "rb").read()
    head = raw[: raw.index(b"<AppendedData")].decode()
    time = float(re.search(r'Name="TIME"[^>]*>\s*([^<]+)<', head).group(1))
    names = re.findall(r'Name="([^"]+)" NumberOfComponents="(\d+)" format="append" offset="(\d+)"', head)
    start = raw.index(b"_", raw.index(b"<AppendedData")) + 1
    out = {}
    for name, ncomp, off in names:
        o = start + int(off)
        nbytes = struct.unpack("<I", raw[o : o + 4])[0]
        arr = np.frombuffer(raw[o + 4 : o + 4 + nbytes], dtype="<f4")
        out[name.strip()] = arr.reshape(-1, int(ncomp)) if int(ncomp) > 1 else arr
    return time, out


def key(d):
    return list(zip(d["tag01"].astype(int), d["tag02"].astype(int), d["tag03"].astype(int)))


def main():
    base = sys.argv[1]
    for v in ("pristine", "patched"):
        for np_ in (1, 2, 3):
            d = os.path.join(base, f"out_{v}_stokes_2d_np{np_}")
            files = sorted(glob.glob(os.path.join(d, "par*.vtu")))
            if len(files) < 2:
                print(f"{v} np{np_}: no output")
                continue
            # ppiclF stamps a file with the time passed to IntegrateParticle,
            # i.e. the start of the step; the state in it is one dt later
            t1, a = read_vtu(files[0])
            t2, b = read_vtu(files[-1])
            t1, t2 = t1 + DT, t2 + DT
            ka, kb = key(a), key(b)
            ia = {k: i for i, k in enumerate(ka)}
            common = [k for k in kb if k in ia]
            ya = np.array([a["y02"][ia[k]] for k in common])
            yb = np.array([b["y02"][kb.index(k)] for k in common])
            vb = np.array([b["y04"][kb.index(k)] for k in common])
            v_ex = -G * TAU * (1 - np.exp(-t2 / TAU))
            disp_ex = -G * TAU * ((t2 - t1) - TAU * (np.exp(-t1 / TAU) - np.exp(-t2 / TAU)))
            err_v = np.abs(vb - v_ex).max()
            err_y = np.abs((yb - ya) - disp_ex).max()
            nbad = int((np.abs((yb - ya) - disp_ex) > 1e-5).sum())
            print(f"{v:8s} np{np_}: {len(kb)} particles ({len(common)} matched), t {t1:.4f}->{t2:.4f}, "
                  f"max |v - exact| {err_v:.2e}, max |dy - exact| {err_y:.2e}, particles off by > 1e-5: {nbad}")


if __name__ == "__main__":
    main()
