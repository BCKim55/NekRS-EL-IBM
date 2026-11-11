#!/usr/bin/env python3
"""Video S3 tracking panel stacked over sap-jet runs drawn from their EL/body
snapshots alone (no fluid checkpoints needed; matplotlib, no ParaView).

  python3 scripts/pilobolus_jet_parcel_movie.py --video sporangium_tracking.mp4 \
      --run "label=case dir" [--run ...] --out movie.mp4 [--dt 40e-6] [--window -4 12.5 -1 3]

Experiment frames: the 1100x792 tracking video, image panel = x -4..12.5 mm,
y -1..3 mm at cols 248..880, rows 96..253; frame i is t = (i-15)*20 us.
Run frames every --dt from the snapshots (written every 8 us), parcels
coloured by u, the body as its hemisphere outline (dome towards +x).
"""
import argparse
import glob
import os
import shutil
import subprocess
import tempfile

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

ap = argparse.ArgumentParser()
ap.add_argument("--video", required=True)
ap.add_argument("--run", action="append", required=True, help="label=case dir (holds snapshots/)")
ap.add_argument("--out", required=True)
ap.add_argument("--dt", type=float, default=40e-6)
ap.add_argument("--tend", type=float, default=1.25e-3)
ap.add_argument("--window", type=float, nargs=4, default=[-4.0, 12.5, -1.0, 3.0], help="lab window in mm: x0 x1 y0 y1")
ap.add_argument("--umax", type=float, default=15.0)
ap.add_argument("--fps", type=int, default=8)
args = ap.parse_args()

W, H = 1268, 316
x0, x1, y0, y1 = args.window
tmp = tempfile.mkdtemp(prefix="jetmovie_")
exp = os.path.join(tmp, "exp")
os.makedirs(exp)
subprocess.run(["ffmpeg", "-v", "error", "-i", args.video, os.path.join(exp, "e_%04d.png")], check=True)
expf = sorted(glob.glob(os.path.join(exp, "e_*.png")))
font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 22)


def load_run(d):
    body = np.genfromtxt(os.path.join(d, "snapshots", "body.csv"), delimiter=",", names=True)
    els = {}
    for f in glob.glob(os.path.join(d, "snapshots", "el_*.csv")):
        els[int(os.path.basename(f)[3:8])] = f
    return body, els


runs = [(s.split("=", 1)[0], load_run(s.split("=", 1)[1])) for s in args.run]
dpi = 100
times = np.arange(0, args.tend + 1e-9, args.dt)
n = 0
for t in times:
    i = int(round(t / 20e-6)) + 15
    if i >= len(expf):
        break
    e = Image.open(expf[i]).crop((248, 96, 881, 254)).resize((W, H), Image.BICUBIC)
    d = ImageDraw.Draw(e)
    d.text((12, 8), "experiment (Video S3)   t = %.3f ms" % (t * 1e3), fill=(200, 0, 0), font=font)
    rows = [e]
    for label, (body, els) in runs:
        step = int(round(t / 0.8e-6))
        k = min(els, key=lambda s: abs(s - step))
        b = body[np.argmin(np.abs(body["step"] - k))]
        fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi)
        ax = fig.add_axes([0, 0, 1, 1])
        a = np.genfromtxt(els[k], delimiter=",", names=True)
        if a.size:
            a = np.atleast_1d(a)
            ax.scatter(a["x"] * 1e3, a["y"] * 1e3, c=a["u"], s=1.5, vmin=0, vmax=args.umax, cmap="Blues", linewidths=0)
        R = b["radius"] * 1e3
        th = np.linspace(-np.pi / 2, np.pi / 2, 60)
        xs = np.concatenate([[b["x"] * 1e3], b["x"] * 1e3 + R * np.cos(th), [b["x"] * 1e3]])
        ys = np.concatenate([[b["y"] * 1e3 - R], b["y"] * 1e3 + R * np.sin(th), [b["y"] * 1e3 + R]])
        ax.fill(xs, ys, color="0.25", lw=0)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.axvline(0, color="0.7", lw=0.5)
        ax.axhline(0, color="0.7", lw=0.5)
        ax.text(0.01, 0.9, "%s   t = %.3f ms   x %.2f mm  v %.2f m/s   parcels %d" % (label, b["t"] * 1e3, b["x"] * 1e3, b["u"], a.size),
                transform=ax.transAxes, fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.canvas.draw()
        im = Image.frombuffer("RGBA", fig.canvas.get_width_height(), fig.canvas.buffer_rgba(), "raw", "RGBA", 0, 1).convert("RGB")
        plt.close(fig)
        rows.append(im.resize((W, H)))
    c = Image.new("RGB", (W, (H + 4) * len(rows)), (255, 255, 255))
    for r, im in enumerate(rows):
        c.paste(im, (0, r * (H + 4)))
    c.save(os.path.join(tmp, "c_%04d.png" % n))
    n += 1
subprocess.run(["ffmpeg", "-v", "error", "-y", "-framerate", str(args.fps), "-i", os.path.join(tmp, "c_%04d.png"),
                "-c:v", "libx264", "-pix_fmt", "yuv420p", "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", args.out], check=True)
shutil.rmtree(tmp)
print("frames", n, "->", args.out)
