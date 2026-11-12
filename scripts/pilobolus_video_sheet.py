"""Sap sheet metrics from the Video S3 tracking movie (sporangium_tracking.mp4, make_video.py).

  ffmpeg -i sporangium_tracking.mp4 frames/f_%03d.png
  python3 scripts/pilobolus_video_sheet.py frames trajectory.csv --out docs/data/pilobolus_video_sheet.csv

The movie's top panel shows Video S3 (50 000 fps, 65 um/px) in a lab window x -4..12.5, y -1..3 mm
about the sporangium's first tracked position (frame 123 = t 0). Per frame: background-subtracted dark
pixels (pre-launch mean), the tracked sporangium (r 0.3 mm) and the thread (|y| < 0.15 mm) removed,
columns with > 0.4 mm vertical extent kept as the sheet. Reports the sheet's darkness-weighted x, its
3rd/97th y percentiles (the lower rim leaves the field at y -0.9), the forward lean of the arms
(|y| > 0.6) relative to the core, and the pixel count (drops when the sheet breaks into droplets).
"""
import argparse, glob, os
import numpy as np
from PIL import Image

ap = argparse.ArgumentParser()
ap.add_argument("frames")
ap.add_argument("traj")
ap.add_argument("--out", default="pilobolus_video_sheet.csv")
ap.add_argument("--overlay", default=None, help="montage png of the detected sheet (first 24 frames)")
args = ap.parse_args()
fs = sorted(glob.glob(os.path.join(args.frames, "f_*.png")))
X0, X1, Y0, Y1 = 248, 881, 96, 254  # panel axes box in the movie frame
sx = (X1 - X0) / 16.5
sy = (Y1 - Y0) / 4.0


def gray(i):
    return np.asarray(Image.open(fs[i]).convert("L"), dtype=np.float32)[Y0:Y1, X0:X1]


bg = np.mean([gray(i) for i in range(2, 12)], axis=0)  # t < -100 us
traj = np.genfromtxt(args.traj, delimiter=",", names=True)
H, W = bg.shape
ys, xs = np.mgrid[0:H, 0:W]
xmm = -4 + xs / sx
ymm = 3 - ys / sy
rows, over = [], []
for i in range(15, 50):  # movie frame i = video frame 108 + i
    fr = 108 + i
    t = (fr - 123) * 20e-6
    d = bg - gray(i)
    m = d > 25
    k = np.argmin(np.abs(traj["frame"] - fr))
    xs_, ys_ = (traj["x_mm"][k], traj["y_mm"][k]) if abs(traj["frame"][k] - fr) < 0.5 else (np.nan, np.nan)
    m &= (xmm > -0.3) & (xmm < 4.0)
    if not np.isnan(xs_):
        m &= ((xmm - xs_) ** 2 + (ymm - ys_) ** 2) > 0.3**2
    m &= np.abs(ymm) > 0.15
    col = np.zeros(W, bool)
    for c in range(W):
        yy = ymm[m[:, c], c]
        if yy.size and (yy.max() - yy.min()) > 0.4:
            col[c] = True
    ms = m & col[None, :]
    over.append((i, ms))
    if ms.sum() < 5:
        rows.append((fr, t, xs_, np.nan, np.nan, np.nan, np.nan, np.nan, 0))
        continue
    w, X, Y = d[ms], xmm[ms], ymm[ms]
    xc = (w * X).sum() / w.sum()
    ytop, ybot = np.percentile(Y, 97), np.percentile(Y, 3)
    arm, core = np.abs(Y) > 0.6, np.abs(Y) < 0.35
    lean = (X[arm].mean() if arm.any() else np.nan) - (X[core].mean() if core.any() else np.nan)
    rows.append((fr, t, xs_, xc, ytop, ybot, ytop - ybot, lean, int(ms.sum())))
np.savetxt(args.out, np.array(rows, dtype=float), delimiter=",",
           header="frame,t_s,x_spor_mm,x_sheet_mm,y_top_mm,y_bot_mm,height_mm,arm_lean_mm,npix", comments="")
print("frame  t_us  x_spor  x_sheet  y_top  y_bot  height  lean  npix")
for r in rows:
    print("%d %5.0f %7.2f %7.2f %6.2f %6.2f %6.2f %6.2f %5d" % (r[0], r[1] * 1e6, *r[2:]))
if args.overlay:
    tiles = []
    for i, ms in over[:24]:
        px = np.asarray(Image.open(fs[i]).convert("RGB").crop((X0, Y0 - 36, X1, Y1))).copy()
        px[36:][ms] = [255, 0, 0]
        tiles.append(Image.fromarray(px))
    w, h = tiles[0].size
    M = Image.new("RGB", (4 * w, 6 * h), "white")
    for k, tl in enumerate(tiles):
        M.paste(tl, ((k % 4) * w, (k // 4) * h))
    M.save(args.overlay)
