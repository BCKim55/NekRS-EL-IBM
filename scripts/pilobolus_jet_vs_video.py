"""Stack the Video S3 tracking panel (sporangium_tracking.mp4 from pilobolus_track.py/make_video.py) over the two simulations in the same lab window.
Experiment frames: 1100x792 tracking video, image panel = x -4..12.5 mm, y -1..3 mm at
cols 248..880, rows 96..253; frame i is t = (i-15)*20 us. Simulation snapshots every 8 us."""
import glob, os, subprocess, sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont
S = sys.argv[1]
W = 1268                                     # composite width = sim frame width
exp = sorted(glob.glob(os.path.join(S, "exp", "e_*.png")))
sims = {k: sorted(glob.glob(os.path.join(S, "frames_%s" % k, "f_*.png"))) for k in ("bridges", "nobridges")}
font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 22)
out = os.path.join(S, "cmp"); os.makedirs(out, exist_ok=True)
for f in glob.glob(os.path.join(out, "*.png")): os.remove(f)
n = 0
for i, ef in enumerate(exp):
    t = (i - 15) * 20e-6
    if t < 0 or t > 1.25e-3: continue
    e = Image.open(ef).crop((248, 96, 881, 254)).resize((W, 316), Image.BICUBIC)
    d = ImageDraw.Draw(e); d.text((12, 8), "experiment (Video S3)   t = %.3f ms" % (t * 1e3), fill=(200, 0, 0), font=font)
    rows = [e]
    for k in ("bridges", "nobridges"):
        j = min(int(round(t / 8e-6)), len(sims[k]) - 1)
        if t > 8e-6 * (len(sims[k]) - 1) + 4e-6:
            im = Image.new("RGB", (W, 316), (40, 40, 40)); ImageDraw.Draw(im).text((12, 8), "simulation, bridges %s: run ends at %.2f ms" % ("on" if k == "bridges" else "off", 8e-6 * (len(sims[k]) - 1) * 1e3), fill=(220, 220, 220), font=font)
        else:
            im = Image.open(sims[k][j]).convert("RGB")
        rows.append(im)
    c = Image.new("RGB", (W, 316 * 3 + 8), (255, 255, 255))
    for r, im in enumerate(rows): c.paste(im, (0, r * 320))
    c.save(os.path.join(out, "c_%04d.png" % n)); n += 1
print("frames", n)
subprocess.check_call(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "12", "-i", os.path.join(out, "c_%04d.png"),
                       "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", os.path.join(S, "pilobolus_jet_vs_videoS3.mp4")])
print("wrote", os.path.join(S, "pilobolus_jet_vs_videoS3.mp4"))
