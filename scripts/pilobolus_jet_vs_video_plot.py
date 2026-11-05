"""Video S3 tracking (trajectory.csv from pilobolus_track.py) against the sap-jet runs:
x(t), speed and acceleration. The body speed in the runs is prescribed from the same
tracking, so the physics is in the acceleration panel: measured deceleration vs the
fluid-force acceleration a_f (air + sap) the runs report.
    python3 scripts/pilobolus_jet_vs_video_plot.py trajectory.csv --run "bridges on=<dir>" --run "bridges off=<dir>" --out fig.png
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pilobolus_jet_compare import read_log

ap = argparse.ArgumentParser()
ap.add_argument("csv")
ap.add_argument("--run", action="append", default=[], help="label=dir (dir holds the *.out log)")
ap.add_argument("--x0", type=float, default=0.25e-3, help="body centre at t = 0 in the runs [m]")
ap.add_argument("--out", default="pilobolus_jet_vs_videoS3.png")
a = ap.parse_args()

E = np.genfromtxt(a.csv, delimiter=",", names=True)
te, xe, ve, ae = E["t_us"] * 1e-3, E["x_mm"], E["speed_mps"], E["accel_g"] * 9.81
k = 5
ae_s = np.convolve(ae, np.ones(k) / k, mode="same")  # 5-frame (100 us) moving average

fig, ax = plt.subplots(1, 3, figsize=(13, 3.8), constrained_layout=True)
ax[0].plot(te, xe, "k.", ms=4, label="Video S3 (tracked centroid)")
ax[1].plot(te, ve, "k.", ms=4, label="Video S3")
ax[2].plot(te, ae, color="0.7", lw=1, label="Video S3, per frame")
ax[2].plot(te, ae_s, "k-", lw=1.5, label="Video S3, 100 µs mean")
for i, r in enumerate(a.run):
    lab, d = r.split("=", 1)
    L = read_log(d)
    t = L["t"] * 1e3
    c = "C%d" % i
    ax[0].plot(t, (L["x"] - a.x0) * 1e3, c, lw=1.5, label="run, %s" % lab)
    ax[1].plot(t, L["v"], c, lw=1.5, label="run (prescribed), %s" % lab)
    ax[2].plot(t, L["a_f"], c, lw=1.5, label="run a_f (air + sap), %s" % lab)
    if i == 0:
        ax[2].plot(t, L["a"], "C3--", lw=1, label="run, prescribed a")
ax[0].set(xlabel="t [ms]", ylabel="x − x(0) [mm]", title="position")
ax[1].set(xlabel="t [ms]", ylabel="speed [m/s]", title="body speed (runs prescribe the tracked speed)")
ax[2].set(xlabel="t [ms]", ylabel="a [m/s²]", title="acceleration", ylim=(-3000, 3000))
ax[2].axhline(0, color="k", lw=0.5)
for x in ax:
    x.grid(alpha=0.3); x.legend(fontsize=7)
fig.savefig(a.out, dpi=130)
# in-flight numbers (0.4 ms on, where the launch transient is over)
m = te >= 0.4
print("Video S3 mean decel 0.4-1.24 ms: %.0f m/s^2" % -ae[m].mean())
for r in a.run:
    lab, d = r.split("=", 1); L = read_log(d); mm = L["t"] * 1e3 >= 0.4
    print("%s: mean a_f %.0f m/s^2 = %.0f %% of the measured decel" % (lab, L["a_f"][mm].mean(), 100 * L["a_f"][mm].mean() / ae[m].mean()))
print("wrote", a.out)
