#!/usr/bin/env pvpython
"""Side-view movie of a pilobolus sap-jet run, camera following the body:
z = 0 slice coloured by |u| (nearest earlier checkpoint), EL parcels as black
points (every snapshot), the hemisphere as its IBM markers.

    python3 scripts/snapshots_to_vtk.py <run>      # el.pvd / body.pvd
    pvpython scripts/paraview_jet_movie.py <run> [--out movie.mp4] [--umax 15]
        [--window -5 1.5 -1.6 1.6] [--fps 12] [--title "bridges on"]
"""
import argparse
import glob
import os
import subprocess

from paraview.simple import *  # noqa: F401,F403

ap = argparse.ArgumentParser()
ap.add_argument("run")
ap.add_argument("--out", default=None)
ap.add_argument("--umax", type=float, default=15.0)
ap.add_argument("--window", type=float, nargs=4, default=[-5.0, 1.5, -1.6, 1.6], help="mm relative to the body: x0 x1 y0 y1")
ap.add_argument("--fps", type=int, default=12)
ap.add_argument("--title", default="")
ap.add_argument("--size", type=int, nargs=2, default=[1300, 640])
a = ap.parse_args()
run = os.path.abspath(a.run)
snap = os.path.join(run, "snapshots")
out = a.out or os.path.join(run, "jet_movie.mp4")
frames = os.path.join(snap, "jet_frames")
os.makedirs(frames, exist_ok=True)
for f in glob.glob(os.path.join(frames, "*.png")):
    os.remove(f)

nek = OpenDataFile(glob.glob(os.path.join(run, "*.nek5000"))[0])
nek.PointArrays = ["velocity"]
el = OpenDataFile(os.path.join(snap, "el.pvd"))
body = OpenDataFile(os.path.join(snap, "body.pvd"))
bodyT = el.TimestepValues  # body.csv has the same steps
# body x(t) from body.csv: step, t, x, ...
import csv

bx = {}
with open(os.path.join(snap, "body.csv")) as fh:
    rd = csv.reader(fh)
    next(rd)
    for row in rd:
        bx[float(row[1])] = float(row[2])

view = CreateView("RenderView")
view.ViewSize = a.size
view.Background = [1, 1, 1]
view.OrientationAxesVisibility = 0
sl = Slice(Input=nek)
sl.SliceType = "Plane"
sl.SliceType.Normal = [0, 0, 1]
sl.SliceType.Origin = [0, 0, 1e-7]
sl.Triangulatetheslice = 0
calc = Calculator(Input=sl)
calc.ResultArrayName = "umag"
calc.Function = "mag(velocity)"
d = Show(calc, view)
ColorBy(d, ("POINTS", "umag"))
lut = GetColorTransferFunction("umag")
lut.ApplyPreset("Viridis (matplotlib)", True)
lut.RescaleTransferFunction(0.0, a.umax)
bar = GetScalarBar(lut, view)
bar.Title = "|u|  [m/s]"
bar.ComponentTitle = ""
bar.TitleColor = [0, 0, 0]
bar.LabelColor = [0, 0, 0]
bar.WindowLocation = "Lower Right Corner"
bar.ScalarBarLength = 0.3
d.SetScalarBarVisibility(view, True)
e = Show(el, view)
e.SetRepresentationType("Points")
e.ColorArrayName = ["POINTS", ""]
e.AmbientColor = [0, 0, 0]
e.DiffuseColor = [0, 0, 0]
e.PointSize = 2.5
b = Show(body, view)
b.SetRepresentationType("Points")
b.ColorArrayName = ["POINTS", ""]
b.AmbientColor = [0.85, 0.25, 0.1]
b.DiffuseColor = [0.85, 0.25, 0.1]
b.PointSize = 3
txt = Text(Text="")
td = Show(txt, view)
td.Color = [0, 0, 0]
td.FontSize = 22
td.WindowLocation = "Upper Left Corner"
view.CameraParallelProjection = 1
view.CameraViewUp = [0, 1, 0]
x0, x1, y0, y1 = [w * 1e-3 for w in a.window]
scale = 0.5 * (y1 - y0) * max(1.0, (x1 - x0) / (y1 - y0) * a.size[1] / a.size[0])

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
nekT = list(nek.TimestepValues)
for i, t in enumerate(el.TimestepValues):
    scene.AnimationTime = t
    tn = max([s for s in nekT if s <= t + 1e-12] or [nekT[0]])
    nek.UpdatePipeline(tn)
    calc.UpdatePipeline(tn)
    xb = bx.get(t, min(bx, key=lambda k: abs(k - t)) and bx[min(bx, key=lambda k: abs(k - t))])
    cx, cy = xb + 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    view.CameraFocalPoint = [cx, cy, 0]
    view.CameraPosition = [cx, cy, 0.05]
    view.CameraParallelProjection = 1
    view.CameraParallelScale = scale  # Show/Render reset it on the first frame
    txt.Text = "%s   t = %.3f ms   x_body = %.2f mm   fluid at t = %.2f ms" % (a.title, t * 1e3, xb * 1e3, tn * 1e3)
    Render(view)
    SaveScreenshot(os.path.join(frames, "f_%04d.png" % i), view, ImageResolution=a.size)
print("frames:", i + 1)
subprocess.check_call(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(a.fps), "-i", os.path.join(frames, "f_%04d.png"),
                       "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", out])
print("wrote", out)
