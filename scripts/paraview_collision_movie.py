#!/usr/bin/env pvpython
"""Movie of cases/el_ibm_collision: domain outline, midplane slice coloured by
the particle volume fraction phi_p (checkpoint scalar00 = reader array s1),
EL particles and the IBM sphere as sphere glyphs.

    python3 scripts/snapshots_to_vtk.py cases/el_ibm_collision
    /Applications/ParaView-5.10.1.app/Contents/bin/pvpython scripts/paraview_collision_movie.py \
        cases/el_ibm_collision [--out snapshots/collision_paraview.mp4] [--fps 12]

Frames go to <case>/snapshots/frames/ and ffmpeg assembles the mp4. The
particle snapshots are every 4 steps, the fluid checkpoints every 40, so the
slice updates once per 10 frames.
"""
import argparse
import glob
import os
import subprocess

from paraview.simple import *  # noqa: F401,F403

ap = argparse.ArgumentParser()
ap.add_argument("case", nargs="?", default="cases/el_ibm_collision")
ap.add_argument("--out", default=None)
ap.add_argument("--fps", type=int, default=12)
a = ap.parse_args()
case = os.path.abspath(a.case)
snap = os.path.join(case, "snapshots")
out = a.out or os.path.join(snap, "collision_paraview.mp4")
frames = os.path.join(snap, "frames")
os.makedirs(frames, exist_ok=True)
for f in glob.glob(os.path.join(frames, "frame_*.png")):
    os.remove(f)

nek = OpenDataFile(os.path.join(case, "collision.nek5000"))
nek.PointArrays = ["velocity", "s1"]
el = OpenDataFile(os.path.join(snap, "el.pvd"))
body = OpenDataFile(os.path.join(snap, "body.pvd"))

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1280, 960]
view.Background = [1, 1, 1]
try:
    view.UseColorPaletteForBackground = 0
except Exception:
    pass
view.OrientationAxesVisibility = 1

# whole flow domain as an outline
ol = Outline(Input=nek)
dol = Show(ol, view)
dol.DiffuseColor = dol.AmbientColor = [0.2, 0.2, 0.2]
dol.LineWidth = 1.5

# midplane z = 0.5 coloured by phi_p
sl = Slice(Input=nek)
sl.SliceType.Origin = [0.5, 0.5, 0.5]
sl.SliceType.Normal = [0, 0, 1]
dsl = Show(sl, view)
dsl.Ambient = 1.0  # flat colours on the slice, no shading
dsl.Diffuse = 0.0
# edge of the (white) slice so it reads as a plane
sle = Outline(Input=sl)
dsle = Show(sle, view)
dsle.DiffuseColor = dsle.AmbientColor = [0.6, 0.6, 0.6]
ColorBy(dsl, ("POINTS", "s1"))
lut = GetColorTransferFunction("s1")
# white where there are no particles, blue -> dark red in the cloud
lut.RGBPoints = [0.0, 1.0, 1.0, 1.0, 0.006, 0.55, 0.7, 0.95, 0.015, 0.1, 0.25, 0.75, 0.025, 0.55, 0.0, 0.1]
lut.ColorSpace = "Lab"
lut.RescaleTransferFunction(0.0, 0.025)
dsl.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view)
bar.Title = "particle volume fraction phi_p"
bar.ComponentTitle = ""
bar.TitleColor = bar.LabelColor = [0, 0, 0]
bar.TitleFontSize = 16
bar.LabelFontSize = 14

# EL particles (true diameter) and the IBM sphere (true radius)
g = Glyph(Input=el, GlyphType="Sphere")
g.GlyphType.ThetaResolution = g.GlyphType.PhiResolution = 12
g.ScaleArray = ["POINTS", "diameter"]
g.ScaleFactor = 1.0
g.GlyphMode = "All Points"
dg = Show(g, view)
dg.ColorArrayName = ["POINTS", ""]
dg.DiffuseColor = [0.85, 0.25, 0.1]

gb = Glyph(Input=body, GlyphType="Sphere")
gb.GlyphType.ThetaResolution = gb.GlyphType.PhiResolution = 48
gb.ScaleArray = ["POINTS", "radius"]
gb.ScaleFactor = 2.0
gb.GlyphMode = "All Points"
dgb = Show(gb, view)
dgb.DiffuseColor = [0.75, 0.75, 0.8]
dgb.Opacity = 0.85

# time label
ann = AnnotateTimeFilter(Input=el)
ann.Format = "t = {time:.2f}"
dann = Show(ann, view)
dann.Color = [0, 0, 0]
dann.FontSize = 18
dann.WindowLocation = "Upper Left Corner"

view.CameraPosition = [1.9, 1.6, 2.6]
view.CameraFocalPoint = [0.5, 0.5, 0.5]
view.CameraViewUp = [0, 1, 0]
view.CameraViewAngle = 30

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
times = list(el.TimestepValues)
for i, t in enumerate(times):
    scene.AnimationTime = t
    Render(view)
    SaveScreenshot(os.path.join(frames, f"frame_{i:04d}.png"), view, ImageResolution=view.ViewSize)
print(f"{len(times)} frames in {frames}")

cmd = ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(a.fps), "-i", os.path.join(frames, "frame_%04d.png"),
       "-c:v", "libx264", "-pix_fmt", "yuv420p", "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", out]
subprocess.run(cmd, check=True)
print("wrote", out)
