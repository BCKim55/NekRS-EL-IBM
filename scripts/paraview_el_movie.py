#!/usr/bin/env pvpython
"""Movie of an EL-IBM case: domain outline, a midplane slice coloured by the
particle volume fraction phi_p (checkpoint scalar00 = reader array s1), the
EL parcels (true diameter) and the IBM body as sphere glyphs.

    python3 scripts/snapshots_to_vtk.py <case>
    /Applications/ParaView-5.10.1.app/Contents/bin/pvpython scripts/paraview_el_movie.py <case> \\
        [--out <case>/snapshots/movie.mp4] [--fps 12] [--phimax 0.025] [--normal z]
        [--time-unit 1 --time-label t] [--view xy|iso] [--zoom 1]

Defaults reproduce the el_ibm_collision movie (iso view, t in code units).
cases/pilobolus_launch/sap: --phimax 0.1 --time-unit 1e3 --time-label "t (ms)"
--view xy --zoom 3 (the flight path fills the frame, the far field is cut).
The slice updates whenever a checkpoint exists; the parcels every snapshot.
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
ap.add_argument("--phimax", type=float, default=0.025)
ap.add_argument("--normal", default="z", choices="xyz")
ap.add_argument("--time-unit", type=float, default=1.0, help="label = time * unit")
ap.add_argument("--time-label", default="t")
ap.add_argument("--view", default="iso", choices=["iso", "xy"])
ap.add_argument("--zoom", type=float, default=1.0, help="> 1 moves the camera closer to the domain centre")
a = ap.parse_args()
case = os.path.abspath(a.case)
snap = os.path.join(case, "snapshots")
out = a.out or os.path.join(snap, "movie.mp4")
frames = os.path.join(snap, "frames")
os.makedirs(frames, exist_ok=True)
for f in glob.glob(os.path.join(frames, "frame_*.png")):
    os.remove(f)

nek5000 = glob.glob(os.path.join(case, "*.nek5000"))
assert len(nek5000) == 1, nek5000
nek = OpenDataFile(nek5000[0])
nek.PointArrays = ["velocity", "s1"]
el = OpenDataFile(os.path.join(snap, "el.pvd"))
body = OpenDataFile(os.path.join(snap, "body.pvd"))

nek.UpdatePipeline()
lo = nek.GetDataInformation().GetBounds()[0::2]
hi = nek.GetDataInformation().GetBounds()[1::2]
ctr = [0.5 * (l + h) for l, h in zip(lo, hi)]
size = max(h - l for l, h in zip(lo, hi))

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1280, 960] if a.view == "iso" else [1600, 700]
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

# midplane slice coloured by phi_p
k = "xyz".index(a.normal)
sl = Slice(Input=nek)
sl.SliceType.Origin = ctr
sl.SliceType.Normal = [float(i == k) for i in range(3)]
dsl = Show(sl, view)
dsl.Ambient = 1.0  # flat colours on the slice, no shading
dsl.Diffuse = 0.0
sle = Outline(Input=sl)
dsle = Show(sle, view)
dsle.DiffuseColor = dsle.AmbientColor = [0.6, 0.6, 0.6]
ColorBy(dsl, ("POINTS", "s1"))
lut = GetColorTransferFunction("s1")
# white where there are no particles, blue -> dark red in the cloud
pm = a.phimax
lut.RGBPoints = [0.0, 1.0, 1.0, 1.0, 0.24 * pm, 0.55, 0.7, 0.95, 0.6 * pm, 0.1, 0.25, 0.75, pm, 0.55, 0.0, 0.1]
lut.ColorSpace = "Lab"
lut.RescaleTransferFunction(0.0, pm)
dsl.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view)
bar.Title = "particle volume fraction phi_p"
bar.ComponentTitle = ""
bar.TitleColor = bar.LabelColor = [0, 0, 0]
bar.TitleFontSize = 16
bar.LabelFontSize = 14

# EL parcels (true diameter) and the IBM body (true radius)
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
ann.Scale = a.time_unit
ann.Format = a.time_label + " = {time:.2f}"
dann = Show(ann, view)
dann.Color = [0, 0, 0]
dann.FontSize = 18
dann.WindowLocation = "Upper Left Corner"

Render(view)  # the first render resets the camera; set ours after it
if a.view == "iso":
    view.CameraPosition = [c + f * size / a.zoom for c, f in zip(ctr, (1.4, 1.1, 2.1))]
    view.CameraViewUp = [0, 1, 0]
else:  # xy: look down the slice normal, slightly tilted so the glyphs read as spheres
    up = [0, 1, 0] if k != 1 else [0, 0, 1]
    pos = list(ctr)
    pos[k] += 2.2 * size / a.zoom
    pos[(k + 1) % 3] += 0.15 * size / a.zoom
    view.CameraPosition = pos
    view.CameraViewUp = up
view.CameraFocalPoint = ctr
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
