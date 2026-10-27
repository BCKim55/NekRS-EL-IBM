#!/usr/bin/env pvpython
"""ParaView set-up for cases/el_ibm_collision: fluid midplane slice (nekRS
checkpoints), EL particles and the IBM sphere as sphere glyphs.

    python3 scripts/snapshots_to_vtk.py cases/el_ibm_collision
    /Applications/ParaView-5.10.1.app/Contents/bin/pvpython scripts/paraview_collision.py cases/el_ibm_collision

Writes <case>/collision.pvsm (File > Load State in the ParaView GUI, then play
the animation; the reader shows the phi_p checkpoint field as "s1") and <case>/snapshots/paraview.png. The nekRS checkpoints are
every 40 steps, the particle snapshots every 4, so the particles move more
smoothly than the slice.
"""
import os
import sys

from paraview.simple import *  # noqa: F401,F403

case = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "cases/el_ibm_collision")
nek = OpenDataFile(os.path.join(case, "collision.nek5000"))
# the VisIt Nek5000 reader loads nothing by default; s1 = scalar00 = projected particle volume fraction phi_p
nek.PointArrays = ["pressure", "velocity", "s1"]
RenameSource("fluid (nekRS checkpoints)", nek)
el = OpenDataFile(os.path.join(case, "snapshots", "el.pvd"))
RenameSource("EL particles", el)
body = OpenDataFile(os.path.join(case, "snapshots", "body.pvd"))
RenameSource("IBM body", body)

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1000, 800]
sl = Slice(Input=nek)
sl.SliceType.Origin = [0.5, 0.5, 0.5]
sl.SliceType.Normal = [0, 0, 1]
RenameSource("midplane z = 0.5", sl)
d = Show(sl, view)
ColorBy(d, ("POINTS", "velocity", "Magnitude"))
d.SetScalarBarVisibility(view, True)
# second slice coloured by phi_p (hidden; toggle the eye icon in the GUI)
sl2 = Slice(Input=nek)
sl2.SliceType.Origin = [0.5, 0.5, 0.5]
sl2.SliceType.Normal = [0, 0, 1]
RenameSource("midplane: particle volume fraction phi_p (s1)", sl2)
d2 = Show(sl2, view)
ColorBy(d2, ("POINTS", "s1"))
Hide(sl2, view)

g = Glyph(Input=el, GlyphType="Sphere")
g.GlyphType.ThetaResolution = g.GlyphType.PhiResolution = 12
g.ScaleArray = ["POINTS", "diameter"]
g.ScaleFactor = 1.0
g.GlyphMode = "All Points"
RenameSource("EL spheres (d_p)", g)
de = Show(g, view)
ColorBy(de, ("POINTS", "contact"))

gb = Glyph(Input=body, GlyphType="Sphere")
gb.GlyphType.ThetaResolution = gb.GlyphType.PhiResolution = 48
gb.ScaleArray = ["POINTS", "radius"]
gb.ScaleFactor = 2.0  # Sphere glyph radius 0.5 x scale -> diameter = 2 radius
gb.GlyphMode = "All Points"
RenameSource("IBM sphere", gb)
db = Show(gb, view)
db.DiffuseColor = [0.85, 0.85, 0.85]

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
scene.AnimationTime = 0.6
view.ResetCamera()
Render()
SaveScreenshot(os.path.join(case, "snapshots", "paraview.png"), view)
SaveState(os.path.join(case, "collision.pvsm"))
print("wrote", os.path.join(case, "collision.pvsm"), "and snapshots/paraview.png")
