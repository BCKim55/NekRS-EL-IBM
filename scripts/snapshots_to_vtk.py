#!/usr/bin/env python3
"""Convert a case's snapshots/ (el_XXXXX.csv, body.csv) to ParaView time series.

    python3 scripts/snapshots_to_vtk.py cases/el_ibm_collision [--dp 0.01]

Writes snapshots/vtk/el_XXXXX.vtp + snapshots/el.pvd (EL particles: velocity,
contact, id, diameter) and snapshots/vtk/body_XXXXX.vtp + snapshots/body.pvd
(IBM body: velocity, radius). Plain XML polydata, no VTK library needed.
In ParaView open <case>.nek5000 (fluid), el.pvd and body.pvd, and show the
two point sets with Glyph (sphere, scale by diameter / 2 radius).
"""
import argparse
import csv
import glob
import os
import re


def vtp(path, points, scalars, vectors):
    n = len(points)
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n'
                f' <PolyData>\n  <Piece NumberOfPoints="{n}" NumberOfVerts="{n}">\n   <Points>\n'
                '    <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        f.write("".join(f"     {x} {y} {z}\n" for x, y, z in points))
        f.write("    </DataArray>\n   </Points>\n   <PointData>\n")
        for name, vals in scalars.items():
            f.write(f'    <DataArray type="Float64" Name="{name}" format="ascii">\n'
                    + "".join(f"     {v}\n" for v in vals) + "    </DataArray>\n")
        for name, vals in vectors.items():
            f.write(f'    <DataArray type="Float64" Name="{name}" NumberOfComponents="3" format="ascii">\n'
                    + "".join(f"     {u} {v} {w}\n" for u, v, w in vals) + "    </DataArray>\n")
        f.write('   </PointData>\n   <Verts>\n    <DataArray type="Int64" Name="connectivity" format="ascii">\n'
                + "".join(f"     {i}\n" for i in range(n))
                + '    </DataArray>\n    <DataArray type="Int64" Name="offsets" format="ascii">\n'
                + "".join(f"     {i + 1}\n" for i in range(n))
                + "    </DataArray>\n   </Verts>\n  </Piece>\n </PolyData>\n</VTKFile>\n")


def pvd(path, entries):
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n <Collection>\n')
        for t, file in entries:
            f.write(f'  <DataSet timestep="{t}" file="{file}"/>\n')
        f.write(" </Collection>\n</VTKFile>\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--dp", type=float, default=0.01, help="EL particle diameter when the csv has no dp column")
    a = ap.parse_args()
    snap = os.path.join(a.case, "snapshots")
    out = os.path.join(snap, "vtk")
    os.makedirs(out, exist_ok=True)

    body = {int(r["step"]): r for r in csv.DictReader(open(os.path.join(snap, "body.csv")))}
    el_entries, body_entries = [], []
    for path in sorted(glob.glob(os.path.join(snap, "el_*.csv"))):
        step = int(re.search(r"el_(\d+)\.csv", path).group(1))
        rows = list(csv.DictReader(open(path)))
        t = float(body[step]["t"]) if step in body else step
        name = f"vtk/el_{step:05d}.vtp"
        vtp(os.path.join(snap, name), [(r["x"], r["y"], r["z"]) for r in rows],
            {"id": [r["id"] for r in rows], "contact": [r["contact"] for r in rows],
             "diameter": [r.get("dp", a.dp) for r in rows], "attached": [r.get("attached", 0) for r in rows]},
            {"velocity": [(r["u"], r["v"], r["w"]) for r in rows]})
        el_entries.append((t, name))
        if step in body:
            b = body[step]
            name = f"vtk/body_{step:05d}.vtp"
            vtp(os.path.join(snap, name), [(b["x"], b["y"], b["z"])], {"radius": [b["radius"]]},
                {"velocity": [(b["u"], b["v"], b["w"])]})
            body_entries.append((t, name))
    pvd(os.path.join(snap, "el.pvd"), el_entries)
    pvd(os.path.join(snap, "body.pvd"), body_entries)
    print(f"wrote {len(el_entries)} EL and {len(body_entries)} body snapshots: {snap}/el.pvd, {snap}/body.pvd")


if __name__ == "__main__":
    main()
