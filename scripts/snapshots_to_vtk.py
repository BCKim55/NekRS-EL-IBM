#!/usr/bin/env python3
"""Convert a case's snapshots/ (el_XXXXX.csv, body.csv) to ParaView time series.

    python3 scripts/snapshots_to_vtk.py cases/el_ibm_collision [--dp 0.01]

Writes snapshots/vtk/el_XXXXX.vtp + snapshots/el.pvd (EL particles: velocity,
contact, id, diameter) and snapshots/vtk/body_XXXXX.vtp + snapshots/body.pvd
(IBM body: velocity, radius). Plain XML polydata, no VTK library needed.
In ParaView open <case>.nek5000 (fluid), el.pvd and body.pvd, and show the
two point sets with Glyph (sphere, scale by diameter / 2 radius).
The EL parcels are point particles and stay points. --body-surface also
writes body_surf.pvd with the body as a triangulated sphere (or a cylinder
along z of --cylinder-length), which displays as a solid with no filter.
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


def unit_sphere(nth=8, nph=12):
    """UV sphere: points on the unit sphere and triangles (index triples)."""
    import math
    pts = [(0.0, 0.0, 1.0)]
    for i in range(1, nth):
        th = math.pi * i / nth
        for j in range(nph):
            ph = 2 * math.pi * j / nph
            pts.append((math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)))
    pts.append((0.0, 0.0, -1.0))
    last = len(pts) - 1
    tri = []
    for j in range(nph):
        tri.append((0, 1 + j, 1 + (j + 1) % nph))
        tri.append((last, 1 + (nth - 2) * nph + (j + 1) % nph, 1 + (nth - 2) * nph + j))
    for i in range(nth - 2):
        for j in range(nph):
            a, b = 1 + i * nph + j, 1 + i * nph + (j + 1) % nph
            c, d = a + nph, b + nph
            tri += [(a, c, b), (b, c, d)]
    return pts, tri


def unit_hemisphere(nth=12, nph=32):
    """Solid hemisphere of radius 1: dome z >= 0, flat face at z = 0 (sphere centre at the origin)."""
    import math
    pts, tri = [(0.0, 0.0, 1.0)], []
    for i in range(1, nth + 1):
        th = 0.5 * math.pi * i / nth
        for j in range(nph):
            ph = 2 * math.pi * j / nph
            pts.append((math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)))
    for j in range(nph):
        tri.append((0, 1 + j, 1 + (j + 1) % nph))
    for i in range(nth - 1):
        for j in range(nph):
            a, b = 1 + i * nph + j, 1 + i * nph + (j + 1) % nph
            tri += [(a, a + nph, b), (b, a + nph, b + nph)]
    c = len(pts)
    pts.append((0.0, 0.0, 0.0))
    rim = 1 + (nth - 1) * nph
    for j in range(nph):
        tri.append((c, rim + (j + 1) % nph, rim + j))
    return pts, tri


def oriented(shape, axis):
    """Rotate a unit shape so its z axis points along axis (unit vector)."""
    import math
    pts, tri = shape
    z = axis
    seed = (1.0, 0.0, 0.0) if abs(z[0]) < 0.9 else (0.0, 1.0, 0.0)
    x = (seed[1] * z[2] - seed[2] * z[1], seed[2] * z[0] - seed[0] * z[2], seed[0] * z[1] - seed[1] * z[0])
    n = math.sqrt(sum(v * v for v in x))
    x = tuple(v / n for v in x)
    y = (z[1] * x[2] - z[2] * x[1], z[2] * x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0])
    return [tuple(x[k] * p[0] + y[k] * p[1] + z[k] * p[2] for k in range(3)) for p in pts], tri


def unit_cylinder(nph=48):
    """Cylinder of radius 1 from z = -0.5 to 0.5, with end caps."""
    import math
    pts, tri = [], []
    for z in (-0.5, 0.5):
        for j in range(nph):
            ph = 2 * math.pi * j / nph
            pts.append((math.cos(ph), math.sin(ph), z))
    c0, c1 = len(pts), len(pts) + 1
    pts += [(0.0, 0.0, -0.5), (0.0, 0.0, 0.5)]
    for j in range(nph):
        a, b = j, (j + 1) % nph
        tri += [(a, b, a + nph), (b, b + nph, a + nph), (c0, b, a), (c1, a + nph, b + nph)]
    return pts, tri


def vtp_surfaces(path, centres, sizes, shape, scalars, vectors):
    """Triangulated spheres (shape = unit sphere, sizes = diameters) or cylinders
    (sizes = (diameter, length)) at the centres, point data copied per solid."""
    spts, stri = shape
    pts, tris, sc, vc = [], [], {k: [] for k in scalars}, {k: [] for k in vectors}
    for n, (c, sz) in enumerate(zip(centres, sizes)):
        base = len(pts)
        rx = ry = 0.5 * float(sz[0]) if isinstance(sz, (tuple, list)) else 0.5 * float(sz)
        rz = float(sz[1]) if isinstance(sz, (tuple, list)) else rx
        pts += [(float(c[0]) + rx * x, float(c[1]) + ry * y, float(c[2]) + rz * z) for x, y, z in spts]
        tris += [(base + a, base + b, base + c_) for a, b, c_ in stri]
        for k in scalars:
            sc[k] += [scalars[k][n]] * len(spts)
        for k in vectors:
            vc[k] += [vectors[k][n]] * len(spts)
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n'
                f' <PolyData>\n  <Piece NumberOfPoints="{len(pts)}" NumberOfPolys="{len(tris)}">\n   <Points>\n'
                '    <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        f.write("".join(f"     {x} {y} {z}\n" for x, y, z in pts))
        f.write("    </DataArray>\n   </Points>\n   <PointData>\n")
        for name, vals in sc.items():
            f.write(f'    <DataArray type="Float64" Name="{name}" format="ascii">\n'
                    + "".join(f"     {v}\n" for v in vals) + "    </DataArray>\n")
        for name, vals in vc.items():
            f.write(f'    <DataArray type="Float64" Name="{name}" NumberOfComponents="3" format="ascii">\n'
                    + "".join(f"     {u} {v} {w}\n" for u, v, w in vals) + "    </DataArray>\n")
        f.write('   </PointData>\n   <Polys>\n    <DataArray type="Int64" Name="connectivity" format="ascii">\n'
                + "".join(f"     {a} {b} {c}\n" for a, b, c in tris)
                + '    </DataArray>\n    <DataArray type="Int64" Name="offsets" format="ascii">\n'
                + "".join(f"     {3 * (i + 1)}\n" for i in range(len(tris)))
                + "    </DataArray>\n   </Polys>\n  </Piece>\n </PolyData>\n</VTKFile>\n")


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
    ap.add_argument("--body-surface", action="store_true", help="also write body_surf.pvd (solid body)")
    ap.add_argument("--cylinder-length", type=float, default=0, help="> 0: the body is a cylinder along z")
    a = ap.parse_args()
    cyl = unit_cylinder()
    surf_body = []
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
            if a.body_surface:
                name = f"vtk/body_surf_{step:05d}.vtp"
                d = 2 * float(b["radius"])
                shape = int(float(b.get("shape", 0)))
                if a.cylinder_length > 0:
                    vtp_surfaces(os.path.join(snap, name), [(b["x"], b["y"], b["z"])], [(d, a.cylinder_length)], cyl,
                                 {}, {"velocity": [(b["u"], b["v"], b["w"])]})
                elif shape == 2:
                    # hemisphere: mesh about the sphere centre, 3R/8 behind the centre of mass along the axis
                    ax = (float(b["ax"]), float(b["ay"]), float(b["az"]))
                    cen = tuple(float(b[k]) - 0.375 * 0.5 * d * ax[i] for i, k in enumerate("xyz"))
                    vtp_surfaces(os.path.join(snap, name), [cen], [d], oriented(unit_hemisphere(), ax),
                                 {}, {"velocity": [(b["u"], b["v"], b["w"])]})
                else:
                    vtp_surfaces(os.path.join(snap, name), [(b["x"], b["y"], b["z"])], [d], unit_sphere(16, 32),
                                 {}, {"velocity": [(b["u"], b["v"], b["w"])]})
                surf_body.append((t, name))
    pvd(os.path.join(snap, "el.pvd"), el_entries)
    pvd(os.path.join(snap, "body.pvd"), body_entries)
    print(f"wrote {len(el_entries)} EL and {len(body_entries)} body snapshots: {snap}/el.pvd, {snap}/body.pvd")
    if a.body_surface:
        pvd(os.path.join(snap, "body_surf.pvd"), surf_body)
        print(f"wrote solid body: {snap}/body_surf.pvd")


if __name__ == "__main__":
    main()
