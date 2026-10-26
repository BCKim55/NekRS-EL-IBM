#!/usr/bin/env python3
"""Figure for docs/design.md section 5 (validation against the originals).

    python3 scripts/plot_validation.py <nek5000 twoway lines> <nekRS twoway lines> <sphere lines> --out docs/img/validation.png

Panels: (a) ppiclF stokes_2d regression, pristine vs patched, max error per
rank count (numbers from scripts/ppiclf_regression.py); (b) settling-box
budget, Nek5000 + original ppiclF vs nekRS + ppiclF; (c) fixed-sphere drag
history at Re 10 against the unbounded correlations and the periodic-array
(Sangani-Acrivos) Stokes bound for the 4d box.
"""
import argparse
import re

import numpy as np

# scripts/ppiclf_regression.py output: max |v - exact| for np 1, 2, 3
REG = {"original ppiclF": [3.74e-8, 1.84e-3, 1.84e-3], "patched ppiclF": [3.74e-8, 3.74e-8, 3.74e-8]}


def budget(path):
    rows = []
    for l in open(path, errors="ignore"):
        m = re.search(r"t\s+([\d.eE+-]+)\s+Pf\s+([\d.eE+-]+)\s+Pp\s+([\d.eE+-]+)\s+sumf\s+([\d.eE+-]+)", l)
        if m:
            rows.append([float(v) for v in m.groups()])
    return np.array(rows)


def sphere(path):
    rows = []
    for l in open(path, errors="ignore"):
        m = re.search(r"t\s+([\d.eE+-]+).*Cd\s+([\d.eE+-]+)", l)
        if m:
            rows.append([float(v) for v in m.groups()])
    return np.array(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("nek5000")
    ap.add_argument("nekrs")
    ap.add_argument("sphere")
    ap.add_argument("--out", default="validation.png")
    a = ap.parse_args()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))

    ax = axs[0]
    x = np.arange(3)
    for i, (lab, e) in enumerate(REG.items()):
        ax.bar(x + (i - 0.5) * 0.35, e, 0.35, label=lab)
    ax.set_yscale("log")
    ax.set_xticks(x, ["1 rank", "2 ranks", "3 ranks"])
    ax.set_ylabel("max |v - exact| over all particles")
    ax.set_title("(a) ppiclF stokes_2d: library before / after patches")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    ax = axs[1]
    A, B = budget(a.nek5000), budget(a.nekrs)
    for col, lab, ls in ((1, "fluid momentum $P_f$", "-"), (2, "particle momentum $P_p$", "--"), (3, r"total drag $\sum f$", ":")):
        ax.plot(A[:, 0], A[:, col], "k" + ls, label=f"Nek5000 + original ppiclF: {lab}")
        ax.plot(B[:, 0], B[:, col], "C3" + ls, lw=1.2, alpha=0.8, label=f"nekRS + ppiclF (this repo): {lab}")
    ax.set_xlabel("t")
    ax.set_ylabel("momentum / drag")
    ax.set_title("(b) settling box, two-way (4 ranks)")
    ax.legend(fontsize=6.5)
    ax.grid(alpha=0.3)

    ax = axs[2]
    S = sphere(a.sphere)
    ax.plot(S[:, 0], S[:, 1], "C0-", label="this repo's IBM, NekIBM simple set-up (4d box)")
    ax.axhline(4.15, color="k", ls="--", lw=0.8, label="Schiller-Naumann, unbounded (4.15)")
    ax.axhline(4.30, color="k", ls=":", lw=0.8, label="standard drag curve, unbounded (4.30)")
    ax.axhline(4.15 * 1.53, color="gray", ls="-.", lw=0.8, label="Stokes-limit periodic-array bound, 4d spacing (x1.53)")
    ax.set_ylim(3, 8)
    ax.set_xlabel("t")
    ax.set_ylabel("$C_d$")
    ax.set_title("(c) fixed sphere, Re = 10")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
