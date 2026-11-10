#!/usr/bin/env python3
"""Window averages of the cap accelerations in a pilobolus sap-jet run.

a_f: fluid force on the cap / m_b, a_c: contact force of the sap (spring-dashpot
+ liquid bridges) / m_b, a_tot = a_f + a_c.  Averages over the time windows of
the bridge-force hand-off (0-0.34 ms injection, 0.4-0.6, 0.6-1.25 ms flight),
against the tracked in-flight deceleration, and a time-history figure.

  python3 scripts/pilobolus_jet_accel_table.py --run label=dir [...] [--png out.png]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from pilobolus_jet_compare import read_log  # noqa: E402

WINDOWS = [(0.0, 0.34e-3), (0.4e-3, 0.6e-3), (0.6e-3, 1.25e-3)]
MEASURED = (-1592.0, 103.0)  # m/s2, Video S3 track 0.4-1.25 ms


def window_mean(t, y, lo, hi):
    m = (t > lo) & (t <= hi) & np.isfinite(y)
    return y[m].mean() if m.any() else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True, help="label=dir (dir holds the *.out log)")
    ap.add_argument("--png", default="")
    args = ap.parse_args()
    runs = [(s.split("=", 1)[0], read_log(s.split("=", 1)[1])) for s in args.run]

    print(f"window averages [m/s2]; measured in-flight {MEASURED[0]:.0f} +- {MEASURED[1]:.0f}")
    hdr = f"{'run':24s} {'window [ms]':12s} {'a_f':>9s} {'a_c':>9s} {'a_f+a_c':>9s} {'a_f/meas':>9s} {'tot/meas':>9s}  {'Fc [uN]':>8s} {'Fbr [uN]':>8s}"
    print(hdr)
    for label, d in runs:
        for lo, hi in WINDOWS:
            af, ac, at = (window_mean(d["t"], d[k], lo, hi) for k in ("a_f", "a_c", "a_tot"))
            fc, fbr = (window_mean(d["t"], d[k], lo, hi) * 1e6 for k in ("Fc", "Fbr"))
            print(f"{label:24s} {lo*1e3:4.2f}-{hi*1e3:<7.2f} {af:9.0f} {ac:9.0f} {at:9.0f} {af/MEASURED[0]:9.2f} {at/MEASURED[0]:9.2f}  {fc:8.1f} {fbr:8.1f}")

    if args.png:
        fig, ax = plt.subplots(3, 1, figsize=(7, 8), sharex=True)
        for label, d in runs:
            t = d["t"] * 1e3
            ax[0].plot(t, d["a_f"], label=label)
            ax[1].plot(t, d["a_c"], label=label)
            ax[2].plot(t, d["a_tot"], label=label)
        for a, name in zip(ax, ("a_f (air)", "a_c (sap contact)", "a_f + a_c")):
            a.axhspan(MEASURED[0] - MEASURED[1], MEASURED[0] + MEASURED[1], color="0.8", label="measured 0.4-1.25 ms" if a is ax[0] else None)
            a.axhline(0, color="k", lw=0.5)
            a.set_ylabel(name + "  [m/s$^2$]")
            for lo, hi in WINDOWS:
                a.axvline(lo * 1e3, color="0.6", lw=0.5, ls=":")
                a.axvline(hi * 1e3, color="0.6", lw=0.5, ls=":")
        ax[1].set_yscale("symlog", linthresh=1e3)
        ax[2].set_yscale("symlog", linthresh=1e3)
        ax[2].set_xlabel("t  [ms]")
        ax[0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(args.png, dpi=150)
        print("wrote", args.png)


if __name__ == "__main__":
    main()
