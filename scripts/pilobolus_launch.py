#!/usr/bin/env python3
"""Post-process cases/pilobolus_launch logs ("launch:" lines).

    python3 scripts/pilobolus_launch.py --prescribed build/pil_prescribed.log
        [--free build/pil_free.log] [--sap build/pil_sap.log] [--out docs/img/pilobolus_launch.png]

Prints, for the prescribed run, the fluid force on the sporangium and the
deceleration it alone would give (a_f = F/(rho_b V)) during the coasting
phase (t > t_rise), against the measured ~1300 m/s2 (Yafetto et al. 2008),
and the drag coefficient against Schiller-Naumann at the instantaneous Re;
for the free-flight run the actual deceleration of the dry sphere; for the
sap run the stuck mass and the sporangium's force with the droplets on it.
Figure: v(t), a(t), F(t), C_d(t) of the runs that were given.
"""
import argparse
import re

import numpy as np

RHO_F, MU, RHO_B, D = 1.2, 1.8e-5, 1150.0, 3.5e-4
T_RISE, DECEL_MEAS = 2.2e-4, 1300.0


def parse(path):
    rows, free = [], []
    pat = re.compile(r"^launch: step (\d+)\s+t ([\d.eE+-]+)\s+x ([\d.eE+-]+)\s+v ([\d.eE+-]+)\s+a ([\d.eE+-]+)\s+Fm ([\d.eE+-]+) ([\d.eE+-]+) ([\d.eE+-]+)"
                     r"\s+F ([\d.eE+-]+)\s+Cd ([\d.eE+-]+)\s+a_f ([\d.eE+-]+)\s+slip ([\d.eE+-]+)(.*)$")
    extra = re.compile(r"drops (\d+)\s+stuck (\d+)\s+m_stuck ([\d.eE+-]+)\s+P_sap ([\d.eE+-]+)")
    for l in open(path, errors="ignore"):
        m = pat.match(l)
        if not m:
            continue
        r = [float(v) for v in m.groups()[:12]]
        tail = m.group(13)
        e = extra.search(tail)
        r += [float(v) for v in e.groups()] if e else [np.nan] * 4
        rows.append(r)
        free.append("free" in tail)
    return np.array(rows), np.array(free)


def schiller_naumann(v):
    re_ = RHO_F * np.abs(v) * D / MU
    return 24 / np.maximum(re_, 1e-9) * (1 + 0.15 * re_**0.687)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prescribed", required=True)
    ap.add_argument("--free")
    ap.add_argument("--sap")
    ap.add_argument("--out", default="pilobolus_launch.png")
    a = ap.parse_args()
    V = np.pi / 6 * D**3
    runs = {}
    for name, path in (("prescribed", a.prescribed), ("free", a.free), ("sap", a.sap)):
        if path:
            runs[name] = parse(path)

    P, _ = runs["prescribed"]
    t, v, af, cd = P[:, 1], P[:, 3], P[:, 10], P[:, 9]
    coast = t > T_RISE + 2e-5
    if coast.any():
        print(f"prescribed, coasting (t > t_rise): fluid-only deceleration a_f = {-af[coast].mean():.0f} m/s2 "
              f"(range {-af[coast].max():.0f} to {-af[coast].min():.0f}), measured {DECEL_MEAS:.0f}: "
              f"the sap must supply {DECEL_MEAS + af[coast].mean():.0f} m/s2 ({(DECEL_MEAS + af[coast].mean()) / DECEL_MEAS:.0%}); "
              f"C_d {cd[coast].mean():.3f} vs Schiller-Naumann {schiller_naumann(v[coast]).mean():.3f} at Re {RHO_F * v[coast].mean() * D / MU:.0f}")
    if "free" in runs:
        F, fr = runs["free"]
        tf, vf, acc = F[:, 1], F[:, 3], F[:, 4]
        m = fr & (tf > T_RISE + 5e-5)
        if m.any():
            print(f"free flight: v {vf[fr][0]:.3f} -> {vf[-1]:.3f} m/s over {tf[-1] - tf[fr][0]:.2e} s, "
                  f"mean deceleration {-(vf[-1] - vf[fr][0]) / (tf[-1] - tf[fr][0]):.0f} m/s2 (measured {DECEL_MEAS:.0f})")
    if "sap" in runs:
        S, _ = runs["sap"]
        print(f"sap: {int(S[-1, 12])} droplets in the domain at the end, {int(S[-1, 13])} stuck "
              f"(mass {S[-1, 14]:.3e} kg = {S[-1, 14] / (RHO_B * V):.1%} of the sporangium), "
              f"force on the sporangium at the end {S[-1, 8]:.3e} N (single phase {P[-1, 8]:.3e} N)")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5))
    for name, (R, fr) in runs.items():
        tt = R[:, 1] * 1e3
        axs[0, 0].plot(tt, R[:, 3], label=name)
        axs[0, 1].plot(tt, R[:, 4], label=name)
        axs[1, 0].plot(tt, R[:, 8] * 1e6, label=f"{name}: F on body")
        if name == "prescribed":
            axs[1, 0].plot(tt, R[:, 5] * 1e6, "--", lw=0.8, label="prescribed: marker sum F_m")
        axs[1, 1].plot(tt, R[:, 9], label=name)
    axs[0, 1].axhline(-DECEL_MEAS, color="C3", ls=":", label="measured coasting deceleration")
    axs[1, 1].plot(t * 1e3, schiller_naumann(v), "k--", lw=0.8, label="Schiller-Naumann at Re(t)")
    axs[1, 1].set_ylim(0, 5)
    for ax, yl in zip(axs.flat, ("v (m/s)", "a (m/s$^2$)", "x force (µN)", "$C_d$")):
        ax.set_xlabel("t (ms)")
        ax.set_ylabel(yl)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    axs[0, 0].set_title("sporangium velocity")
    axs[0, 1].set_title("acceleration (prescribed) / measured (free)")
    axs[1, 0].set_title("fluid force on the sporangium")
    axs[1, 1].set_title("drag coefficient")
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
