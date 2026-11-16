# Representing a free liquid with pairwise-interacting Lagrangian parcels

Literature survey (2026-11-16) behind `cohesion_model` 1 and 2 in
`src/ppiclf_user_elibm.f` (design doc 7.14–7.15). Target: water,
σ = 0.072 N/m, ρ = 1000, in a DEM neighbour-list code (ppiclF).

## 1. Pairwise-force SPH (PF-SPH)

* Tartakovsky & Meakin 2005, Phys. Rev. E 72, 026301. Pair force added to
  the SPH momentum equation, `F_ij = s_ij cos(1.5 π r/(κh)) e_ij` for
  r ≤ κh: repulsive for r < κh/3, attractive beyond. σ is measured (virial
  pressure of a drop vs 2σ/R), s_fs/s_ff sets the contact angle.
* Tartakovsky & Panchenko 2016, J. Comput. Phys. 305, 1119,
  doi:10.1016/j.jcp.2015.08.037. `F_ij = −s φ(r) r/|r|`,
  `φ(r) = r[−A exp(−r²/2r_a²) + exp(−r²/2r_b²)]`, r_a = r_b/2 = Δ/2.
  Analytic 3D relation `σ = (π/8) n² s ∫ g(r) φ(r) r⁴ dr` (single liquid,
  reduction of their Eq. 34); 2D version in the paper.
* Whebell, Moroney, Turner, Pethiyagoda, McCue 2024/25 (arXiv 2412.03810,
  J. Fluid Mech.). Resolution-independent form `F_i = (H/m_i) Σ s_ij
  f(|x_ij|/H) x_ij/|x_ij|`, `f(τ) = 1 + ατ² − 2(2+α)τ³ + (3+α)τ⁴`,
  α_ff = −23.5 (zero crossing at 1.1 Δx, the FCC rest spacing), α_fs = −11,
  s in N/m. Laplace calibration: **σ = 30.96 s_ff**, independent of Δx from
  20 to 80 µm. Water: ρ 1000, µ 0.89 mPa s, c₀ 60 m/s, H 120 µm.
  → this is the profile of `cohesion_model = 2`.
* Akinci, Akinci, Teschner 2013, ACM TOG 32(6),
  doi:10.1145/2508363.2508395. Cohesion kernel `C(r) = 32/(πh⁹)
  {(h−r)³r³, 2r > h; 2(h−r)³r³ − h⁶/64, 2r ≤ h}` (repulsive below h/2),
  plus a curvature term from per-particle normals `n_i = h Σ (m_j/ρ_j) ∇W`,
  `F = −γ m_i (n_i − n_j)`, and the neighbour-deficiency factor
  `K_ij = 2ρ₀/(ρ_i + ρ_j)`. Cohesion-only models cluster into "cobwebs";
  the n_i − n_j term fixes that. γ is tuned, no σ relation.
* Becker & Teschner 2007 (SCA): kernel-weighted cohesion only, clusters.
* Nugent & Posch 2000, Phys. Rev. E 62, 4968: van der Waals cohesive
  pressure as a long-range pair force with twice the repulsive smoothing
  length; σ emergent, checked by Laplace and drop oscillation.

## 2. Double density relaxation

Clavet, Beaudoin, Poulin 2005 (SCA, doi:10.2312/SCA/SCA05/219-228).
Position-based: `ρ_i = Σ (1 − r/h)²`, `P = k(ρ − ρ₀)`,
`ρ_near = Σ (1 − r/h)³`, `P_near = k_near ρ_near`; displacement
`D_ij = dt² [P (1 − r/h) + P_near (1 − r/h)²] r̂`. The far term pulls when
ρ < ρ₀ (free surface attracted inward), the near term stops collapse: a
rest spacing without a lattice. No σ calibration; implicit, not a force
law — not usable in a velocity-Verlet DEM loop.

## 3. MPS potential force

Kondo, Koshizuka, Takimoto 2007 (Trans. JSCES, doi:10.11421/jsces.2007.20070021).
`F(r) = C (r − r_min)(r − r_s)`: repulsive below the initial spacing
r_min, attractive to r_s. C from the work to separate two half-spaces
= 2σ per area; ISIJ Int. 56(7) 2016 writes `E = 2σd²/(h³ Σφ)`. Simple
closed form, but the half-space derivation assumes a 3D bulk and
r_s ≈ 1.5–2 r_min.

## 4. DEM representations of liquid

* Capillary-bridge DEM: Willett et al. 2000 (Langmuir 16, 9396)
  `F = 2πRγ cosθ/(1 + 2.1S⁺ + 10S⁺²)`, `S⁺ = a/(2√(V/R))`, rupture
  `a_c = (1 + θ/2)(V*^{1/3} + V*^{2/3}/10)`, **valid for θ < 50° and
  V* = V/R³ < 0.1**; Rabinovich et al. 2005 (Langmuir 21, 10992); Lian,
  Thornton, Adams 1993 (JCIS 161, 138). All describe a small pendular
  bridge pinned on two *solid* spheres: force bounded by 2πRγ, scales with
  the solid radius, needs a solid contact angle, and per-contact adhesion
  does not sum to an interface tension. When the parcels are the liquid
  (V* ~ 1, no solid) the model is structurally wrong — the reason the
  bridge runs of 7.11–7.12 are a granular chain, not a sheet.
* Many-body DPD: Warren 2003, Phys. Rev. E 68, 066702
  `F = A w_c e + B(ρ̄_i + ρ̄_j) w_d e`, A < 0, B > 0; EOS
  `p = ρkT + (π/30)(A + 2Br_d⁴ρ)ρ²`; σ measured (pendant drop vs
  Young–Laplace). Arienti et al. 2011, J. Chem. Phys. 134, 204114: MDPD
  threads follow Rayleigh–Plateau and the pinch-off power laws. Closest
  published analogue of "DEM with a density-dependent pair force".
* No DEM paper found that models free liquid as cohesive parcels calibrated
  to σ (two searches); treat as a gap.

## 5. Wet-collision and wall-impact rules in E–L spray codes

* O'Rourke 1981: collision probability Poisson, `n̄ = n₂π(r₁+r₂)² v_rel dt/V`,
  offset `b = (r₁+r₂)√Y`, `b_crit = (r₁+r₂)√min(1, 2.4 f(γ)/We)`,
  `f(γ) = γ³ − 2.4γ² + 2.7γ`; b < b_crit coalesce, else grazing bounce.
* Ashgriz & Poo 1990, J. Fluid Mech. 221, 183: We = ρu²d_s/σ,
  B = 2b/(d_s + d_l); head-on water drops coalesce below We ≈ 15–20 and
  separate reflexively above; stretching separation
  `B = C_a We^{−1/2}[1 + (µ/σ)(ρd_s/σ)^{1/2}][Δ⁻³ − 2.4Δ⁻² + 2.7Δ⁻¹]^{1/2}`,
  C_a ≈ 1.14 (Jiang 1992); bouncing only at high B (Estrade 1999).
* Wall impact: Mundo, Sommerfeld, Tropea 1995 `K = We^{0.5} Re^{0.25}`,
  deposit for K < 57.7; Bai & Gosman 1995 wetted wall rebound→spread We = 5,
  spread→splash `We = 1320 La^{−0.18}`, La = ρσD/µ²; Stanton & Rutland 1996
  `E² = (ρV_n²D/σ)/(min(h₀/D,1) + δ/D)`, stick for E < 16, splash E > 57.7.
  → parcels hitting the cap at 15 m/s (We 80) spread, do not bounce:
  e → 0 and tangential damping (`tangential_damping`).

## 6. Calibration tests (water, our scales)

| test | target | our numbers |
|---|---|---|
| Laplace | ΔP = 2σ/R (sphere), σ/R (cylinder) | 75 µm jet: 1.92 kPa |
| drop oscillation l = 2 | ω = √(8σ/ρR³), damping 5µ/ρR² | R 37.5 µm: period ≈ 60 µs |
| Rayleigh–Plateau | kR = 0.697, ω_max = 0.343 √(σ/ρR³) | R 37.5 µm: growth time ≈ 79 µs |
| Taylor–Culick rim | u = √(2σ/ρh) | h 1, 3, 7 µm: 12.0, 6.9, 4.5 m/s |
| thread tension | πdσ − (2σ/d)(πd²/4) = πdσ/2 | d 75 µm: 8.5 µN; d 25 µm: 2.8 µN |

## 7. Assessment for the Pilobolus sap

The resolution is the hard constraint: every model above wants several
particles across a feature, and the Video S3 sheet (1–7 µm) is 4–25×
thinner than one 25 µm parcel. The sheet is therefore a parcel monolayer
whose mass per area is set by the in-plane spacing (h_eff = V_p/Δ²), and
the 3D σ relations do not apply to it: calibrate the monolayer edge
tension to 2σ instead (then Taylor–Culick follows with ρh_eff), and
accept that a 75 µm jet is 3 parcels across. Recommendation adopted:
PF-SPH bond with the Whebell quartic (`cohesion_model = 2`), spring as
the repulsive core, dashpot as the only dissipation, e ≈ 0 inside the
attractive range (O'Rourke: water coalesces for We < 15), s set by the
static lattice sums of `scripts/pf_cohesion_calibration.py`; the Akinci
K_ij neighbour-deficiency factor is the next refinement if the thread
(single file, needs 0.16 σ) and the sheet (0.043 σ) cannot share one s.
