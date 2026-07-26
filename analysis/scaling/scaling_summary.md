# Scaling analysis from erfc τ₅₀

Using τ₅₀ from `erfc(ln τ)` fits. Angular speed ω = 2π/T_rot. Optical q = (4π/λ) sin(θ/2) with λ = 650 nm.

Notes predict τ_c = Δθ_c(q)/ω, so 1/τ_c vs ω is linear through the origin with slope 1/Δθ_c. Ideal centered rotation further gives τ⁻¹ ∝ Ω sinθ.

## Plots

- `inv_tau50_vs_omega.png` — 1/τ₅₀ vs ω by angle (solid = free intercept, dashed = origin)
- `tau50_vs_inv_omega.png` — τ₅₀ vs 1/ω
- `slopes_vs_q_and_sintheta.png` — slopes vs q and vs sinθ
- `omega_tau50_vs_omega.png` — ω·τ₅₀ vs ω (should be flat at fixed angle)

## Per-angle slopes of 1/τ₅₀ = m·ω (+ b)

| Angle | q (µm⁻¹) | sinθ | m (origin) | m (free) | intercept b (s⁻¹) | Δθ_c (deg) |
|---|---|---|---|---|---|---|
| small angle (nom. 0°) | 0.000 | 0.0000 | 195.343±20.685 | 160.629±1.949 | 15.4±0.7 | 0.293±0.031 |
| 10° | 1.685 | 0.1736 | 602.292±12.654 | 572.559±2.955 | 11.9±0.9 | 0.095±0.002 |
| 15° | 2.523 | 0.2588 | 931.226±9.185 | 917.827±3.789 | 7.7±1.7 | 0.062±0.001 |
| 20° | 3.357 | 0.3420 | 1244.903±12.649 | 1266.580±3.757 | -9.7±1.2 | 0.046±0.000 |

Δθ_c = 1/m_origin (characteristic diffuser rotation angle to decorrelate).

Numeric slopes: `angle_slopes.csv`
