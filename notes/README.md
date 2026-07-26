# Notes — scratchpad, not results

Working theory notes and predictions for the rotating-diffuser speckle experiment.
**These are ideas and hypotheses, not established truth**, and not a record of what was
measured. They are rough material that may later be edited into papers.

Unlike the raw data (`../test*.txt`, `../test*.png`) and the lab notebook
(`../lab_notebook.md`), these documents are freely editable.

## Contents

### `diffuser_speckle_correlations/main.tex`
*Rotating-Diffuser Speckle as a Spatiotemporal Heating Field* — a tutorial building the
ideal centered-rotation model. Derives the jinc-squared form

    g2(θ,τ) = 1 + β [2 J₁(X)/X]²,   X = (2πD/λ) sinθ sin(Ωτ/2)

with τ₁/ₑ = 1.91499 λ / (π D Ω sinθ), predicts rotation echoes at multiples of the
revolution period, and argues that Ī(Q) plus single-point g₂ do **not** determine the full
sample-plane spectrum S_I(k,ω). Also covers finite detector size, off-axis illumination,
and the connection to thermal transport.

### `rotating_diffuser_extra/main.tex`
*Initial Scaling Analysis of Dynamic Speckle from a Rotating Diffuser* — an analysis
specification aimed directly at the measured curves. Centers on collapse tests:
G(τ) = F(ωτ), then F(Rωτ), then F(Rωτ/ℓ_c(q)); stretched/compressed exponential fits;
the empirical decorrelation length ℓ_c(q) = R ω τ_c. Part II sketches applications in
random-feature optical computing, information transmission, and random-transient-grating
thermal transport.

## Relation to the July 24 2026 data

The July 24 run measured four angles (0°, 10°, 15°, 20°) at six rotation periods each,
which directly supports the τ_c ∝ 1/ω and sinθ scaling tests. The notebook observation
that decay persists at nominal 0° is exactly the case the ideal model says should *not*
decorrelate, so it points to alignment error, wobble, or boiling.

Note that the two documents use conflicting angle symbols (θ vs φ) and different scattering
vector conventions; reconcile before merging any of this into a paper.
