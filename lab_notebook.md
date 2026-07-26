# Lab notebook — rotating ground-glass diffuser autocorrelations

**Date:** July 24, 2026  
**Operators:** GG and EL

---

## Setup

Rotating diffuser (??? grit) placed in Gana’s rotating mount and using diode laser **3.31 V, 0.033 A**, iris adjusted to keep counts near **1 MHz**. Diffuser plane placed at diffraction plane / diffractometer center of rotation. Autocorrelations acquired at 4 different angles: **20°, 15°, 10°, 0°**.

Note: laser beam was focused approx into the diffuser plate like before with DLS from liquid samples.

### Apparatus photos

*Add apparatus photos here — top and side views.*

---

## Acquisition

Ten second acquisition times were used for each and appeared sufficient at these count rates to get fairly smooth curves.

We used a simplified variety of the multiangle_multitau.m program called **simple_multitau.m**, which simply acquires an autocorrelation trace and saves the data and an image on semilog-x axes.

---

## Rotation periods

The period of rotation was measured using a stopwatch for one revolution at a series of different rotation drives.

Period is written in the filename.

**Example:** `test025_7.81s_15deg.txt` → period **7.81 s**, **15°** off the direct forward beam.

We took **4 angles** and **6 rotation speeds** each. The unlabeled angles are all **20°**.

---

## Observations

- Correlation time gets smaller (faster) as the period goes down (faster rotation).
- Seems most of the time to also get faster with higher angle.
- An exception is “zero” degrees, where in fact we may not have been perfectly aligned with the very small iris required for low count rates exactly on-axis.
