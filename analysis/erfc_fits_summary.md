# erfc(ln τ) batch fits

Model: \(g_2(\tau)=A\,\mathrm{erfc}\big((\ln\tau-B)/C_w\big)+D\)

- \(\tau_{50}=e^{B}\) (midpoint on the log-τ sigmoid)
- \(\beta=2A\) (full step height; high plateau \(\approx D+2A\))
- Skip first 10 lags (afterpulsing); fit window τ ≤ 0.1 s
- 1σ uncertainties from residual bootstrap (400 resamples); formal OLS errors are much smaller and ignore systematic residual structure
- `test03`: first Octave dump only; no period in filename

## Parameters

| File | Angle | T_rot (s) | τ₅₀ (ms) | β | C_w | D | RMS | Plot |
|---|---|---|---|---|---|---|---|---|
| `test016_3.78s_00deg` | small angle (nominal 0°) | 3.78 | 3.320±0.046 | 0.1122±0.0007 | 1.805±0.027 | 0.9956±0.0005 | 0.0022 | `test016_3.78s_00deg_erfc.png` |
| `test017_6.31s_00deg` | small angle (nominal 0°) | 6.31 | 5.441±0.082 | 0.1100±0.0007 | 1.792±0.028 | 0.9971±0.0006 | 0.0021 | `test017_6.31s_00deg_erfc.png` |
| `test018_9.77s_00deg` | small angle (nominal 0°) | 9.77 | 9.335±0.153 | 0.1167±0.0009 | 1.921±0.027 | 0.9938±0.0008 | 0.0017 | `test018_9.77s_00deg_erfc.png` |
| `test019_14.23s_00deg` | small angle (nominal 0°) | 14.23 | 13.872±0.317 | 0.1123±0.0012 | 1.988±0.032 | 0.9961±0.0011 | 0.0016 | `test019_14.23s_00deg_erfc.png` |
| `test020_30.93s_00deg` | small angle (nominal 0°) | 30.93 | 16.935±0.164 | 0.1077±0.0006 | 1.616±0.015 | 1.0055±0.0005 | 0.0009 | `test020_30.93s_00deg_erfc.png` |
| `test021_33.59s_00deg` | small angle (nominal 0°) | 33.59 | 33.705±0.870 | 0.1172±0.0018 | 1.705±0.027 | 0.9932±0.0017 | 0.0010 | `test021_33.59s_00deg_erfc.png` |
| `test011_3.75s_10deg` | 10° | 3.75 | 1.007±0.008 | 0.1871±0.0007 | 0.882±0.013 | 1.0016±0.0004 | 0.0037 | `test011_3.75s_10deg_erfc.png` |
| `test010_6.80s_10deg` | 10° | 6.80 | 1.778±0.014 | 0.1895±0.0008 | 0.925±0.016 | 1.0008±0.0005 | 0.0040 | `test010_6.80s_10deg_erfc.png` |
| `test012_10.90s_10deg` | 10° | 10.90 | 2.986±0.021 | 0.1922±0.0008 | 0.923±0.014 | 0.9987±0.0005 | 0.0038 | `test012_10.90s_10deg_erfc.png` |
| `test014_15.51s_10deg` | 10° | 15.51 | 4.252±0.029 | 0.1962±0.0007 | 0.934±0.012 | 1.0025±0.0005 | 0.0034 | `test014_15.51s_10deg_erfc.png` |
| `test015_19.93s_10deg` | 10° | 19.93 | 5.173±0.040 | 0.1837±0.0008 | 0.935±0.014 | 0.9946±0.0006 | 0.0035 | `test015_19.93s_10deg_erfc.png` |
| `test013_46.73s_10deg` | 10° | 46.73 | 11.074±0.088 | 0.1812±0.0009 | 0.967±0.014 | 0.9908±0.0007 | 0.0032 | `test013_46.73s_10deg_erfc.png` |
| `test022_3.85s_15deg` | 15° | 3.85 | 0.659±0.004 | 0.2066±0.0008 | 0.884±0.014 | 1.0047±0.0004 | 0.0035 | `test022_3.85s_15deg_erfc.png` |
| `test023_4.65s_15deg` | 15° | 4.65 | 0.786±0.005 | 0.2069±0.0007 | 0.901±0.012 | 1.0047±0.0004 | 0.0034 | `test023_4.65s_15deg_erfc.png` |
| `test024_5.88s_15deg` | 15° | 5.88 | 1.030±0.007 | 0.2063±0.0007 | 0.949±0.012 | 1.0047±0.0004 | 0.0036 | `test024_5.88s_15deg_erfc.png` |
| `test025_7.81s_15deg` | 15° | 7.81 | 1.305±0.009 | 0.2031±0.0007 | 0.937±0.013 | 1.0044±0.0004 | 0.0036 | `test025_7.81s_15deg_erfc.png` |
| `test026_11.45s_15deg` | 15° | 11.45 | 2.003±0.014 | 0.2072±0.0007 | 0.950±0.012 | 1.0033±0.0005 | 0.0034 | `test026_11.45s_15deg_erfc.png` |
| `test027_17.50s_15deg` | 15° | 17.50 | 3.039±0.021 | 0.2033±0.0007 | 0.907±0.014 | 1.0050±0.0005 | 0.0037 | `test027_17.50s_15deg_erfc.png` |
| `test027_37.33s_15deg` | 15° | 37.33 | 6.032±0.051 | 0.2019±0.0010 | 0.919±0.016 | 1.0025±0.0008 | 0.0043 | `test027_37.33s_15deg_erfc.png` |
| `test08_3.03s` | 20° (unlabeled) | 3.03 | 0.379±0.002 | 0.2470±0.0007 | 0.906±0.013 | 1.0010±0.0003 | 0.0031 | `test08_3.03s_erfc.png` |
| `test04_5.70s` | 20° (unlabeled) | 5.70 | 0.713±0.003 | 0.2404±0.0006 | 0.927±0.009 | 0.9993±0.0003 | 0.0028 | `test04_5.70s_erfc.png` |
| `test05_6.30s` | 20° (unlabeled) | 6.30 | 0.807±0.003 | 0.2467±0.0006 | 0.922±0.008 | 0.9998±0.0003 | 0.0027 | `test05_6.30s_erfc.png` |
| `test06_10.74s` | 20° (unlabeled) | 10.74 | 1.357±0.007 | 0.2407±0.0006 | 0.927±0.009 | 0.9998±0.0004 | 0.0029 | `test06_10.74s_erfc.png` |
| `test07_29.09s` | 20° (unlabeled) | 29.09 | 3.893±0.023 | 0.2430±0.0007 | 0.956±0.011 | 1.0002±0.0005 | 0.0036 | `test07_29.09s_erfc.png` |
| `test09_47.59s` | 20° (unlabeled) | 47.59 | 6.214±0.043 | 0.2602±0.0010 | 0.868±0.013 | 0.9975±0.0008 | 0.0046 | `test09_47.59s_erfc.png` |
| `test03` | 20° (unlabeled) | — | 0.718±0.003 | 0.2455±0.0006 | 0.893±0.009 | 0.9991±0.0003 | 0.0028 | `test03_erfc.png` |

## Overviews

- `erfc_fits/overview_tau50_vs_period.png` — τ₅₀ vs T_rot
- `erfc_fits/overview_normalized_overlays.png` — contrast-normalized data+fits

Full numeric table: `erfc_fits/erfc_fit_parameters.csv`
