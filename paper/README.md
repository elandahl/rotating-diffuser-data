# Paper draft — Overleaf package

This folder is a self-contained LaTeX project for the rotating-diffuser
autocorrelation draft. Upload the **contents of this `paper/` directory**
(or the whole folder as a zip) to Overleaf.

## Contents

| File | Role |
|------|------|
| `main.tex` | Article source |
| `refs.bib` | Bibliography |
| `figures/` | PDF figures included by `main.tex` |
| `make_figures.py` | Regenerates figures from `../analysis/` (optional; not needed on Overleaf) |

## Open in Overleaf

1. Zip this directory, e.g. from the repo root:
   ```bash
   cd paper && zip -r ../rotating-diffuser-paper.zip main.tex refs.bib figures/*.pdf
   ```
2. Overleaf → New Project → Upload Project → select the zip.
3. Set the main document to `main.tex` if prompted.
4. Recompile. Bibliography uses `natbib` + `unsrtnat`.

Or push this repo to GitHub and use Overleaf → New Project → Import from GitHub, then open the `paper/` folder as the project root (or set Overleaf’s root to `paper/main.tex` if the whole repo is imported).

## Local compile

With [Tectonic](https://tectonic-typesetting.github.io/):

```bash
cd paper
tectonic -X compile main.tex
```

With latexmk (TeX Live / MacTeX):

```bash
cd paper
latexmk -pdf main.tex
```

## Regenerating figures

From the repository root (requires the analysis CSVs already present):

```bash
python3 paper/make_figures.py
```

Raw `test*.txt` data are never modified; figures are derived products only.

## Scope of this draft

The text reports only what the July 24 2026 data set currently supports:
erfc(ln τ) midpoint times, \(1/\tau_{50}\propto\omega\) at fixed angle, and
slopes growing roughly as \(\sin\theta\) for \(10^\circ\)–\(20^\circ\). It is written
for an advanced undergraduate physics audience.
