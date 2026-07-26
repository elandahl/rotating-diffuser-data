# Paper draft — Overleaf package

Self-contained LaTeX project for the rotating-diffuser autocorrelation draft.
Upload **this `paper/` directory** to Overleaf, or import the GitHub repo and
set the project root to `paper/`.

## Contents

| File | Role |
|------|------|
| `main.tex` | Article source |
| `refs.bib` | Bibliography |
| `figures/*.pdf` | Figures included by `main.tex` |
| `main.pdf` | Local compile preview (optional on Overleaf) |

Figure generation lives outside this folder: `../analysis/make_paper_figures.py`.

## Open in Overleaf

**GitHub import (preferred):** New Project → Import from GitHub →
`elandahl/rotating-diffuser-data`, then set the root document to `paper/main.tex`.

**Zip upload:**
```bash
cd paper && zip -r ../rotating-diffuser-paper.zip main.tex refs.bib figures/*.pdf README.md
```

Recompile with `natbib` + `unsrtnat` (default in `main.tex`).

## Local compile

```bash
cd paper
tectonic -X compile main.tex
# or: latexmk -pdf main.tex
```

## Regenerating figures

From the repository root (requires `analysis/erfc_fits/erfc_fit_parameters.csv`):

```bash
python3 analysis/make_paper_figures.py
cd paper && tectonic -X compile main.tex
```

Raw `test*.txt` data are never modified.

## Authors

Gantulga Gankhuyag and E. Landahl — draft based on data acquired 24 July 2026.
