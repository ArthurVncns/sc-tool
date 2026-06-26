# sc_tool

An interactive web application for single-cell RNA-seq analysis, built with Python and Streamlit.

Upload an `.h5ad` dataset and explore it step-by-step: quality control, preprocessing, dimensionality reduction, and visualization — all through a point-and-click interface, no code required.

---

## Features

| Step | What it does |
|---|---|
| **Upload** | Load any `.h5ad` file (AnnData format from Scanpy, Cell Ranger, etc.) |
| **QC** | Compute per-cell metrics (genes, counts, mitochondrial / ribosomal / hemoglobin %) and filter low-quality cells with interactive sliders |
| **Preprocessing** | Two coherent workflows: Standard (CPM + log1p) or Pearson Residuals (scTransform-like) |
| **HVG selection** | Seurat dispersion, LOESS / Seurat v3, or Pearson residuals variance ranking |
| **Dimensionality reduction** | PCA (elbow plot) → UMAP (tunable neighbors and min-distance) |

---

## Requirements

- Python 3.12+
- [uv](https://docs.astral.sh/uv/) (recommended) or pip

---

## Setup

```bash
# Clone the repository
git clone <repo-url>
cd sc_tool

# Install dependencies
uv sync --extra dev
```

---

## Run

```bash
streamlit run app.py
```

The app opens in your browser at `http://localhost:8501`.

---

## Test

```bash
pytest
```

50 tests across IO, QC, preprocessing, and dimensionality reduction.

---

## Project structure

```
sc_tool/
├── app.py                   # Home page and entry point
├── pages/
│   ├── 01_Upload.py         # Dataset upload and overview
│   ├── 02_QC.py             # Quality control and cell filtering
│   ├── 03_Preprocessing.py  # Normalization, log transform, HVG selection
│   └── 04_Reduction.py      # PCA and UMAP
├── analysis/                # Pure Python analysis logic (no Streamlit)
│   ├── io.py                # AnnData loading
│   ├── qc.py                # QC metrics and filtering
│   ├── preprocessing.py     # Normalization, log transform, HVG
│   └── reduction.py         # PCA and UMAP
├── ui/
│   ├── state.py             # Centralized session_state helpers
│   └── plots.py             # Plotly figure builders
└── tests/
```

The `analysis/` layer is fully independent of Streamlit and can be imported into notebooks or scripts directly.

---

## Preprocessing workflows

### Workflow 1 — Standard Scanpy
Normalize counts per cell → log1p → select highly variable genes.
- HVG options: Seurat (dispersion) or LOESS / Seurat v3 (variance on raw counts)

### Workflow 2 — Pearson Residuals
Fit a regularized negative binomial model per gene → replace counts with Pearson residuals → select HVGs ranked by residual variance.
- More robust for datasets with large differences in sequencing depth.
- HVG selection uses raw counts via `adata.layers['counts']`, not the normalized residuals.

---

## Tech stack

- [Streamlit](https://streamlit.io/) — UI framework
- [Scanpy](https://scanpy.readthedocs.io/) — single-cell analysis algorithms
- [AnnData](https://anndata.readthedocs.io/) — data container
- [Plotly](https://plotly.com/python/) — interactive visualizations
- [scikit-misc](https://has2k1.github.io/scikit-misc/) — LOESS regression for Seurat v3 HVG
