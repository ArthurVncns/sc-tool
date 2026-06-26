# sc_tool

An interactive web application for single-cell RNA-seq analysis, built with Python and Streamlit.

Upload an `.h5ad` dataset and explore it step-by-step: quality control, preprocessing, dimensionality reduction, clustering, and cell type annotation — all through a point-and-click interface, no code required.

---

## Features

| Step | What it does |
|---|---|
| **Upload** | Load any `.h5ad` file (AnnData format from Scanpy, Cell Ranger, etc.) |
| **QC** | Per-cell metrics (genes, counts, mt / ribo / hb %) with interactive filtering sliders and live cell-count preview |
| **Preprocessing** | Two mutually exclusive workflows: Standard (CPM + log1p) or Pearson Residuals (scTransform-like) |
| **HVG selection** | Seurat dispersion, LOESS / Seurat v3 (raw counts), or Pearson residuals variance ranking |
| **Dimensionality reduction** | PCA with elbow plot → UMAP with tunable neighbors and min-distance |
| **Clustering** | Leiden algorithm with adjustable resolution; UMAP colored by cluster |
| **Marker genes** | Wilcoxon rank-sum test per cluster; top-N marker table |
| **Cell type annotation** | Three methods: Manual (cluster-by-cluster), CellTypist (pre-trained models), Marker gene scoring (custom JSON dictionary) |

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

61 tests across IO, QC, preprocessing, dimensionality reduction, and annotation.

---

## Project structure

```
sc_tool/
├── app.py                    # Home page and entry point
├── pages/
│   ├── 01_Upload.py          # Dataset upload and overview
│   ├── 02_QC.py              # Quality control and cell filtering
│   ├── 03_Preprocessing.py   # Normalization, log transform, HVG selection
│   ├── 04_Reduction.py       # PCA and UMAP
│   └── 05_Annotation.py      # Clustering, marker genes, cell type annotation
├── analysis/                 # Pure Python analysis logic (no Streamlit)
│   ├── io.py                 # AnnData loading
│   ├── qc.py                 # QC metrics and filtering
│   ├── preprocessing.py      # Normalization, log transform, HVG
│   ├── reduction.py          # PCA and UMAP
│   └── annotation.py         # Clustering, marker genes, cell type annotation
├── ui/
│   ├── state.py              # Centralized session_state helpers
│   ├── plots.py              # Plotly figure builders
│   └── theme.py              # Glass UI theme (CSS injection + navigation)
└── tests/
```

The `analysis/` layer is fully independent of Streamlit and can be imported into notebooks or scripts directly.

---

## Preprocessing workflows

### Workflow 1 — Standard Scanpy
Normalize counts per cell → log1p → select highly variable genes.
- HVG options: Seurat (dispersion on log-normalized data) or LOESS / Seurat v3 (variance on raw counts)

### Workflow 2 — Pearson Residuals
Fit a regularized negative binomial model per gene → replace counts with Pearson residuals → select HVGs ranked by residual variance.
- More robust for datasets with large differences in sequencing depth.
- HVG selection uses raw counts via `adata.layers['counts']`.
- CellTypist annotation is not available with this workflow (requires log-normalized data).

---

## Cell type annotation methods

| Method | How it works | Best for |
|---|---|---|
| **Manual** | Enter a cell type name per cluster based on marker genes | When you know your tissue well |
| **CellTypist** | Pre-trained logistic regression models (auto-downloaded) | PBMC, lung, brain, and other common tissues |
| **Marker gene scoring** | Score cells against a custom JSON marker dictionary | Any species or tissue without a pre-trained model |

For PBMC data, use CellTypist with the **Immune — broad** model and majority voting enabled.

---

## Tech stack

- [Streamlit](https://streamlit.io/) — UI framework
- [Scanpy](https://scanpy.readthedocs.io/) — single-cell analysis algorithms
- [AnnData](https://anndata.readthedocs.io/) — data container
- [Plotly](https://plotly.com/python/) — interactive visualizations
- [CellTypist](https://www.celltypist.org/) — automated cell type classification
- [scikit-misc](https://has2k1.github.io/scikit-misc/) — LOESS regression for Seurat v3 HVG
- [igraph](https://python.igraph.org/) — Leiden clustering
