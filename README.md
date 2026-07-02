# sc_tool

An interactive web application for single-cell RNA-seq analysis, built with Python and Streamlit.

Upload a dataset and explore it step-by-step: quality control, doublet detection, preprocessing, batch correction, dimensionality reduction, clustering, cell type annotation, and differential expression — all through a point-and-click interface, no code required.

---

## Features

| Step | What it does |
|---|---|
| **Upload** | Load `.h5ad`, 10x HDF5 (`.h5`), 10x MTX folder, or Seurat RDS files |
| **Quality Control** | Per-cell metrics (genes, counts, mt / ribo / hb %), doublet detection (Scrublet), interactive filtering with live preview |
| **Preprocessing** | Two workflows: Standard (CPM + log1p) or Pearson Residuals; HVG selection via Seurat, LOESS / Seurat v3, or residual variance |
| **Dimensionality Reduction** | PCA with elbow plot, optional Harmony batch correction, UMAP with tunable neighbors and min-distance |
| **Annotation** | Leiden clustering, Wilcoxon marker genes, cell type assignment via Manual / CellTypist / Marker gene scoring |
| **Differential Expression** | Wilcoxon or t-test, volcano plot, filtered gene table, CSV export |
| **Export** | Download processed `.h5ad`, cell/gene metadata CSVs, DE results, and a plain-text analysis summary |

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

86 tests across IO, QC, preprocessing, dimensionality reduction, annotation, and differential expression.

---

## Project structure

```
sc_tool/
├── app.py                      # Navigation controller (st.navigation)
├── pages/
│   ├── home.py                 # Overview page
│   ├── 01_Upload.py            # Dataset upload (h5ad, h5, MTX, RDS)
│   ├── 02_QC.py                # Quality control, doublet detection, filtering
│   ├── 03_Preprocessing.py     # Normalization, log transform, HVG selection
│   ├── 04_Reduction.py         # PCA, Harmony, UMAP
│   ├── 05_Annotation.py        # Clustering, marker genes, cell type annotation
│   ├── 06_DE.py                # Differential expression
│   └── 07_Export.py            # Data and summary export
├── analysis/                   # Pure Python analysis logic (no Streamlit)
│   ├── io.py                   # Multi-format data loading
│   ├── qc.py                   # QC metrics and cell filtering
│   ├── doublets.py             # Scrublet doublet detection
│   ├── preprocessing.py        # Normalization, log transform, HVG
│   ├── reduction.py            # PCA, Harmony, UMAP
│   ├── annotation.py           # Clustering, marker genes, annotation
│   └── de.py                   # Differential expression
├── ui/
│   ├── state.py                # Centralized session_state helpers
│   ├── plots.py                # Plotly figure builders
│   └── theme.py                # Glass UI theme and navigation buttons
└── tests/
```

The `analysis/` layer is fully independent of Streamlit and can be imported into notebooks or scripts directly.

---

## Preprocessing workflows

### Workflow 1 — Standard Scanpy
Normalize counts per cell (CPM) → log1p → select highly variable genes.
- HVG options: Seurat (dispersion) or LOESS / Seurat v3 (variance on raw counts)

### Workflow 2 — Pearson Residuals
Fit a regularized negative binomial model per gene → Pearson residuals → HVGs ranked by residual variance.
- More robust for datasets with large differences in sequencing depth.
- CellTypist annotation requires Workflow 1 (log-normalized data).

---

## Cell type annotation methods

| Method | How it works | Best for |
|---|---|---|
| **Manual** | Enter a cell type name per cluster based on marker genes | When you know your tissue |
| **CellTypist** | Pre-trained logistic regression models, auto-downloaded | PBMC, lung, brain, and other common tissues |
| **Marker gene scoring** | Score cells against a custom JSON marker dictionary | Any species or tissue |

For PBMC data, use CellTypist with the **Immune — broad** model and majority voting enabled.

---

## Supported input formats

| Format | How to load |
|---|---|
| AnnData `.h5ad` | File uploader |
| 10x Genomics HDF5 `.h5` | File uploader |
| 10x Genomics MTX | Local folder path (contains `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`) |
| Seurat RDS | File uploader — requires R, Seurat, and `pip install rpy2` |

---

## Tech stack

- [Streamlit](https://streamlit.io/) — UI framework
- [Scanpy](https://scanpy.readthedocs.io/) — single-cell analysis algorithms
- [AnnData](https://anndata.readthedocs.io/) — data container
- [Plotly](https://plotly.com/python/) — interactive visualizations
- [CellTypist](https://www.celltypist.org/) — automated cell type classification
- [harmonypy](https://github.com/slowkow/harmonypy) — batch correction
- [scikit-misc](https://has2k1.github.io/scikit-misc/) — LOESS regression for Seurat v3 HVG
- [igraph](https://python.igraph.org/) — Leiden clustering
