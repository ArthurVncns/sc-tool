from pathlib import Path

import streamlit as st

from analysis import io
from ui import state
from ui.theme import apply_theme, next_page_button

st.set_page_config(page_title="sc_tool — Upload", layout="wide")
apply_theme()

st.title("Upload Dataset")

# ── Format selector ────────────────────────────────────────────────────────────
fmt = st.radio(
    "File format",
    options=[
        "AnnData (.h5ad)",
        "10x HDF5 (.h5)",
        "10x MTX folder",
        "Seurat RDS (.rds)",
    ],
    horizontal=True,
    key="upload_format",
)

st.divider()

# ── Format-specific loading ────────────────────────────────────────────────────

def _store(adata, file_key: str) -> None:
    """Save adata to session state and record the file key."""
    state.set_adata(adata)
    state.set_file_key(file_key)


# ── AnnData (.h5ad) ────────────────────────────────────────────────────────────
if fmt == "AnnData (.h5ad)":
    uploaded = st.file_uploader("Upload .h5ad file", type=["h5ad"])
    if uploaded is not None:
        file_key = f"h5ad_{uploaded.name}_{uploaded.size}"
        if state.get_file_key() != file_key:
            with st.spinner(f"Loading {uploaded.name}…"):
                _store(io.load_h5ad(uploaded), file_key)

# ── 10x HDF5 (.h5) ────────────────────────────────────────────────────────────
elif fmt == "10x HDF5 (.h5)":
    st.caption("Cell Ranger HDF5 output — typically named `filtered_feature_bc_matrix.h5`.")
    uploaded = st.file_uploader("Upload .h5 file", type=["h5"])
    if uploaded is not None:
        file_key = f"h5_{uploaded.name}_{uploaded.size}"
        if state.get_file_key() != file_key:
            with st.spinner(f"Loading {uploaded.name}…"):
                try:
                    _store(io.load_h5(uploaded), file_key)
                except Exception as e:
                    st.error(f"Failed to load H5 file: {e}")

# ── 10x MTX folder ────────────────────────────────────────────────────────────
elif fmt == "10x MTX folder":
    st.caption(
        "Folder must contain `matrix.mtx.gz`, `barcodes.tsv.gz`, and "
        "`features.tsv.gz` (CellRanger 3+) or `genes.tsv.gz` (CellRanger 2.x)."
    )
    folder_input = st.text_input(
        "Path to 10x folder",
        placeholder="/path/to/filtered_feature_bc_matrix",
    )
    if folder_input:
        file_key = f"mtx_{folder_input}"
        if state.get_file_key() != file_key:
            if st.button("Load", type="primary"):
                with st.spinner("Loading 10x MTX…"):
                    try:
                        _store(io.load_10x_mtx(Path(folder_input)), file_key)
                    except (FileNotFoundError, ValueError) as e:
                        st.error(str(e))

# ── Seurat RDS (.rds) ──────────────────────────────────────────────────────────
elif fmt == "Seurat RDS (.rds)":
    st.caption(
        "Converts a Seurat object to AnnData via **rpy2**. "
        "Requires R and the Seurat package to be installed on this machine. "
        "Extracts raw counts, cell metadata, and any stored PCA / UMAP embeddings."
    )
    st.info(
        "**Requirements:** `pip install rpy2` and, in R: `install.packages('Seurat')`"
    )
    uploaded = st.file_uploader("Upload .rds file", type=["rds"])
    if uploaded is not None:
        file_key = f"rds_{uploaded.name}_{uploaded.size}"
        if state.get_file_key() != file_key:
            with st.spinner(f"Converting {uploaded.name} — this may take a moment…"):
                try:
                    _store(io.load_seurat_rds(uploaded), file_key)
                except ImportError as e:
                    st.error(str(e))
                except ValueError as e:
                    st.error(str(e))
                except Exception as e:
                    st.error(f"Conversion failed: {e}")

# ── Dataset summary (shown for all formats once data is loaded) ────────────────
if state.has_adata():
    adata = state.get_adata()

    st.divider()
    st.subheader("Dataset Overview")

    col1, col2 = st.columns(2)
    col1.metric("Cells (obs)", f"{adata.n_obs:,}")
    col2.metric("Genes (var)", f"{adata.n_vars:,}")

    col_obs, col_var = st.columns(2)

    with col_obs:
        st.subheader("Cell metadata (obs)")
        if len(adata.obs.columns) > 0:
            st.dataframe(adata.obs.head(20), use_container_width=True)
        else:
            st.info("No cell metadata columns found.")

    with col_var:
        st.subheader("Gene metadata (var)")
        if len(adata.var.columns) > 0:
            st.dataframe(adata.var.head(20), use_container_width=True)
        else:
            st.info("No gene metadata columns found.")

    if adata.obsm or adata.uns:
        st.subheader("Stored results")
        col_obsm, col_uns = st.columns(2)
        with col_obsm:
            if adata.obsm:
                st.write("**Embeddings (obsm):**", list(adata.obsm.keys()))
        with col_uns:
            if adata.uns:
                st.write("**Unstructured metadata (uns):**", list(adata.uns.keys()))

    next_page_button("QC", "pages/02_QC.py")
