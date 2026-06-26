"""Centralized helpers for reading and writing Streamlit session state.

All session_state keys live here so there is one authoritative place to see
every piece of state the app holds. This prevents typos and makes future
refactoring easier.
"""

import anndata as ad
import streamlit as st

_ADATA_KEY = "adata"
_FILE_KEY = "uploaded_file_key"


def get_adata() -> ad.AnnData | None:
    """Return the currently loaded AnnData, or None if none has been loaded."""
    return st.session_state.get(_ADATA_KEY)


def set_adata(adata: ad.AnnData) -> None:
    st.session_state[_ADATA_KEY] = adata


def has_adata() -> bool:
    return _ADATA_KEY in st.session_state


def get_file_key() -> str | None:
    """Return the key identifying the last loaded file, or None."""
    return st.session_state.get(_FILE_KEY)


def set_file_key(key: str) -> None:
    st.session_state[_FILE_KEY] = key


def qc_is_computed() -> bool:
    """Return True if QC metrics have already been added to the current adata."""
    adata = get_adata()
    return adata is not None and "n_genes_by_counts" in adata.obs.columns


def mt_detected() -> bool:
    """Return True if mitochondrial QC metrics were computed for the current adata."""
    adata = get_adata()
    return adata is not None and "pct_counts_mt" in adata.obs.columns


def ribo_detected() -> bool:
    """Return True if ribosomal QC metrics were computed for the current adata."""
    adata = get_adata()
    return adata is not None and "pct_counts_ribo" in adata.obs.columns


def hb_detected() -> bool:
    """Return True if hemoglobin QC metrics were computed for the current adata."""
    adata = get_adata()
    return adata is not None and "pct_counts_hb" in adata.obs.columns


# --- Preprocessing state ---

def is_normalized() -> bool:
    """Return True if normalization has been applied (either workflow).

    Checks adata.layers['counts'], written by both normalize_cpm() and
    normalize_pearson() before modifying adata.X.
    """
    adata = get_adata()
    return adata is not None and "counts" in adata.layers


def active_workflow() -> str | None:
    """Return 'standard', 'pearson', or None if preprocessing has not started."""
    adata = get_adata()
    if adata is None:
        return None
    return adata.uns.get("sc_tool_workflow")


def is_log_transformed() -> bool:
    """Return True if log1p has been applied (Workflow 1 only)."""
    adata = get_adata()
    return adata is not None and "log1p" in adata.uns


def ready_for_hvg() -> bool:
    """Return True when the data is ready for HVG selection.

    - Workflow 1 (standard): requires log transform to be done.
    - Workflow 2 (pearson):  normalization alone is sufficient.
    """
    workflow = active_workflow()
    if workflow == "standard":
        return is_log_transformed()
    if workflow == "pearson":
        return is_normalized()
    return False


def hvg_selected() -> bool:
    """Return True if highly variable genes have been selected."""
    adata = get_adata()
    return adata is not None and "highly_variable" in adata.var.columns


# --- Dimensionality reduction state ---

def pca_done() -> bool:
    """Return True if PCA has been computed."""
    adata = get_adata()
    return adata is not None and "X_pca" in adata.obsm


def umap_done() -> bool:
    """Return True if a UMAP embedding has been computed."""
    adata = get_adata()
    return adata is not None and "X_umap" in adata.obsm


# --- Annotation state ---

def clustering_done() -> bool:
    """Return True if Leiden clustering has been run."""
    adata = get_adata()
    return adata is not None and "leiden" in adata.obs.columns


def markers_computed() -> bool:
    """Return True if marker genes have been computed."""
    adata = get_adata()
    return adata is not None and "rank_genes_groups" in adata.uns


def annotation_done() -> bool:
    """Return True if manual cell type annotations have been applied."""
    adata = get_adata()
    return adata is not None and "cell_type" in adata.obs.columns


def celltypist_done() -> bool:
    """Return True if CellTypist annotation has been run."""
    adata = get_adata()
    return adata is not None and "celltypist_cell_type" in adata.obs.columns


def marker_score_done() -> bool:
    """Return True if marker gene scoring annotation has been run."""
    adata = get_adata()
    return adata is not None and "marker_score_cell_type" in adata.obs.columns


def any_annotation_done() -> bool:
    """Return True if any annotation method has produced results."""
    return annotation_done() or celltypist_done() or marker_score_done()
