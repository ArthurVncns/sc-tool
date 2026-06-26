import streamlit as st

from ui import state
from ui.theme import apply_theme, next_page_button

st.set_page_config(page_title="sc_tool", layout="wide")
apply_theme()

st.title("sc_tool")
st.caption("Interactive single-cell RNA-seq analysis")

st.markdown("""
**Workflow**

1. **Upload** — load an `.h5ad` dataset
2. **QC** — compute quality metrics and filter low-quality cells
3. **Preprocessing** — normalization and HVG selection
4. **Dimensionality reduction** — PCA and UMAP
5. **Annotation** — cluster cells, find marker genes, assign cell types

Use the sidebar to navigate between steps.
""")

if state.has_adata():
    adata = state.get_adata()
    st.divider()
    st.subheader("Current dataset")

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Cells", f"{adata.n_obs:,}")
    col2.metric("Genes", f"{adata.n_vars:,}")
    col3.metric("QC computed", "Yes" if state.qc_is_computed() else "No")
    col3.metric("Preprocessed", "Yes" if state.hvg_selected() else "No")
    col4.metric("UMAP", "Yes" if state.umap_done() else "No")

    if state.annotation_done():
        n_types = adata.obs["cell_type"].nunique()
    elif state.celltypist_done():
        n_types = adata.obs["celltypist_cell_type"].nunique()
    elif state.marker_score_done():
        n_types = adata.obs["marker_score_cell_type"].nunique()
    else:
        n_types = None
    col4.metric("Cell types", str(n_types) if n_types is not None else "—")

    if not state.qc_is_computed():
        next_page_button("QC", "pages/02_QC.py")
    elif not state.hvg_selected():
        next_page_button("Preprocessing", "pages/03_Preprocessing.py")
    elif not state.umap_done():
        next_page_button("Dimensionality Reduction", "pages/04_Reduction.py")
    elif not state.any_annotation_done():
        next_page_button("Annotation", "pages/05_Annotation.py")
else:
    st.info("No dataset loaded. Start by uploading a file on the **Upload** page.")
    next_page_button("Upload", "pages/01_Upload.py")
