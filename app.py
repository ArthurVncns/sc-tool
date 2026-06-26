import streamlit as st

from ui import state

st.set_page_config(page_title="sc_tool", layout="wide")

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
    # Show cell type count from whichever annotation method was run
    if state.annotation_done():
        n_types = adata.obs["cell_type"].nunique()
    elif state.celltypist_done():
        n_types = adata.obs["celltypist_cell_type"].nunique()
    elif state.marker_score_done():
        n_types = adata.obs["marker_score_cell_type"].nunique()
    else:
        n_types = None
    col4.metric("Cell types", str(n_types) if n_types is not None else "—")
else:
    st.info("No dataset loaded. Start by uploading a file on the **Upload** page.")
