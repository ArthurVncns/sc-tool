import streamlit as st

from ui import state

st.set_page_config(page_title="sc_tool", layout="wide")

st.title("sc_tool")
st.caption("Interactive single-cell RNA-seq analysis")

st.markdown("""
**Workflow**

1. **Upload** — load an `.h5ad` dataset
2. **QC** — compute quality metrics and filter low-quality cells
3. **Preprocessing** — normalize, log-transform, select variable genes *(coming soon)*
4. **Dimensionality reduction** — PCA and UMAP *(coming soon)*
5. **Clustering** — identify cell populations *(coming soon)*

Use the sidebar to navigate between steps.
""")

if state.has_adata():
    adata = state.get_adata()
    st.divider()
    st.subheader("Current dataset")
    col1, col2, col3 = st.columns(3)
    col1.metric("Cells", f"{adata.n_obs:,}")
    col2.metric("Genes", f"{adata.n_vars:,}")
    col3.metric("QC computed", "Yes" if state.qc_is_computed() else "No")
else:
    st.info("No dataset loaded. Start by uploading a file on the **Upload** page.")
