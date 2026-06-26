import streamlit as st

from analysis import io
from ui import state
from ui.theme import apply_theme, next_page_button

st.set_page_config(page_title="sc_tool — Upload", layout="wide")
apply_theme()

st.title("Upload Dataset")

uploaded_file = st.file_uploader(
    "Upload an AnnData file (.h5ad)",
    type=["h5ad"],
    help="Standard AnnData format produced by Scanpy, Cell Ranger, etc.",
)

if uploaded_file is not None:
    # name + size together identify a unique file well enough for v0.1.
    # This avoids reloading the same file on every Streamlit rerun.
    file_key = f"{uploaded_file.name}_{uploaded_file.size}"

    if state.get_file_key() != file_key:
        with st.spinner(f"Loading {uploaded_file.name}..."):
            adata = io.load_h5ad(uploaded_file)
            state.set_adata(adata)
            state.set_file_key(file_key)

if state.has_adata():
    adata = state.get_adata()

    st.divider()
    st.header("Dataset Overview")

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
