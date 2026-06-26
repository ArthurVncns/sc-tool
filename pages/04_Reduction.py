import pandas as pd
import streamlit as st

from analysis import reduction
from ui import plots, state

st.set_page_config(page_title="sc_tool — Reduction", layout="wide")

st.title("Dimensionality Reduction")

# --- Guards ---
if not state.has_adata():
    st.info("Please upload a dataset on the **Upload** page first.")
    st.stop()

if not state.hvg_selected():
    st.warning("Please complete **Preprocessing** (including HVG selection) before running PCA.")
    st.stop()

adata = state.get_adata()
st.caption(f"Current dataset: **{adata.n_obs:,} cells** × {adata.n_vars:,} genes")

# --- Feedback message ---
if "_reduction_result" in st.session_state:
    st.success(st.session_state.pop("_reduction_result"))

# =========================================================================
# Step 1 — PCA
# =========================================================================
st.subheader("Step 1 — PCA")
st.write(
    "Reduces the High Variable Genes dimensions to a smaller set of principal components "
    "that capture the main axes of variation. Use the elbow plot below to choose "
    "how many components contain meaningful signal."
)

n_hvg = int(adata.var["highly_variable"].sum())
max_comps = max(10, min(adata.n_obs - 1, n_hvg - 1, 100))

n_comps = st.slider(
    "Number of principal components",
    min_value=10,
    max_value=max_comps,
    value=min(50, max_comps),
    step=5,
    help=f"Capped at min(n_cells − 1, n_HVGs − 1) = {max_comps}.",
)

if st.button("Run PCA", type="primary"):
    with st.spinner("Computing PCA..."):
        adata = reduction.run_pca(adata, n_comps=n_comps)
        state.set_adata(adata)
    st.session_state["_reduction_result"] = (
        f"PCA complete — {n_comps} components computed. "
        "Any previous UMAP has been cleared."
        if not state.umap_done()
        else f"PCA complete — {n_comps} components computed."
    )
    st.rerun()

if state.pca_done():
    adata = state.get_adata()
    n_computed = adata.obsm["X_pca"].shape[1]
    st.success(f"PCA computed — {n_computed} components.")
    st.plotly_chart(plots.pca_elbow(adata), use_container_width=True)
    st.caption(
        "The elbow plot shows how much variance each PC explains. "
        "Look for where the curve flattens — that is usually a good cutoff for **n_pcs** below."
    )

# =========================================================================
# Step 2 — UMAP
# =========================================================================
st.divider()
st.subheader("Step 2 — UMAP")

if not state.pca_done():
    st.info("Run PCA first.")
else:
    adata = state.get_adata()
    n_pcs_available = adata.obsm["X_pca"].shape[1]

    st.write(
        "Projects cells into 2D using a neighborhood graph built on PCA coordinates. "
        "The key parameters are explained below."
    )

    col1, col2, col3 = st.columns(3)

    with col1:
        n_pcs = st.slider(
            "PCs used for neighborhood graph",
            min_value=5,
            max_value=n_pcs_available,
            value=min(40, n_pcs_available),
            step=5,
            help=(
                "Informed by the elbow plot above — choose the PC where the curve flattens. "
                "Using too many PCs adds noise; too few loses biological signal."
            ),
        )

    with col2:
        n_neighbors = st.slider(
            "Number of neighbors",
            min_value=5,
            max_value=100,
            value=15,
            step=1,
            help=(
                "How many nearest neighbors each cell connects to. "
                "Low (5–10) → preserves fine local structure. "
                "High (30–50) → emphasizes global relationships."
            ),
        )

    with col3:
        min_dist = st.slider(
            "Min distance",
            min_value=0.0,
            max_value=1.0,
            value=0.1,
            step=0.05,
            help=(
                "Minimum distance between points in 2D UMAP space. "
                "Low values (0.0–0.1) → tight, compact clusters. "
                "High values (0.5–1.0) → more evenly spread layout."
            ),
        )

    label = "Re-run UMAP" if state.umap_done() else "Run UMAP"
    if st.button(label, type="primary"):
        with st.spinner("Building neighborhood graph and computing UMAP..."):
            adata = reduction.run_umap(
                adata,
                n_neighbors=n_neighbors,
                n_pcs=n_pcs,
                min_dist=min_dist,
            )
            state.set_adata(adata)
        st.session_state["_reduction_result"] = (
            f"UMAP complete (n_neighbors={n_neighbors}, "
            f"n_pcs={n_pcs}, min_dist={min_dist})."
        )
        st.rerun()

    if state.umap_done():
        adata = state.get_adata()

        # Build color options: numeric and low-cardinality categorical obs columns
        obs_cols: list[str] = []
        for col in adata.obs.columns:
            if pd.api.types.is_numeric_dtype(adata.obs[col]):
                obs_cols.append(col)
            elif adata.obs[col].nunique() <= 50:
                obs_cols.append(col)

        color_options = ["(none)"] + obs_cols
        color_by_choice = st.selectbox(
            "Color cells by",
            options=color_options,
            help="Select any cell metadata column to use as the color dimension.",
        )
        color_by = None if color_by_choice == "(none)" else color_by_choice

        st.plotly_chart(plots.umap_scatter(adata, color_by=color_by), use_container_width=True)
