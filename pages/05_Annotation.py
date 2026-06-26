import json

import streamlit as st

from analysis import annotation
from ui import plots, state

st.set_page_config(page_title="sc_tool — Annotation", layout="wide")

st.title("Cell Type Annotation")

# --- Guards ---
if not state.has_adata():
    st.info("Please upload a dataset on the **Upload** page first.")
    st.stop()

if not state.umap_done():
    st.warning("Please complete **Dimensionality Reduction** (PCA + UMAP) first.")
    st.stop()

adata = state.get_adata()
st.caption(f"Current dataset: **{adata.n_obs:,} cells** × {adata.n_vars:,} genes")

# --- Feedback message ---
if "_annotation_result" in st.session_state:
    st.success(st.session_state.pop("_annotation_result"))

# =========================================================================
# Step 1 — Clustering
# =========================================================================
st.subheader("Step 1 — Leiden Clustering")
st.write(
    "Groups cells into clusters based on the neighborhood graph built during UMAP. "
    "The resolution parameter controls granularity: higher values produce more, "
    "smaller clusters."
)

resolution = st.slider(
    "Resolution",
    min_value=0.1, max_value=2.0, value=0.5, step=0.05,
    help="Typical range: 0.3–1.0 for most datasets. Start at 0.5 and adjust.",
)

label = "Re-cluster" if state.clustering_done() else "Cluster"
if st.button(label, type="primary"):
    with st.spinner("Clustering..."):
        adata = annotation.run_clustering(state.get_adata(), resolution=resolution)
        state.set_adata(adata)
    n_clusters = adata.obs["leiden"].nunique()
    st.session_state["_annotation_result"] = (
        f"Clustering complete — {n_clusters} clusters (resolution={resolution})."
    )
    st.rerun()

if state.clustering_done():
    adata = state.get_adata()
    n_clusters = adata.obs["leiden"].nunique()
    st.success(f"{n_clusters} clusters found.")
    st.plotly_chart(plots.umap_scatter(adata, color_by="leiden"), use_container_width=True)

# =========================================================================
# Step 2 — Marker Genes
# =========================================================================
st.divider()
st.subheader("Step 2 — Marker Genes")

if not state.clustering_done():
    st.info("Run clustering first.")
else:
    adata = state.get_adata()
    st.write(
        "Identifies genes significantly more expressed in each cluster "
        "vs. all other cells (Wilcoxon rank-sum test). "
        "Use these to identify the cell type of each cluster."
    )

    n_top = st.slider("Top markers to show per cluster", min_value=3, max_value=20, value=5)

    label = "Re-compute markers" if state.markers_computed() else "Find marker genes"
    if st.button(label, type="primary"):
        with st.spinner("Running Wilcoxon test per cluster..."):
            adata = annotation.find_marker_genes(adata, n_genes=25)
            state.set_adata(adata)
        st.session_state["_annotation_result"] = "Marker genes computed."
        st.rerun()

    if state.markers_computed():
        adata = state.get_adata()
        st.dataframe(annotation.get_top_markers_df(adata, n_genes=n_top), use_container_width=True)
        st.caption(
            "Each column is a cluster. Rank 1 = strongest marker. "
            "Use these to fill in the annotation form below."
        )

# =========================================================================
# Step 3 — Cell Type Annotation
# =========================================================================
st.divider()
st.subheader("Step 3 — Cell Type Annotation")

if not state.markers_computed():
    st.info("Find marker genes first.")
    st.stop()

adata = state.get_adata()

tab_manual, tab_celltypist, tab_scoring = st.tabs([
    "Manual",
    "CellTypist — automated",
    "Marker gene scoring",
])

# ─────────────────────────────────────────────────────────────────
# Tab 1: Manual annotation
# ─────────────────────────────────────────────────────────────────
with tab_manual:
    st.write(
        "Assign a cell type name to each cluster based on the marker genes above. "
        "Leave a field blank to mark the cluster as 'Unknown'."
    )

    clusters = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x))
    n_cols = min(4, len(clusters))

    with st.form("annotation_form"):
        cols = st.columns(n_cols)
        annotations: dict[str, str] = {}

        for i, cluster in enumerate(clusters):
            with cols[i % n_cols]:
                existing = ""
                if "cell_type" in adata.obs.columns:
                    val = adata.obs.loc[adata.obs["leiden"] == cluster, "cell_type"].iloc[0]
                    existing = "" if str(val) == "Unknown" else str(val)
                annotations[cluster] = st.text_input(
                    f"Cluster {cluster}", value=existing, placeholder="e.g. T cells",
                )

        if st.form_submit_button("Apply annotations", type="primary"):
            filled = {k: v.strip() if v.strip() else "Unknown" for k, v in annotations.items()}
            adata = annotation.apply_annotations(adata, filled)
            state.set_adata(adata)
            n_typed = sum(1 for v in filled.values() if v != "Unknown")
            st.session_state["_annotation_result"] = (
                f"Manual annotations applied — {n_typed}/{len(clusters)} clusters labeled."
            )
            st.rerun()

    if state.annotation_done():
        st.plotly_chart(
            plots.umap_scatter(state.get_adata(), color_by="cell_type"),
            use_container_width=True,
        )

# ─────────────────────────────────────────────────────────────────
# Tab 2: CellTypist
# ─────────────────────────────────────────────────────────────────
with tab_celltypist:
    st.write(
        "Automated annotation using pre-trained logistic regression models "
        "from the [CellTypist](https://www.celltypist.org/) project (Wellcome Sanger Institute). "
        "Models are downloaded automatically on first use (~100 MB each)."
    )
    st.caption(
        "**Note on Azimuth:** Azimuth (Satija lab) is an excellent reference-based tool "
        "but is only available as an R package or web app "
        "(https://azimuth.hubmapconsortium.org/). "
        "CellTypist is the recommended Python equivalent."
    )

    if state.active_workflow() == "pearson":
        st.warning(
            "CellTypist expects log-normalized data (Standard workflow). "
            "You used Pearson Residuals — results may be less accurate."
        )

    model_label = st.selectbox(
        "Model",
        options=list(annotation.CELLTYPIST_MODELS.keys()),
        help="Choose a model matching your tissue type.",
    )
    model_name = annotation.CELLTYPIST_MODELS[model_label]

    majority_voting = st.checkbox(
        "Majority voting",
        value=True,
        help=(
            "Refines cell-level predictions using cluster-level voting. "
            "Recommended — uses the Leiden clusters from Step 1."
        ),
    )

    label = "Re-run CellTypist" if state.celltypist_done() else "Run CellTypist"
    if st.button(label, type="primary"):
        with st.spinner(f"Running CellTypist with {model_name}..."):
            try:
                adata = annotation.run_celltypist(
                    state.get_adata(),
                    model_name=model_name,
                    majority_voting=majority_voting,
                )
                state.set_adata(adata)
                n_types = adata.obs["celltypist_cell_type"].nunique()
                st.session_state["_annotation_result"] = (
                    f"CellTypist complete — {n_types} cell types predicted using {model_name}."
                )
            except Exception as e:
                st.error(f"CellTypist failed: {e}")
        st.rerun()

    if state.celltypist_done():
        adata = state.get_adata()
        n_types = adata.obs["celltypist_cell_type"].nunique()
        st.success(f"{n_types} cell types predicted.")
        st.plotly_chart(
            plots.umap_scatter(adata, color_by="celltypist_cell_type"),
            use_container_width=True,
        )

# ─────────────────────────────────────────────────────────────────
# Tab 3: Marker gene scoring
# ─────────────────────────────────────────────────────────────────
with tab_scoring:
    st.write(
        "Score each cell against a custom marker gene dictionary. "
        "Each cell is assigned to the cell type with the highest average marker gene expression. "
        "Works for any species or tissue — no pre-trained model required."
    )

    default_json = json.dumps(annotation.PBMC_MARKERS, indent=2)
    marker_input = st.text_area(
        "Marker gene dictionary (JSON)",
        value=default_json,
        height=280,
        help=(
            'Format: {"cell type": ["GENE1", "GENE2", ...], ...}. '
            "Gene names must match your dataset's var_names exactly."
        ),
    )

    label = "Re-run scoring" if state.marker_score_done() else "Run marker gene scoring"
    if st.button(label, type="primary"):
        try:
            marker_dict = json.loads(marker_input)
            if not isinstance(marker_dict, dict):
                raise ValueError("Input must be a JSON object.")
        except (json.JSONDecodeError, ValueError) as e:
            st.error(f"Invalid JSON: {e}")
        else:
            with st.spinner("Scoring cells against marker genes..."):
                try:
                    adata = annotation.score_marker_genes(state.get_adata(), marker_dict)
                    state.set_adata(adata)
                    n_types = adata.obs["marker_score_cell_type"].nunique()
                    st.session_state["_annotation_result"] = (
                        f"Marker scoring complete — {n_types} cell types assigned."
                    )
                except ValueError as e:
                    st.error(str(e))
            st.rerun()

    if state.marker_score_done():
        adata = state.get_adata()
        n_types = adata.obs["marker_score_cell_type"].nunique()
        st.success(f"{n_types} cell types assigned by marker gene scoring.")
        st.plotly_chart(
            plots.umap_scatter(adata, color_by="marker_score_cell_type"),
            use_container_width=True,
        )
