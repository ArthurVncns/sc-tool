import streamlit as st

from analysis import preprocessing
from ui import plots, state

st.set_page_config(page_title="sc_tool — Preprocessing", layout="wide")

st.title("Preprocessing")

# --- Guards ---
if not state.has_adata():
    st.info("Please upload a dataset on the **Upload** page first.")
    st.stop()

if not state.qc_is_computed():
    st.warning("Please complete **QC** before preprocessing.")
    st.stop()

adata = state.get_adata()
st.caption(f"Current dataset: **{adata.n_obs:,} cells** × {adata.n_vars:,} genes")

# --- Feedback message ---
if "_preprocessing_result" in st.session_state:
    st.success(st.session_state.pop("_preprocessing_result"))

# =========================================================================
# Workflow selection
# =========================================================================
workflow = state.active_workflow()

if workflow is None:
    # No workflow started yet — present the choice
    st.subheader("Choose a preprocessing workflow")
    st.write(
        "Select one of the two workflows below. "
        "The choice is locked once normalization starts."
    )

    col1, col2 = st.columns(2, gap="large")

    with col1:
        st.markdown("#### Workflow 1 — Standard Scanpy")
        st.markdown(
            "Normalize counts per cell → log1p → select highly variable genes.\n\n"
            "**HVG options:** Seurat (dispersion) or LOESS (Seurat v3, raw counts).\n\n"
            "Best for most standard analyses."
        )

    with col2:
        st.markdown("#### Workflow 2 — Pearson Residuals")
        st.markdown(
            "Fit a negative binomial model per gene, replace counts with Pearson residuals → "
            "select HVGs ranked by residual variance.\n\n"
            "**HVG:** ranked by residual variance (no flavor choice needed).\n\n"
            "Better for datasets with large depth differences."
        )

    workflow_choice = st.radio(
        "Active workflow",
        ["Workflow 1 — Standard Scanpy", "Workflow 2 — Pearson Residuals"],
        label_visibility="collapsed",
    )
    chosen = "standard" if "Standard" in workflow_choice else "pearson"

    st.divider()

    if chosen == "standard":
        st.subheader("Step 1 — Normalize (CPM)")
        target_sum = st.number_input(
            "Target sum (counts per cell)",
            min_value=100, max_value=1_000_000, value=10_000, step=1_000,
        )
        if st.button("Normalize", type="primary"):
            with st.spinner("Normalizing..."):
                adata = preprocessing.normalize_cpm(adata, target_sum=float(target_sum))
                state.set_adata(adata)
            st.session_state["_preprocessing_result"] = (
                f"CPM normalization complete (target = {int(target_sum):,}). "
                "Raw counts saved to `adata.layers['counts']`."
            )
            st.rerun()
    else:
        st.subheader("Step 1 — Pearson Residuals Normalization")
        st.write(
            "Fits a regularized negative binomial model per gene and replaces counts "
            "with Pearson residuals. Genes with zero counts in all cells are removed first."
        )
        if st.button("Normalize", type="primary"):
            with st.spinner("Computing Pearson residuals..."):
                n_before = adata.n_vars
                adata = preprocessing.normalize_pearson(adata)
                n_removed = n_before - adata.n_vars
                state.set_adata(adata)
            msg = "Pearson residuals normalization complete."
            if n_removed > 0:
                msg += f" {n_removed:,} zero-count genes removed."
            st.session_state["_preprocessing_result"] = msg
            st.rerun()

    st.stop()  # Don't render the pipeline steps until a workflow is started


# =========================================================================
# Workflow 1 — Standard Scanpy
# =========================================================================
def _render_switch_button() -> None:
    """Show a collapsible switch-workflow control with a confirmation step."""
    with st.expander("Switch workflow"):
        st.warning(
            "This will restore raw counts and clear normalization, HVG selection, "
            "PCA, and UMAP. QC results are preserved."
        )
        if st.button("Confirm — switch workflow", type="secondary"):
            adata = preprocessing.reset_preprocessing(state.get_adata())
            state.set_adata(adata)
            st.rerun()


if workflow == "standard":
    st.info("Active workflow: **Standard Scanpy** (CPM + log1p + HVG)")
    _render_switch_button()

    # Step 1 — done
    st.subheader("Step 1 — CPM normalization")
    st.success("Complete. Raw counts preserved in `adata.layers['counts']`.")

    # Step 2 — Log transform
    st.divider()
    st.subheader("Step 2 — Log transformation")

    if state.is_log_transformed():
        st.success("Complete. All values are on a log1p scale.")
    else:
        st.write("Applies log(x + 1) to stabilize variance before HVG selection, PCA, and UMAP.")
        if st.button("Log-transform", type="primary"):
            with st.spinner("Applying log1p..."):
                adata = preprocessing.log_transform(state.get_adata())
                state.set_adata(adata)
            st.session_state["_preprocessing_result"] = "Log transformation complete."
            st.rerun()

    # Step 3 — HVG
    if state.is_log_transformed():
        st.divider()
        st.subheader("Step 3 — Highly Variable Genes")
        adata = state.get_adata()

        if state.hvg_selected():
            n_hvg = int(adata.var["highly_variable"].sum())
            st.success(f"**{n_hvg:,}** highly variable genes selected.")
            st.plotly_chart(plots.hvg_scatter(adata), use_container_width=True)
            st.caption("Adjust settings below and re-select if needed.")

        flavor_choice = st.radio(
            "Selection method",
            ["Seurat (dispersion on log-normalized data)",
             "LOESS / Seurat v3 (variance on raw counts)"],
            help=(
                "**Seurat**: uses normalized dispersion. Fast and widely used.\n\n"
                "**LOESS / Seurat v3**: fits a LOESS curve to the mean-variance "
                "relationship of raw counts. More accurate for datasets with "
                "variable sequencing depth."
            ),
        )
        chosen_flavor = "seurat_v3" if "LOESS" in flavor_choice else "seurat"

        n_top_genes = st.slider(
            "Number of highly variable genes",
            min_value=500, max_value=min(6_000, adata.n_vars),
            value=min(2_000, adata.n_vars), step=100,
        )

        label = "Re-select HVGs" if state.hvg_selected() else "Select HVGs"
        if st.button(label, type="primary"):
            with st.spinner("Selecting highly variable genes..."):
                adata = preprocessing.select_hvg(adata, n_top_genes=n_top_genes, flavor=chosen_flavor)
                state.set_adata(adata)
            method_label = "LOESS (seurat_v3)" if chosen_flavor == "seurat_v3" else "Seurat"
            st.session_state["_preprocessing_result"] = (
                f"{int(adata.var['highly_variable'].sum()):,} HVGs selected ({method_label})."
            )
            st.rerun()


# =========================================================================
# Workflow 2 — Pearson Residuals
# =========================================================================
elif workflow == "pearson":
    st.info("Active workflow: **Pearson Residuals**")
    _render_switch_button()

    # Step 1 — done
    st.subheader("Step 1 — Pearson residuals normalization")
    st.success(
        f"Complete. `adata.X` contains Pearson residuals. "
        f"Raw counts preserved in `adata.layers['counts']`."
    )

    # Step 2 — HVG
    st.divider()
    st.subheader("Step 2 — Highly Variable Genes")
    adata = state.get_adata()

    st.write(
        "Selects genes whose Pearson residuals show the highest variance across cells. "
        "Raw counts are used for ranking, not the normalized residuals — "
        "because normalization has already equalized variance across genes."
    )

    if state.hvg_selected():
        n_hvg = int(adata.var["highly_variable"].sum())
        st.success(f"**{n_hvg:,}** highly variable genes selected by residual variance.")
        st.plotly_chart(plots.hvg_scatter(adata), use_container_width=True)
        st.caption("Adjust settings below and re-select if needed.")

    n_top_genes = st.slider(
        "Number of highly variable genes",
        min_value=500, max_value=min(6_000, adata.n_vars),
        value=min(2_000, adata.n_vars), step=100,
    )

    label = "Re-select HVGs" if state.hvg_selected() else "Select HVGs"
    if st.button(label, type="primary"):
        with st.spinner("Selecting highly variable genes..."):
            adata = preprocessing.select_hvg_pearson(adata, n_top_genes=n_top_genes)
            state.set_adata(adata)
        st.session_state["_preprocessing_result"] = (
            f"{int(adata.var['highly_variable'].sum()):,} HVGs selected by residual variance."
        )
        st.rerun()
