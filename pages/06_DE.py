import numpy as np
import pandas as pd
import streamlit as st

from analysis import de
from ui import plots, state
from ui.theme import apply_theme, next_page_button

st.set_page_config(page_title="sc_tool — Differential Expression", layout="wide")
apply_theme()

st.title("Differential Expression")

# --- Guard ---
if not state.has_adata():
    st.info("Please upload a dataset on the **Upload** page first.")
    st.stop()

adata = state.get_adata()
st.caption(f"Current dataset: **{adata.n_obs:,} cells** × {adata.n_vars:,} genes")

if not state.is_log_transformed():
    st.warning(
        "Data does not appear to be log-normalized. "
        "DE analysis is most meaningful on log-normalized data — "
        "complete **Preprocessing** first for best results."
    )

# --- Feedback message ---
if "_de_result" in st.session_state:
    st.success(st.session_state.pop("_de_result"))

# =========================================================================
# Step 1 — Define groups
# =========================================================================
st.subheader("Step 1 — Define Groups")

source = st.radio(
    "Group source",
    options=["Existing metadata variable", "Load custom metadata file (CSV / TSV)"],
    horizontal=True,
    key="de_source",
)

# ── Custom metadata upload ─────────────────────────────────────────────────
if source == "Load custom metadata file (CSV / TSV)":
    st.write(
        "Upload a CSV or TSV file with one row per cell. "
        "One column must match the cell barcodes in your dataset; "
        "another column should contain the group labels."
    )

    uploaded_meta = st.file_uploader(
        "Metadata file", type=["csv", "tsv", "txt"], key="de_meta_file"
    )

    if uploaded_meta is not None:
        try:
            meta_df = de.load_metadata_file(uploaded_meta)
        except Exception as exc:
            st.error(f"Could not parse file: {exc}")
            st.stop()

        st.dataframe(meta_df.head(5), use_container_width=True)

        col1, col2, col3 = st.columns(3)
        with col1:
            id_col = st.selectbox("Cell barcode column", meta_df.columns.tolist())
        with col2:
            label_col = st.selectbox(
                "Group label column",
                [c for c in meta_df.columns if c != id_col],
            )
        with col3:
            obs_key = st.text_input("Add to obs as", value="custom_group")

        if st.button("Add to cell metadata", type="primary"):
            n_before = adata.obs[obs_key].notna().sum() if obs_key in adata.obs else 0
            adata = de.add_obs_column(adata, meta_df, id_col, label_col, obs_key)
            state.set_adata(adata)
            n_matched = adata.obs[obs_key].notna().sum()
            st.session_state["_de_result"] = (
                f"Added '{obs_key}' to cell metadata — "
                f"{n_matched:,} / {adata.n_obs:,} cells matched."
            )
            st.rerun()
    else:
        st.stop()

# ── Group selector ─────────────────────────────────────────────────────────
# Offer all obs columns that have between 2 and 50 unique non-null values.
adata = state.get_adata()
groupable = [
    col for col in adata.obs.columns
    if 2 <= adata.obs[col].dropna().nunique() <= 50
]

if not groupable:
    st.warning(
        "No suitable grouping variable found in cell metadata. "
        "Upload a custom metadata file or run clustering / annotation first."
    )
    st.stop()

st.divider()
groupby = st.selectbox("Metadata variable", groupable, key="de_groupby")

unique_values = sorted(adata.obs[groupby].dropna().astype(str).unique().tolist())
n_per_group = adata.obs[groupby].value_counts()

col1, col2 = st.columns(2)
with col1:
    group = st.selectbox(
        "Group A (test)",
        unique_values,
        help="The group you want to find marker genes for.",
    )
with col2:
    reference_options = ["Rest of cells"] + [v for v in unique_values if v != group]
    reference_choice = st.selectbox(
        "Reference (Group B)",
        reference_options,
        help="'Rest of cells' compares Group A against all other cells.",
    )
    reference = "rest" if reference_choice == "Rest of cells" else reference_choice

# Show group sizes
n_group = int((adata.obs[groupby].astype(str) == group).sum())
n_ref = adata.n_obs - n_group if reference == "rest" else int((adata.obs[groupby].astype(str) == reference).sum())
st.caption(f"Group A: **{n_group:,} cells** — Reference: **{n_ref:,} cells**")

# =========================================================================
# Step 2 — Test settings
# =========================================================================
st.divider()
st.subheader("Step 2 — Statistical Test")

col1, col2 = st.columns([2, 1])
with col1:
    method = st.radio(
        "Method",
        options=["wilcoxon", "t-test"],
        horizontal=True,
        help=(
            "**Wilcoxon** (recommended): non-parametric rank-sum test. "
            "Makes no assumptions about gene expression distributions — "
            "robust and appropriate for most scRNA-seq datasets. "
            "Produces p-values and fold changes for the volcano plot.\n\n"
            "**t-test**: faster parametric test, assumes approximately normal "
            "distributions. Can overestimate significance on zero-inflated "
            "scRNA-seq data. Good for quick explorations on large datasets."
        ),
    )
with col2:
    n_genes = st.slider("Genes to rank", min_value=50, max_value=500, value=200, step=50)

label = "Re-run DE" if state.de_done() else "Run Differential Expression"
if st.button(label, type="primary"):
    with st.spinner(f"Running {method} test — {group} vs {reference}…"):
        adata = de.run_de(
            state.get_adata(),
            groupby=groupby,
            group=group,
            reference=reference,
            method=method,
            n_genes=n_genes,
        )
        state.set_adata(adata)
    results = de.get_de_results(adata)
    n_sig = int(((results["pvals_adj"] < 0.05) & (results["logfoldchanges"].abs() > 0.5)).sum())
    st.session_state["_de_result"] = (
        f"DE complete — **{n_sig}** significant genes "
        f"(adj. p < 0.05, |log FC| > 0.5) out of {len(results):,} tested."
    )
    st.rerun()

# =========================================================================
# Step 3 — Results
# =========================================================================
if state.de_done():
    adata = state.get_adata()
    de_info = adata.uns["sc_tool_de"]

    # Only show results if they match the current comparison
    if de_info["groupby"] != groupby or de_info["group"] != group:
        st.info(
            f"Showing cached results from a previous run "
            f"({de_info['group']} vs {de_info['reference']} "
            f"on '{de_info['groupby']}'). "
            "Re-run to update."
        )

    st.divider()
    st.subheader("Step 3 — Results")

    results = de.get_de_results(adata)
    sig_mask = (results["pvals_adj"] < 0.05) & (results["logfoldchanges"].abs() > 0.5)
    n_up = int(((results["pvals_adj"] < 0.05) & (results["logfoldchanges"] > 0.5)).sum())
    n_down = int(((results["pvals_adj"] < 0.05) & (results["logfoldchanges"] < -0.5)).sum())

    c1, c2, c3 = st.columns(3)
    c1.metric(f"Up in {de_info['group']}", n_up)
    c2.metric(f"Up in {de_info['reference']}", n_down)
    c3.metric("Not significant", len(results) - n_up - n_down)

    # Volcano plot
    st.plotly_chart(
        plots.volcano_plot(results, de_info["group"], de_info["reference"]),
        use_container_width=True,
    )

    # Significant gene table
    st.subheader("Significant genes")
    pval_cutoff = st.slider(
        "Adjusted p-value cutoff",
        min_value=0.001, max_value=0.1, value=0.05, step=0.001, format="%.3f",
    )
    lfc_cutoff = st.slider(
        "Min |log fold change|",
        min_value=0.0, max_value=2.0, value=0.5, step=0.1,
    )

    sig = results[
        (results["pvals_adj"] < pval_cutoff) & (results["logfoldchanges"].abs() > lfc_cutoff)
    ].copy()
    sig = sig.sort_values("pvals_adj").reset_index(drop=True)
    sig["logfoldchanges"] = sig["logfoldchanges"].round(3)
    sig["pvals_adj"] = sig["pvals_adj"].apply(lambda x: f"{x:.2e}")

    st.dataframe(
        sig.rename(columns={
            "names": "Gene",
            "logfoldchanges": "Log FC",
            "pvals_adj": "Adj. p-value",
            "scores": "Score",
        })[["Gene", "Log FC", "Adj. p-value", "Score"]],
        use_container_width=True,
        hide_index=True,
    )

    # Download
    csv = results.to_csv(index=False)
    st.download_button(
        "Download full results (CSV)",
        data=csv,
        file_name=f"de_{de_info['group']}_vs_{de_info['reference']}.csv",
        mime="text/csv",
    )

    next_page_button("Export", "pages/07_Export.py")
