import tempfile
from pathlib import Path

import streamlit as st

from ui import state
from ui.theme import apply_theme


def _build_summary(adata) -> str:
    """Generate a plain-text analysis summary from the current AnnData state."""
    lines = ["# sc_tool — Analysis Summary", ""]

    lines += [
        "## Dataset",
        f"- Cells: {adata.n_obs:,}",
        f"- Genes: {adata.n_vars:,}",
        "",
    ]

    lines.append("## Quality Control")
    if "n_genes_by_counts" in adata.obs.columns:
        lines.append("- QC metrics: computed")
    if "doublet_score" in adata.obs.columns:
        n_d = int(adata.obs["predicted_doublet"].sum())
        lines.append(f"- Doublet detection: {n_d:,} predicted doublets")
    lines.append("")

    workflow = adata.uns.get("sc_tool_workflow")
    if workflow:
        lines.append("## Preprocessing")
        lines.append(
            f"- Workflow: {'Standard (CPM + log1p)' if workflow == 'standard' else 'Pearson Residuals'}"
        )
        if "highly_variable" in adata.var.columns:
            lines.append(f"- Highly variable genes: {int(adata.var['highly_variable'].sum()):,}")
        lines.append("")

    if "X_pca" in adata.obsm:
        lines.append("## Dimensionality Reduction")
        lines.append(f"- PCA: {adata.obsm['X_pca'].shape[1]} components")
        if "X_pca_harmony" in adata.obsm:
            lines.append(f"- Harmony: applied (batch: '{adata.uns.get('harmony_batch_key', '?')}')")
        if "X_umap" in adata.obsm:
            lines.append("- UMAP: computed")
        lines.append("")

    if "leiden" in adata.obs.columns:
        lines.append("## Clustering")
        lines.append(f"- Leiden: {adata.obs['leiden'].nunique()} clusters")
        lines.append("")

    ann_cols = [
        ("cell_type", "Manual"),
        ("celltypist_cell_type", "CellTypist"),
        ("marker_score_cell_type", "Marker gene scoring"),
    ]
    present = [(col, label) for col, label in ann_cols if col in adata.obs.columns]
    if present:
        lines.append("## Cell Type Annotation")
        for col, label in present:
            lines.append(f"- {label}: {adata.obs[col].nunique()} cell types")
        lines.append("")

    if "sc_tool_de" in adata.uns:
        de = adata.uns["sc_tool_de"]
        lines.append("## Differential Expression")
        lines.append(f"- Comparison: {de['group']} vs {de['reference']}")
        lines.append(f"- Variable: {de['groupby']}")
        lines.append(f"- Method: {de['method']}")
        lines.append("")

    return "\n".join(lines)


def _adata_to_bytes(adata) -> bytes:
    """Write adata to a temp h5ad file and return the bytes."""
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        adata.write_h5ad(tmp_path)
        return tmp_path.read_bytes()
    finally:
        tmp_path.unlink(missing_ok=True)


# ── Page ───────────────────────────────────────────────────────────────────────
st.set_page_config(page_title="sc_tool — Export", layout="wide")
apply_theme()

st.title("Export & Summary")

if not state.has_adata():
    st.info("Please upload a dataset on the **Upload** page first.")
    st.stop()

adata = state.get_adata()
st.caption(f"Current dataset: **{adata.n_obs:,} cells** × {adata.n_vars:,} genes")

# ── Analysis summary ───────────────────────────────────────────────────────────
st.subheader("Analysis Summary")
summary = _build_summary(adata)
st.markdown(summary)

st.download_button(
    "Download summary (.md)",
    data=summary,
    file_name="sc_tool_summary.md",
    mime="text/markdown",
    type="secondary",
)

# ── Data downloads ─────────────────────────────────────────────────────────────
st.divider()
st.subheader("Download Data")

col1, col2 = st.columns(2)

with col1:
    st.markdown("**Processed AnnData (.h5ad)**")
    st.caption("Full dataset with all computed embeddings and metadata.")
    if st.button("Prepare .h5ad", type="primary"):
        with st.spinner("Writing h5ad file…"):
            st.session_state["_h5ad_bytes"] = _adata_to_bytes(adata)
        st.rerun()
    if "_h5ad_bytes" in st.session_state:
        st.download_button(
            "Download .h5ad",
            data=st.session_state["_h5ad_bytes"],
            file_name="sc_tool_processed.h5ad",
            mime="application/octet-stream",
        )

with col2:
    st.markdown("**Cell metadata (obs) — CSV**")
    st.caption(f"{len(adata.obs.columns)} columns × {adata.n_obs:,} cells.")
    st.download_button(
        "Download cell metadata",
        data=adata.obs.to_csv(),
        file_name="cell_metadata.csv",
        mime="text/csv",
    )

col3, col4 = st.columns(2)

with col3:
    st.markdown("**Gene metadata (var) — CSV**")
    st.caption(f"{len(adata.var.columns)} columns × {adata.n_vars:,} genes.")
    st.download_button(
        "Download gene metadata",
        data=adata.var.to_csv(),
        file_name="gene_metadata.csv",
        mime="text/csv",
    )

with col4:
    if state.de_done():
        from analysis.de import get_de_results
        de_info = adata.uns["sc_tool_de"]
        st.markdown("**DE results — CSV**")
        st.caption(
            f"{de_info['group']} vs {de_info['reference']} "
            f"({de_info['method']})."
        )
        st.download_button(
            "Download DE results",
            data=get_de_results(adata).to_csv(index=False),
            file_name=f"de_{de_info['group']}_vs_{de_info['reference']}.csv",
            mime="text/csv",
        )
    else:
        st.markdown("**DE results — CSV**")
        st.caption("Run differential expression analysis to enable this download.")
