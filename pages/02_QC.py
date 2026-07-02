import streamlit as st

from analysis import doublets, qc
from ui import plots, state
from ui.theme import apply_theme, next_page_button

st.set_page_config(page_title="sc_tool — QC", layout="wide")
apply_theme()

st.title("Quality Control")

# --- Guard: require an uploaded dataset ---
if not state.has_adata():
    st.info("Please upload a dataset on the **Upload** page first.")
    st.stop()

adata = state.get_adata()

# --- Auto-compute QC metrics on first visit ---
if not state.qc_is_computed():
    with st.spinner("Computing QC metrics..."):
        adata = qc.compute_qc_metrics(adata)
        state.set_adata(adata)

    detected = [
        name
        for name, flag in [
            ("mitochondrial", state.mt_detected()),
            ("ribosomal", state.ribo_detected()),
            ("hemoglobin", state.hb_detected()),
        ]
        if flag
    ]
    if not detected:
        st.info(
            "No mitochondrial, ribosomal, or hemoglobin genes detected. "
            "Expected prefixes: MT-/mt- (mt), RPS/RPL (ribo), HB (hb)."
        )

# Always read fresh from session_state so we work with the latest version.
adata = state.get_adata()

# --- Filter result message (set in session_state before st.rerun()) ---
if "_filter_result" in st.session_state:
    st.success(st.session_state.pop("_filter_result"))

# --- Violin plots: row 1 (always shown) ---
st.subheader("QC Metric Distributions")

col1, col2 = st.columns(2)
with col1:
    st.plotly_chart(
        plots.qc_violin(
            adata.obs["n_genes_by_counts"].values,
            title="Genes per cell",
            yaxis_label="n_genes_by_counts",
        ),
        use_container_width=True,
    )
with col2:
    st.plotly_chart(
        plots.qc_violin(
            adata.obs["total_counts"].values,
            title="Total counts per cell",
            yaxis_label="total_counts",
        ),
        use_container_width=True,
    )

# --- Violin plots: row 2 (gene set percentages, only if detected) ---
gene_set_metrics = [
    ("pct_counts_mt",   "Mitochondrial %"),
    ("pct_counts_ribo", "Ribosomal %"),
    ("pct_counts_hb",   "Hemoglobin %"),
]
available_metrics = [
    (col_name, title)
    for col_name, title in gene_set_metrics
    if col_name in adata.obs.columns
]

if available_metrics:
    pct_cols = st.columns(len(available_metrics))
    for i, (col_name, title) in enumerate(available_metrics):
        with pct_cols[i]:
            st.plotly_chart(
                plots.qc_violin(
                    adata.obs[col_name].values,
                    title=title,
                    yaxis_label=col_name,
                ),
                use_container_width=True,
            )

# --- Doublet Detection (optional) ---
st.divider()
st.subheader("Doublet Detection (optional)")
st.write(
    "Scrublet simulates artificial doublets and scores each real cell by its "
    "similarity to them. Run before filtering so you can remove doublets in one step."
)

doublet_rate = st.slider(
    "Expected doublet rate",
    min_value=0.01, max_value=0.20, value=0.05, step=0.01,
    help="Typical 10x rates: ~4% for 4k cells, ~8% for 8k cells.",
)

label = "Re-run doublet detection" if state.doublet_detected() else "Detect doublets"
if st.button(label, type="secondary"):
    with st.spinner("Running Scrublet…"):
        adata = doublets.run_scrublet(state.get_adata(), expected_doublet_rate=doublet_rate)
        state.set_adata(adata)
    n_pred = int(adata.obs["predicted_doublet"].sum())
    st.session_state["_filter_result"] = (
        f"Doublet detection complete — **{n_pred:,}** predicted doublets "
        f"({n_pred / adata.n_obs * 100:.1f}% of cells)."
    )
    st.rerun()

if state.doublet_detected():
    adata = state.get_adata()
    n_pred = int(adata.obs["predicted_doublet"].sum())
    st.success(f"{n_pred:,} predicted doublets ({n_pred / adata.n_obs * 100:.1f}%).")
    st.plotly_chart(plots.doublet_histogram(adata), use_container_width=True)

# --- Filtering ---
st.divider()
st.subheader("Cell Filtering")
st.caption(f"Current dataset: **{adata.n_obs:,} cells** × {adata.n_vars:,} genes")

max_genes_in_data = int(adata.obs["n_genes_by_counts"].max())
max_counts_in_data = int(adata.obs["total_counts"].max())

col1, col2 = st.columns(2)

with col1:
    min_genes = st.slider(
        "Min genes per cell",
        min_value=0,
        max_value=max_genes_in_data,
        value=200,
        step=50,
        help="Remove cells with too few detected genes (likely empty droplets).",
    )
    max_genes = st.slider(
        "Max genes per cell",
        min_value=0,
        max_value=max_genes_in_data,
        value=min(5_000, max_genes_in_data),
        step=50,
        help="Remove cells with too many detected genes (likely doublets).",
    )

with col2:
    min_counts = st.slider(
        "Min counts per cell",
        min_value=0,
        max_value=max_counts_in_data,
        value=500,
        step=100,
        help="Remove cells with very few total UMI counts.",
    )
    max_counts = st.slider(
        "Max counts per cell",
        min_value=0,
        max_value=max_counts_in_data,
        value=min(25_000, max_counts_in_data),
        step=100,
        help="Remove cells with unusually high counts (likely doublets).",
    )

max_pct_mt: float | None = None
max_pct_ribo: float | None = None
max_pct_hb: float | None = None

if state.mt_detected() or state.ribo_detected() or state.hb_detected():
    st.markdown("**Gene set filters**")
    pct_slider_cols = st.columns(3)

    with pct_slider_cols[0]:
        if state.mt_detected():
            max_pct_mt = st.slider(
                "Max mitochondrial %",
                min_value=0.0, max_value=100.0, value=20.0, step=0.5,
                help="Remove cells with high mitochondrial content (likely damaged).",
            )

    with pct_slider_cols[1]:
        if state.ribo_detected():
            max_pct_ribo = st.slider(
                "Max ribosomal %",
                min_value=0.0, max_value=100.0, value=50.0, step=0.5,
                help="Remove cells with unusually high ribosomal gene expression.",
            )

    with pct_slider_cols[2]:
        if state.hb_detected():
            max_pct_hb = st.slider(
                "Max hemoglobin %",
                min_value=0.0, max_value=100.0, value=1.0, step=0.1,
                help="Remove cells with hemoglobin expression (red blood cell contamination).",
            )

# --- Doublet removal toggle ---
remove_doublets = False
if state.doublet_detected():
    n_pred = int(adata.obs["predicted_doublet"].sum())
    remove_doublets = st.checkbox(
        f"Remove **{n_pred:,}** predicted doublets",
        value=True,
        help="Cells flagged by Scrublet will be excluded in addition to the filters above.",
    )

# --- Live preview ---
filters = qc.QCFilters(
    min_genes=min_genes,
    max_genes=max_genes,
    min_counts=min_counts,
    max_counts=max_counts,
    max_pct_mt=max_pct_mt,
    max_pct_ribo=max_pct_ribo,
    max_pct_hb=max_pct_hb,
)

n_passing = qc.count_cells_passing_filters(adata, filters, remove_doublets=remove_doublets)
n_removed = adata.n_obs - n_passing
pct_removed = n_removed / adata.n_obs * 100

st.markdown(
    f"After filtering: **{n_passing:,} cells** kept "
    f"({n_removed:,} removed — {pct_removed:.1f}%)"
)

# --- Apply ---
if st.button("Apply Filters", type="primary", disabled=(n_passing == 0)):
    filtered = qc.filter_cells(adata, filters, remove_doublets=remove_doublets)
    state.set_adata(filtered)
    st.session_state["_filter_result"] = (
        f"Dataset filtered. Removed **{n_removed:,} cells** ({pct_removed:.1f}%). "
        f"**{filtered.n_obs:,} cells** remaining."
    )
    st.rerun()

if state.qc_is_computed():
    next_page_button("Preprocessing", "pages/03_Preprocessing.py")
