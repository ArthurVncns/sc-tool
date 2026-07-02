import anndata as ad
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def doublet_histogram(adata: ad.AnnData) -> go.Figure:
    """Histogram of Scrublet doublet scores, split by predicted class."""
    scores = adata.obs["doublet_score"].values
    is_doublet = adata.obs["predicted_doublet"].values.astype(bool)

    threshold = float(scores[is_doublet].min()) if is_doublet.any() else 1.0

    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=scores[~is_doublet], name="Singlets",
        marker_color="steelblue", opacity=0.75, nbinsx=50,
    ))
    fig.add_trace(go.Histogram(
        x=scores[is_doublet], name="Predicted doublets",
        marker_color="coral", opacity=0.75, nbinsx=50,
    ))
    fig.add_vline(
        x=threshold, line_dash="dash", line_color="red", opacity=0.6,
        annotation_text=f"threshold = {threshold:.2f}",
        annotation_position="top right",
    )
    fig.update_layout(
        barmode="overlay",
        title="Doublet Score Distribution",
        xaxis_title="Doublet score",
        yaxis_title="Number of cells",
        height=350,
        margin=dict(l=40, r=20, t=50, b=40),
        legend=dict(itemsizing="constant"),
    )
    return fig


def qc_violin(values: np.ndarray, title: str, yaxis_label: str) -> go.Figure:
    """Return a Plotly violin plot for a single QC metric distribution.

    Includes a box plot and marks outlier points so extreme values
    are visible without cluttering the main distribution.
    """
    fig = go.Figure(
        go.Violin(
            y=values,
            box_visible=True,
            meanline_visible=True,
            points="outliers",
            fillcolor="lightsteelblue",
            line_color="steelblue",
        )
    )
    fig.update_layout(
        title=title,
        yaxis_title=yaxis_label,
        showlegend=False,
        height=350,
        margin=dict(l=40, r=20, t=50, b=40),
    )
    return fig


def hvg_scatter(adata: ad.AnnData) -> go.Figure:
    """Return a Plotly scatter of mean expression vs. variability metric.

    Handles all three HVG methods:
    - seurat:             y = dispersions_norm   (Workflow 1, log-normalized)
    - seurat_v3:          y = variances_norm     (Workflow 1, LOESS on raw counts)
    - pearson_residuals:  y = residual_variances (Workflow 2)

    Each point is a gene, colored by selection status. Gene name shown on hover.
    """
    # Detect which variability column was produced by the chosen flavor/workflow
    if "dispersions_norm" in adata.var.columns:
        y_col = "dispersions_norm"
        y_label = "Normalized dispersion"
        x_label = "Mean expression (log-normalized)"
    elif "variances_norm" in adata.var.columns:
        y_col = "variances_norm"
        y_label = "Normalized variance (LOESS)"
        x_label = "Mean expression (raw counts)"
    elif "residual_variances" in adata.var.columns:
        y_col = "residual_variances"
        y_label = "Residual variance"
        x_label = "Mean expression (raw counts)"
    else:
        raise ValueError("No variability column found in adata.var. Run select_hvg() or select_hvg_pearson() first.")

    var = adata.var[["means", y_col, "highly_variable"]].copy()
    var["gene"] = var.index
    n_selected = int(var["highly_variable"].sum())

    fig = px.scatter(
        var,
        x="means",
        y=y_col,
        color="highly_variable",
        hover_name="gene",
        color_discrete_map={True: "steelblue", False: "lightgray"},
        labels={
            "means": x_label,
            y_col: y_label,
            "highly_variable": "Highly variable",
        },
        title=f"Highly Variable Genes — {n_selected:,} selected out of {len(var):,}",
    )
    fig.update_traces(marker=dict(size=3, opacity=0.6))
    fig.update_layout(height=420, margin=dict(l=40, r=20, t=50, b=40))
    return fig


def pca_elbow(adata: ad.AnnData) -> go.Figure:
    """Return a dual-axis elbow plot: per-PC and cumulative variance explained.

    The left axis shows individual variance per PC (use this to spot the elbow).
    The right axis shows cumulative variance (use this to pick an n_pcs threshold
    that captures a target percentage, e.g. 80%).
    """
    var_ratio = adata.uns["pca"]["variance_ratio"] * 100
    n = len(var_ratio)
    pcs = list(range(1, n + 1))
    cumulative = np.cumsum(var_ratio).tolist()

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=pcs,
        y=var_ratio.tolist(),
        mode="lines+markers",
        name="Per-PC variance (%)",
        marker=dict(size=5, color="steelblue"),
        line=dict(color="steelblue"),
    ))
    fig.add_trace(go.Scatter(
        x=pcs,
        y=cumulative,
        mode="lines",
        name="Cumulative (%)",
        line=dict(color="coral", dash="dash"),
        yaxis="y2",
    ))
    fig.update_layout(
        title="PCA — Variance Explained",
        xaxis_title="Principal component",
        yaxis=dict(title="Variance explained (%)", rangemode="tozero"),
        yaxis2=dict(
            title="Cumulative variance (%)",
            overlaying="y",
            side="right",
            range=[0, 105],
        ),
        legend=dict(x=0.6, y=0.95),
        height=380,
        margin=dict(l=40, r=60, t=50, b=40),
    )
    return fig


def umap_scatter(adata: ad.AnnData, color_by: str | None = None) -> go.Figure:
    """Return an interactive 2D UMAP scatter plot.

    Each point is a cell. Color can be mapped to any numeric or low-cardinality
    categorical column from adata.obs.

    Args:
        adata: AnnData with 'X_umap' in obsm.
        color_by: Column name from adata.obs to use for coloring, or None.
    """
    coords = adata.obsm["X_umap"]
    df = pd.DataFrame({
        "UMAP1": coords[:, 0],
        "UMAP2": coords[:, 1],
    }, index=adata.obs_names)

    title = "UMAP"

    if color_by is not None and color_by in adata.obs.columns:
        df[color_by] = adata.obs[color_by].values
        title = f"UMAP — colored by {color_by}"

        if pd.api.types.is_numeric_dtype(df[color_by]):
            fig = px.scatter(
                df, x="UMAP1", y="UMAP2", color=color_by,
                color_continuous_scale="Viridis",
                title=title,
            )
        else:
            df[color_by] = df[color_by].astype(str)
            fig = px.scatter(
                df, x="UMAP1", y="UMAP2", color=color_by,
                title=title,
            )
    else:
        fig = px.scatter(df, x="UMAP1", y="UMAP2", title=title)

    fig.update_traces(marker=dict(size=3, opacity=0.7))
    fig.update_layout(
        height=550,
        margin=dict(l=40, r=20, t=50, b=40),
        xaxis=dict(title="UMAP1", showgrid=False, zeroline=False),
        yaxis=dict(title="UMAP2", showgrid=False, zeroline=False),
        # 'constant' decouples legend symbol size from the 3px scatter marker size,
        # keeping legend dots at a readable fixed size.
        legend=dict(itemsizing="constant"),
    )
    return fig


def volcano_plot(
    df: pd.DataFrame,
    group: str,
    reference: str,
    lfc_threshold: float = 0.5,
    pval_threshold: float = 0.05,
) -> go.Figure:
    """Return an interactive volcano plot for differential expression results.

    Args:
        df: DataFrame from de.get_de_results() — columns include
            names, logfoldchanges, pvals_adj.
        group: Name of the test group (shown in legend and title).
        reference: Name of the reference group.
        lfc_threshold: Absolute log fold change cutoff for significance colouring.
        pval_threshold: Adjusted p-value cutoff for significance colouring.
    """
    df = df.copy()
    df["neg_log10_pval"] = -np.log10(df["pvals_adj"].clip(lower=1e-300))

    up = (df["pvals_adj"] < pval_threshold) & (df["logfoldchanges"] > lfc_threshold)
    down = (df["pvals_adj"] < pval_threshold) & (df["logfoldchanges"] < -lfc_threshold)

    df["category"] = "Not significant"
    df.loc[up, "category"] = f"Up in {group}"
    df.loc[down, "category"] = f"Up in {reference}"

    color_map = {
        "Not significant": "lightgray",
        f"Up in {group}": "#6D28D9",
        f"Up in {reference}": "#F97316",
    }

    fig = px.scatter(
        df,
        x="logfoldchanges",
        y="neg_log10_pval",
        color="category",
        hover_name="names",
        color_discrete_map=color_map,
        labels={
            "logfoldchanges": "Log fold change",
            "neg_log10_pval": "−log₁₀(adj. p-value)",
            "category": "",
        },
        title=f"Differential Expression — {group} vs {reference}",
    )

    # Significance threshold lines
    fig.add_hline(
        y=-np.log10(pval_threshold),
        line_dash="dash", line_color="red", opacity=0.4,
    )
    fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="#6D28D9", opacity=0.3)
    fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="#F97316", opacity=0.3)

    fig.update_traces(marker=dict(size=4, opacity=0.7))
    fig.update_layout(
        height=500,
        margin=dict(l=40, r=20, t=60, b=40),
        legend=dict(itemsizing="constant"),
    )
    return fig
