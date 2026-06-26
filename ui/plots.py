import anndata as ad
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


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
    )
    return fig
