"""Two coherent preprocessing workflows for scRNA-seq analysis.

Workflow 1 — Standard Scanpy:
    normalize_cpm()  →  log_transform()  →  select_hvg(flavor="seurat"|"seurat_v3")

Workflow 2 — Pearson Residuals:
    normalize_pearson()  →  select_hvg_pearson()

Both workflows save raw integer counts to adata.layers['counts'] and record
the active workflow in adata.uns['sc_tool_workflow'].
"""

import anndata as ad
import scanpy as sc


# ── Workflow 1: Standard Scanpy ───────────────────────────────────────────────

def normalize_cpm(adata: ad.AnnData, target_sum: float = 1e4) -> ad.AnnData:
    """Normalize each cell's counts to a common target sum (Workflow 1).

    Saves raw integer counts to adata.layers['counts'] before modifying adata.X.

    Args:
        adata: AnnData after QC filtering.
        target_sum: Target counts per cell. Default 1e4 (counts per 10k).

    Returns:
        adata modified in place.
    """
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    adata.uns["sc_tool_workflow"] = "standard"
    return adata


def log_transform(adata: ad.AnnData) -> ad.AnnData:
    """Apply log(x + 1) transformation (Workflow 1, after normalize_cpm).

    Sets adata.raw to preserve the full-gene log-normalized matrix for
    differential expression, before downstream HVG subsetting.

    Args:
        adata: AnnData after CPM normalization.

    Returns:
        adata modified in place.

    Raises:
        ValueError: If called during the Pearson Residuals workflow.
    """
    if adata.uns.get("sc_tool_workflow") == "pearson":
        raise ValueError(
            "log_transform() is not part of the Pearson Residuals workflow. "
            "Proceed directly to select_hvg_pearson()."
        )
    sc.pp.log1p(adata)
    adata.raw = adata
    return adata


def select_hvg(
    adata: ad.AnnData,
    n_top_genes: int = 2_000,
    flavor: str = "seurat",
) -> ad.AnnData:
    """Select highly variable genes for Workflow 1.

    Args:
        adata: AnnData after log-normalization.
        n_top_genes: Number of HVGs to select. Clamped to adata.n_vars.
        flavor: HVG algorithm.
            - 'seurat': dispersion-based, operates on log-normalized adata.X.
              Adds 'dispersions_norm' to adata.var.
            - 'seurat_v3': LOESS regression on raw count mean-variance
              relationship via adata.layers['counts']. Adds 'variances_norm'.

    Returns:
        adata with adata.var['highly_variable'] added.

    Raises:
        ValueError: If flavor='seurat_v3' and raw counts are unavailable.
    """
    n_top = min(n_top_genes, adata.n_vars)

    if flavor == "seurat_v3":
        if "counts" not in adata.layers:
            raise ValueError(
                "flavor='seurat_v3' requires raw counts in adata.layers['counts']. "
                "Call normalize_cpm() first."
            )
        sc.pp.highly_variable_genes(
            adata, n_top_genes=n_top, flavor="seurat_v3", layer="counts"
        )
    else:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat")

    return adata


# ── Workflow 2: Pearson Residuals ─────────────────────────────────────────────

def normalize_pearson(adata: ad.AnnData) -> ad.AnnData:
    """Pearson residuals normalization (Workflow 2).

    Fits a regularized negative binomial model per gene and replaces counts
    with Pearson residuals. Genes expressed in no cells are removed first
    to avoid division by zero (mu = 0 in the residual denominator).

    Sets adata.raw to preserve the full-gene residuals for differential
    expression, before HVG subsetting.

    Args:
        adata: AnnData after QC filtering with raw integer counts in adata.X.

    Returns:
        adata modified in place. Gene count may decrease (zero-count genes removed).
    """
    adata.layers["counts"] = adata.X.copy()
    # Genes with zero counts in every cell cause mu=0 → division by zero
    # in the residual formula (X - mu) / sqrt(mu + mu²/θ).
    sc.pp.filter_genes(adata, min_cells=1)
    sc.experimental.pp.normalize_pearson_residuals(adata)
    adata.raw = adata
    adata.uns["sc_tool_workflow"] = "pearson"
    return adata


def select_hvg_pearson(
    adata: ad.AnnData,
    n_top_genes: int = 2_000,
) -> ad.AnnData:
    """Select HVGs ranked by residual variance (Workflow 2).

    Computes Pearson residuals from raw counts (adata.layers['counts']) and
    ranks genes by their residual variance, identifying genes with unexpectedly
    high biological variation given their mean expression.

    Note: raw counts are used for ranking — not the normalized adata.X —
    because the normalization step has already equalized variance across genes.

    Args:
        adata: AnnData after normalize_pearson().
        n_top_genes: Number of HVGs to select. Clamped to adata.n_vars.

    Returns:
        adata with adata.var['highly_variable'] and 'residual_variances' added.

    Raises:
        ValueError: If called without prior Pearson residuals normalization.
    """
    if adata.uns.get("sc_tool_workflow") != "pearson":
        raise ValueError(
            "select_hvg_pearson() requires Pearson residuals normalization. "
            "Call normalize_pearson() first, or use select_hvg() for Workflow 1."
        )
    n_top = min(n_top_genes, adata.n_vars)
    sc.experimental.pp.highly_variable_genes(
        adata, n_top_genes=n_top, flavor="pearson_residuals", layer="counts"
    )
    return adata


# ── Shared utility ────────────────────────────────────────────────────────────

def reset_preprocessing(adata: ad.AnnData) -> ad.AnnData:
    """Reset all preprocessing so a different workflow can be selected.

    Restores adata.X from adata.layers['counts'] (the raw counts saved before
    normalization) and clears normalization, HVG, PCA, and UMAP results.

    QC metrics in adata.obs are preserved.
    adata.layers['counts'] is preserved so it survives the reset.

    Args:
        adata: AnnData with adata.layers['counts'] set by a previous normalize call.

    Returns:
        adata modified in place.

    Raises:
        ValueError: If raw counts are not available to restore from.
    """
    if "counts" not in adata.layers:
        raise ValueError(
            "Cannot reset preprocessing: adata.layers['counts'] not found. "
            "Raw counts must have been saved by normalize_cpm() or normalize_pearson()."
        )

    adata.X = adata.layers["counts"].copy()

    for key in ("sc_tool_workflow", "log1p", "pearson_residuals_normalization",
                "pca", "neighbors", "umap", "hvg"):
        adata.uns.pop(key, None)

    # HVG columns to remove — QC columns (mt, ribo, hb, n_genes_by_counts, etc.) are kept.
    _hvg_cols = {
        "highly_variable", "highly_variable_rank", "highly_variable_nbatches",
        "means", "dispersions", "dispersions_norm",
        "variances", "variances_norm", "residual_variances",
    }
    for col in list(adata.var.columns):
        if col in _hvg_cols:
            del adata.var[col]

    for key in ("X_pca", "X_umap"):
        adata.obsm.pop(key, None)

    adata.raw = None

    return adata
