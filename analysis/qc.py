from dataclasses import dataclass

import anndata as ad
import pandas as pd
import scanpy as sc


@dataclass
class QCFilters:
    """Thresholds used to filter low-quality cells from an AnnData object."""

    min_genes: int = 200
    max_genes: int = 5_000
    min_counts: int = 500
    max_counts: int = 25_000
    max_pct_mt: float | None = None    # None means skip this filter
    max_pct_ribo: float | None = None
    max_pct_hb: float | None = None


def compute_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    """Compute per-cell QC metrics and annotate adata.obs in place.

    Always adds n_genes_by_counts and total_counts.
    Conditionally adds percentage columns for gene sets that are detected:
      - pct_counts_mt   (mitochondrial: MT-/mt- prefix)
      - pct_counts_ribo (ribosomal: RPS/RPL prefix)
      - pct_counts_hb   (hemoglobin: HB prefix, excluding HBP/HBEGF)

    Args:
        adata: Input AnnData object.

    Returns:
        adata modified in place with new obs columns.
    """
    qc_vars: list[str] = []

    mt_mask = adata.var_names.str.startswith(("MT-", "mt-"))
    if mt_mask.any():
        adata.var["mt"] = mt_mask
        qc_vars.append("mt")

    ribo_mask = adata.var_names.str.startswith(("RPS", "RPL", "Rps", "Rpl"))
    if ribo_mask.any():
        adata.var["ribo"] = ribo_mask
        qc_vars.append("ribo")

    # HBP encodes HBEGF (heparin-binding EGF), not a hemoglobin gene — exclude it.
    hb_mask = adata.var_names.str.contains("^HB[^P]", regex=True)
    if hb_mask.any():
        adata.var["hb"] = hb_mask
        qc_vars.append("hb")

    # percent_top=None skips "top-N gene" proportions, which require n_vars >= 500.
    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=None, inplace=True)

    return adata


def _build_filter_mask(adata: ad.AnnData, filters: QCFilters) -> pd.Series:
    obs = adata.obs
    mask = (
        (obs["n_genes_by_counts"] >= filters.min_genes)
        & (obs["n_genes_by_counts"] <= filters.max_genes)
        & (obs["total_counts"] >= filters.min_counts)
        & (obs["total_counts"] <= filters.max_counts)
    )
    if filters.max_pct_mt is not None and "pct_counts_mt" in obs.columns:
        mask &= obs["pct_counts_mt"] <= filters.max_pct_mt
    if filters.max_pct_ribo is not None and "pct_counts_ribo" in obs.columns:
        mask &= obs["pct_counts_ribo"] <= filters.max_pct_ribo
    if filters.max_pct_hb is not None and "pct_counts_hb" in obs.columns:
        mask &= obs["pct_counts_hb"] <= filters.max_pct_hb
    return mask


def filter_cells(adata: ad.AnnData, filters: QCFilters) -> ad.AnnData:
    """Return a filtered copy of adata with low-quality cells removed.

    Args:
        adata: AnnData with QC metrics already in adata.obs.
        filters: Threshold values to apply.

    Returns:
        New AnnData containing only cells that pass all filters.
    """
    return adata[_build_filter_mask(adata, filters)].copy()


def count_cells_passing_filters(adata: ad.AnnData, filters: QCFilters) -> int:
    """Return the number of cells that would pass the given filters.

    Cheaper than filter_cells() when only the count is needed.
    """
    return int(_build_filter_mask(adata, filters).sum())
