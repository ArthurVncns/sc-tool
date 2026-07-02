import anndata as ad
import scanpy as sc


def run_scrublet(
    adata: ad.AnnData,
    expected_doublet_rate: float = 0.05,
    random_state: int = 0,
) -> ad.AnnData:
    """Detect potential doublets using Scanpy's built-in Scrublet implementation.

    Simulates artificial doublets by combining pairs of cells, then assigns
    each real cell a score based on its similarity to the simulated doublets.

    Uses adata.layers['counts'] (raw integer counts) when available — the
    doublet detection must run on unnormalized data. Falls back to adata.X.

    Adds to adata.obs:
        - 'doublet_score':     continuous score [0, 1]; higher = more likely doublet.
        - 'predicted_doublet': boolean flag for cells exceeding the auto-threshold.

    Args:
        adata: AnnData. Raw counts should be in adata.layers['counts'] or adata.X.
        expected_doublet_rate: Fraction of doublets expected in the dataset.
            Typical 10x values: ~0.04 for 4k cells, ~0.08 for 8k cells.
        random_state: Seed for reproducibility.

    Returns:
        adata modified in place.
    """
    # Run on raw counts — doublet simulation requires unnormalized data.
    if "counts" in adata.layers:
        from anndata import AnnData
        adata_raw = AnnData(X=adata.layers["counts"], obs=adata.obs, var=adata.var)
    else:
        adata_raw = adata

    sc.pp.scrublet(
        adata_raw,
        expected_doublet_rate=expected_doublet_rate,
        random_state=random_state,
        verbose=False,
    )

    # Copy results back to the main adata
    adata.obs["doublet_score"] = adata_raw.obs["doublet_score"].values
    adata.obs["predicted_doublet"] = adata_raw.obs["predicted_doublet"].values

    # Persist the count so it remains readable after doublets are filtered out.
    adata.uns["sc_tool_n_doublets"] = int(adata.obs["predicted_doublet"].sum())

    return adata
