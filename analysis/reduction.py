import anndata as ad
import scanpy as sc


def run_pca(adata: ad.AnnData, n_comps: int = 50) -> ad.AnnData:
    """Compute PCA on highly variable genes.

    Reduces the HVG expression matrix to n_comps principal components.
    Results are stored in adata.obsm['X_pca'] and adata.uns['pca'].

    Invalidates any previously computed Harmony correction, neighborhood
    graph, and UMAP, because those were derived from the old PCA coordinates.

    Args:
        adata: AnnData with highly variable genes marked in adata.var.
        n_comps: Number of principal components to compute.
            Automatically clamped so it never exceeds min(n_cells, n_HVGs) - 1.

    Returns:
        adata modified in place.
    """
    n_hvg = (
        int(adata.var["highly_variable"].sum())
        if "highly_variable" in adata.var.columns
        else adata.n_vars
    )
    n_comps = min(n_comps, adata.n_obs - 1, n_hvg - 1)

    sc.pp.pca(adata, n_comps=n_comps, mask_var="highly_variable")

    # Invalidate all downstream results that depend on the PCA embedding.
    adata.obsm.pop("X_pca_harmony", None)
    adata.uns.pop("harmony_batch_key", None)
    adata.uns.pop("neighbors", None)
    adata.obsm.pop("X_umap", None)

    return adata


def run_harmony(adata: ad.AnnData, batch_key: str) -> ad.AnnData:
    """Correct the PCA embedding for batch effects using Harmony.

    Creates adata.obsm['X_pca_harmony'] — a batch-corrected version of
    adata.obsm['X_pca']. After running Harmony, run_umap() will
    automatically use the corrected embedding.

    Invalidates any previously computed neighborhood graph and UMAP.

    Args:
        adata: AnnData with PCA computed in adata.obsm['X_pca'].
        batch_key: obs column identifying batch membership (e.g. 'sample', 'donor').

    Returns:
        adata modified in place.

    Raises:
        ImportError: If harmonypy is not installed.
        KeyError: If batch_key is not in adata.obs.
    """
    try:
        import harmonypy as hm
    except ImportError as exc:
        raise ImportError(
            "Harmony requires harmonypy. Install with: pip install harmonypy"
        ) from exc

    if batch_key not in adata.obs.columns:
        raise KeyError(f"Batch key '{batch_key}' not found in adata.obs.")

    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, batch_key)

    # harmonypy 2.x changed Z_corr orientation relative to 0.x;
    # ensure (n_cells, n_pcs) regardless of version.
    import numpy as np
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != adata.n_obs:
        Z = Z.T
    adata.obsm["X_pca_harmony"] = Z
    adata.uns["harmony_batch_key"] = batch_key

    # Invalidate downstream results that were built on uncorrected coordinates.
    adata.uns.pop("neighbors", None)
    adata.obsm.pop("X_umap", None)

    return adata


def run_umap(
    adata: ad.AnnData,
    n_neighbors: int = 15,
    n_pcs: int = 40,
    min_dist: float = 0.1,
    random_state: int = 0,
) -> ad.AnnData:
    """Build a neighborhood graph, then compute a UMAP embedding.

    Automatically uses the Harmony-corrected embedding (X_pca_harmony) if
    Harmony has been run, otherwise falls back to X_pca with n_pcs components.

    Args:
        adata: AnnData with PCA (and optionally Harmony) already computed.
        n_neighbors: Number of neighbors per cell in the KNN graph.
        n_pcs: PCA components to use. Ignored when Harmony is active.
        min_dist: Minimum distance between points in 2D UMAP space.
        random_state: Seed for reproducibility.

    Returns:
        adata modified in place.
    """
    if "X_pca_harmony" in adata.obsm:
        # Harmony embedding already encodes all components — n_pcs is not applicable.
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca_harmony")
    else:
        n_pcs = min(n_pcs, adata.obsm["X_pca"].shape[1])
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    sc.tl.umap(adata, min_dist=min_dist, random_state=random_state)

    return adata
