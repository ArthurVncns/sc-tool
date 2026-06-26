import anndata as ad
import scanpy as sc


def run_pca(adata: ad.AnnData, n_comps: int = 50) -> ad.AnnData:
    """Compute PCA on highly variable genes.

    Reduces the HVG expression matrix to n_comps principal components.
    Results are stored in adata.obsm['X_pca'] and adata.uns['pca'].

    Invalidates any previously computed neighborhood graph and UMAP, because
    those were derived from the old PCA coordinates.

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

    # Invalidate downstream results that depend on these PCA coordinates.
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
    """Build a neighborhood graph on PCA coordinates, then compute a UMAP embedding.

    Two steps happen internally:
    1. sc.pp.neighbors: builds a K-nearest-neighbor graph using the first
       n_pcs principal components. The graph structure drives the UMAP layout.
    2. sc.tl.umap: projects cells to 2D. Results stored in adata.obsm['X_umap'].

    Args:
        adata: AnnData with PCA already computed in adata.obsm['X_pca'].
        n_neighbors: Number of neighbors per cell in the KNN graph.
            Low values → fine local structure. High values → broader global structure.
        n_pcs: How many PCA components to use for building the graph.
            Informed by the PCA elbow plot — choose where the curve flattens.
        min_dist: Minimum distance between points in 2D UMAP space.
            Low values → tight clusters. High values → even spread.
        random_state: Seed for reproducibility.

    Returns:
        adata modified in place.
    """
    n_pcs_available = adata.obsm["X_pca"].shape[1]
    n_pcs = min(n_pcs, n_pcs_available)

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata, min_dist=min_dist, random_state=random_state)

    return adata
