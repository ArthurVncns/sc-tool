import anndata as ad
import pandas as pd
import scanpy as sc

# Curated list of CellTypist models useful for common tissue types.
# Models are auto-downloaded on first use (~100 MB each, cached locally).
CELLTYPIST_MODELS: dict[str, str] = {
    "Immune — broad (pan-tissue)": "Immune_All_Low.pkl",
    "Immune — detailed (pan-tissue)": "Immune_All_High.pkl",
    "Human Lung Atlas": "Human_Lung_Atlas.pkl",
    "Developing Human Brain": "Developing_Human_Brain.pkl",
    "Mouse Gut (10x)": "Mouse_Gut_10x.pkl",
    "Pan-fetal Human": "Pan_Fetal_Human.pkl",
}

# A ready-to-use PBMC marker dictionary for the marker gene scoring method.
PBMC_MARKERS: dict[str, list[str]] = {
    "T cells": ["CD3D", "CD3E", "CD3G", "TRAC"],
    "CD4+ T cells": ["CD4", "IL7R", "CCR7", "S100A4"],
    "CD8+ T cells": ["CD8A", "CD8B", "GZMK", "GZMB"],
    "B cells": ["CD19", "MS4A1", "CD79A", "CD79B"],
    "NK cells": ["GNLY", "NKG7", "KLRD1", "NCAM1"],
    "Monocytes": ["CD14", "LYZ", "CST3", "FCGR3A"],
    "Dendritic cells": ["FCER1A", "CLEC10A", "HLA-DQA1"],
    "Platelets": ["PPBP", "PF4", "GP1BB"],
}


def run_clustering(adata: ad.AnnData, resolution: float = 0.5) -> ad.AnnData:
    """Cluster cells using the Leiden algorithm on the precomputed neighborhood graph.

    Stores cluster labels as strings in adata.obs['leiden'].
    Invalidates any previously computed marker genes and cell type annotations,
    since those were derived from the old cluster assignments.

    Args:
        adata: AnnData with 'neighbors' in adata.uns (built by sc.pp.neighbors).
        resolution: Controls the coarseness of the clustering.
            Higher values produce more, smaller clusters.

    Returns:
        adata modified in place.
    """
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2)

    # Invalidate results that depended on the previous cluster assignments.
    adata.uns.pop("rank_genes_groups", None)
    if "cell_type" in adata.obs.columns:
        del adata.obs["cell_type"]

    return adata


def find_marker_genes(
    adata: ad.AnnData,
    groupby: str = "leiden",
    n_genes: int = 25,
) -> ad.AnnData:
    """Find marker genes for each cluster using the Wilcoxon rank-sum test.

    Each cluster is compared against all other cells. Results are stored in
    adata.uns['rank_genes_groups'] and can be retrieved with get_top_markers_df().

    Args:
        adata: AnnData with the groupby column in adata.obs.
        groupby: obs column to group cells by. Default 'leiden'.
        n_genes: Number of top marker genes to compute per group.

    Returns:
        adata modified in place.
    """
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method="wilcoxon",
        n_genes=n_genes,
    )
    return adata


def get_top_markers_df(adata: ad.AnnData, n_genes: int = 5) -> pd.DataFrame:
    """Return a DataFrame of the top marker genes per cluster.

    Columns are cluster names; rows are ranked marker genes (Rank 1 = best).

    Args:
        adata: AnnData with 'rank_genes_groups' in adata.uns.
        n_genes: Number of top markers to show per cluster.

    Returns:
        DataFrame with shape (n_genes, n_clusters).
    """
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names  # cluster names as a tuple
    return pd.DataFrame(
        {group: result["names"][group][:n_genes] for group in groups},
        index=[f"Rank {i + 1}" for i in range(n_genes)],
    )


def apply_annotations(
    adata: ad.AnnData,
    annotation_map: dict[str, str],
) -> ad.AnnData:
    """Assign cell type labels to clusters.

    Adds adata.obs['cell_type'] by mapping each cluster label to a cell type name.

    Args:
        adata: AnnData with 'leiden' in adata.obs.
        annotation_map: Mapping from cluster label to cell type name.
            Example: {"0": "T cells", "1": "B cells", "2": "NK cells"}.

    Returns:
        adata modified in place.
    """
    adata.obs["cell_type"] = (
        adata.obs["leiden"].map(annotation_map).astype("category")
    )
    return adata


def run_celltypist(
    adata: ad.AnnData,
    model_name: str = "Immune_All_Low.pkl",
    majority_voting: bool = True,
) -> ad.AnnData:
    """Annotate cell types using a CellTypist pre-trained model.

    The model is downloaded automatically on first use and cached locally.
    Adds 'celltypist_cell_type' to adata.obs.

    CellTypist expects log-normalized data (Standard workflow: CPM + log1p).
    Results may be less accurate on Pearson residuals.

    Args:
        adata: AnnData with log-normalized expression in adata.X.
        model_name: Filename of the CellTypist model (see CELLTYPIST_MODELS).
        majority_voting: If True, refines predictions using cluster-level
            majority voting. Requires 'leiden' in adata.obs.

    Returns:
        adata modified in place.
    """
    import celltypist
    from celltypist import models

    model = models.Model.load(model=model_name)

    over_clustering = (
        "leiden"
        if majority_voting and "leiden" in adata.obs.columns
        else None
    )

    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=majority_voting,
        over_clustering=over_clustering,
    )

    label_col = "majority_voting" if majority_voting else "predicted_labels"
    adata.obs["celltypist_cell_type"] = (
        predictions.predicted_labels[label_col].values
    )
    adata.obs["celltypist_cell_type"] = adata.obs["celltypist_cell_type"].astype(
        "category"
    )
    return adata


def score_marker_genes(
    adata: ad.AnnData,
    marker_dict: dict[str, list[str]],
) -> ad.AnnData:
    """Assign cell types by scoring each cell against custom marker gene lists.

    For each cell type in marker_dict, computes a score using sc.tl.score_genes
    (average expression of the gene set, corrected for background).
    Each cell is assigned to the cell type with the highest score.

    Adds per-cell-type score columns and 'marker_score_cell_type' to adata.obs.

    Args:
        adata: AnnData with log-normalized expression (Standard workflow recommended).
        marker_dict: Cell type name → list of marker genes.
            Genes absent from the dataset are silently skipped.

    Returns:
        adata modified in place.

    Raises:
        ValueError: If none of the provided genes are found in the dataset.
    """
    score_map: dict[str, str] = {}  # cell_type → score column name

    for cell_type, genes in marker_dict.items():
        valid = [g for g in genes if g in adata.var_names]
        if not valid:
            continue
        col = f"score_{cell_type.replace(' ', '_').replace('+', 'pos').replace('-', 'neg')}"
        sc.tl.score_genes(adata, gene_list=valid, score_name=col)
        score_map[cell_type] = col

    if not score_map:
        raise ValueError(
            "None of the provided marker genes were found in the dataset. "
            "Check gene names — they must match the var_names of your AnnData exactly."
        )

    scores = adata.obs[[col for col in score_map.values()]]
    best_col = scores.idxmax(axis=1)
    col_to_type = {col: ct for ct, col in score_map.items()}
    adata.obs["marker_score_cell_type"] = best_col.map(col_to_type).astype("category")

    return adata
