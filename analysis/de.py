"""Differential expression analysis.

Typical workflow:
    1. add_obs_column()   — optional, only when loading custom metadata
    2. run_de()           — runs Wilcoxon / t-test between two groups
    3. get_de_results()   — retrieves the ranked gene table as a DataFrame
"""

import io as _io

import anndata as ad
import pandas as pd
import scanpy as sc


def load_metadata_file(source) -> pd.DataFrame:
    """Parse a CSV or TSV metadata file, auto-detecting the separator.

    Args:
        source: Binary-readable object (e.g. Streamlit UploadedFile).

    Returns:
        DataFrame with one row per cell.
    """
    content = source.read()
    return pd.read_csv(_io.BytesIO(content), sep=None, engine="python")


def add_obs_column(
    adata: ad.AnnData,
    metadata: pd.DataFrame,
    id_column: str,
    label_column: str,
    obs_key: str,
) -> ad.AnnData:
    """Merge an external group label into adata.obs.

    Aligns on adata.obs_names using id_column as the index.
    Cells not present in the metadata receive NaN.

    Args:
        adata: AnnData to annotate.
        metadata: DataFrame loaded by load_metadata_file().
        id_column: Column in metadata that matches adata.obs_names.
        label_column: Column in metadata containing group labels.
        obs_key: Name of the new column in adata.obs.

    Returns:
        adata modified in place.
    """
    indexed = metadata.set_index(id_column)
    adata.obs[obs_key] = indexed[label_column].reindex(adata.obs_names)
    return adata


def run_de(
    adata: ad.AnnData,
    groupby: str,
    group: str,
    reference: str = "rest",
    method: str = "wilcoxon",
    n_genes: int = 200,
) -> ad.AnnData:
    """Run differential expression analysis between two groups of cells.

    Compares cells labelled `group` against cells labelled `reference`
    (or all other cells when reference='rest').

    Uses adata.raw when available (log-normalized, full-gene matrix set
    by log_transform() or normalize_pearson()). Falls back to adata.X.

    Results are stored in two places:
    - adata.uns['rank_genes_groups'] — Scanpy's standard output key
    - adata.uns['sc_tool_de']       — our own cache (a dict of lists),
      so the results survive if rank_genes_groups is later overwritten
      by the annotation page's marker gene computation.

    Args:
        adata: Preprocessed AnnData.
        groupby: obs column that defines the groups (e.g. 'condition').
        group: The group to test (e.g. 'malignant').
        reference: Reference group (e.g. 'healthy') or 'rest'.
        method: Statistical test — 'wilcoxon' (recommended), 't-test', 'logreg'.
        n_genes: Number of top-ranked genes to store.

    Returns:
        adata modified in place.
    """
    use_raw = adata.raw is not None

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=[group],
        reference=reference,
        method=method,
        n_genes=n_genes,
        use_raw=use_raw,
    )

    results_df = sc.get.rank_genes_groups_df(adata, group=group)

    adata.uns["sc_tool_de"] = {
        "groupby": groupby,
        "group": group,
        "reference": reference,
        "method": method,
        "results": results_df.to_dict(orient="list"),
    }

    return adata


def get_de_results(adata: ad.AnnData) -> pd.DataFrame:
    """Return the cached DE results as a tidy DataFrame.

    Columns: names, scores, logfoldchanges, pvals, pvals_adj.
    """
    de = adata.uns["sc_tool_de"]
    return pd.DataFrame(de["results"])
