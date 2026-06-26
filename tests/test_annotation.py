import numpy as np
import pandas as pd
import anndata as ad
import pytest

from analysis.annotation import (
    apply_annotations,
    find_marker_genes,
    get_top_markers_df,
    run_clustering,
)
from analysis.preprocessing import log_transform, normalize_cpm, select_hvg
from analysis.reduction import run_pca, run_umap


def _make_analysis_ready_adata(n_obs: int = 100, n_vars: int = 200, seed: int = 0) -> ad.AnnData:
    """Return an AnnData that has gone through the full pipeline up to UMAP."""
    rng = np.random.default_rng(seed)
    X = rng.poisson(5, size=(n_obs, n_vars)).astype(float)
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"Cell{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(n_vars)]),
    )
    normalize_cpm(adata)
    log_transform(adata)
    select_hvg(adata, n_top_genes=50)
    run_pca(adata, n_comps=20)
    run_umap(adata, n_pcs=20, n_neighbors=10)
    return adata


# --- run_clustering ---

def test_clustering_adds_leiden_column():
    adata = _make_analysis_ready_adata()
    run_clustering(adata, resolution=0.3)
    assert "leiden" in adata.obs.columns


def test_clustering_produces_string_labels():
    adata = _make_analysis_ready_adata()
    run_clustering(adata, resolution=0.3)
    # Leiden labels are strings: "0", "1", "2", ...
    assert all(isinstance(label, str) for label in adata.obs["leiden"].unique())


def test_reclustering_invalidates_markers():
    adata = _make_analysis_ready_adata()
    run_clustering(adata, resolution=0.3)
    find_marker_genes(adata)
    assert "rank_genes_groups" in adata.uns
    run_clustering(adata, resolution=0.5)  # re-cluster
    assert "rank_genes_groups" not in adata.uns


def test_reclustering_invalidates_cell_type():
    adata = _make_analysis_ready_adata()
    run_clustering(adata, resolution=0.3)
    adata.obs["cell_type"] = "T cells"  # mock annotation
    run_clustering(adata, resolution=0.5)
    assert "cell_type" not in adata.obs.columns


# --- find_marker_genes ---

def test_find_marker_genes_adds_uns_key():
    adata = _make_analysis_ready_adata()
    run_clustering(adata, resolution=0.5)
    find_marker_genes(adata)
    assert "rank_genes_groups" in adata.uns


def test_find_marker_genes_has_names_key():
    adata = _make_analysis_ready_adata()
    run_clustering(adata, resolution=0.5)
    find_marker_genes(adata, n_genes=10)
    assert "names" in adata.uns["rank_genes_groups"]


# --- get_top_markers_df ---

def test_get_top_markers_df_shape():
    adata = _make_analysis_ready_adata()
    run_clustering(adata, resolution=0.5)
    find_marker_genes(adata, n_genes=20)
    n_clusters = adata.obs["leiden"].nunique()
    df = get_top_markers_df(adata, n_genes=5)
    assert df.shape == (5, n_clusters)


def test_get_top_markers_df_columns_are_cluster_names():
    adata = _make_analysis_ready_adata()
    run_clustering(adata, resolution=0.5)
    find_marker_genes(adata)
    df = get_top_markers_df(adata, n_genes=3)
    clusters = set(adata.obs["leiden"].unique())
    assert set(df.columns) == clusters


# --- apply_annotations ---

def test_apply_annotations_adds_cell_type_column():
    adata = _make_analysis_ready_adata()
    run_clustering(adata, resolution=0.3)
    clusters = adata.obs["leiden"].unique()
    annotation_map = {c: f"CellType_{c}" for c in clusters}
    apply_annotations(adata, annotation_map)
    assert "cell_type" in adata.obs.columns


def test_apply_annotations_maps_correctly():
    adata = _make_analysis_ready_adata()
    run_clustering(adata, resolution=0.3)
    clusters = adata.obs["leiden"].unique()
    annotation_map = {c: f"Type_{c}" for c in clusters}
    apply_annotations(adata, annotation_map)
    for cluster, cell_type in annotation_map.items():
        mask = adata.obs["leiden"] == cluster
        assert (adata.obs.loc[mask, "cell_type"] == cell_type).all()


def test_apply_annotations_result_is_categorical():
    adata = _make_analysis_ready_adata()
    run_clustering(adata, resolution=0.3)
    clusters = adata.obs["leiden"].unique()
    apply_annotations(adata, {c: "T cells" for c in clusters})
    assert str(adata.obs["cell_type"].dtype) == "category"
