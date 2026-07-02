import numpy as np
import pandas as pd
import anndata as ad
import pytest

from analysis.preprocessing import log_transform, normalize_cpm, select_hvg
from analysis.reduction import run_harmony, run_pca, run_umap


def _make_preprocessed_adata(n_obs: int = 80, n_vars: int = 200, seed: int = 0) -> ad.AnnData:
    """Return a fully preprocessed AnnData ready for PCA."""
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
    return adata


# --- PCA ---

def test_run_pca_adds_obsm_key():
    adata = _make_preprocessed_adata()
    assert "X_pca" not in adata.obsm
    run_pca(adata, n_comps=10)
    assert "X_pca" in adata.obsm


def test_run_pca_shape_matches_n_comps():
    adata = _make_preprocessed_adata()
    run_pca(adata, n_comps=10)
    assert adata.obsm["X_pca"].shape == (adata.n_obs, 10)


def test_run_pca_variance_ratio_in_uns():
    adata = _make_preprocessed_adata()
    run_pca(adata, n_comps=10)
    assert "pca" in adata.uns
    assert "variance_ratio" in adata.uns["pca"]
    assert len(adata.uns["pca"]["variance_ratio"]) == 10


def test_run_pca_clamps_n_comps():
    """n_comps should be silently clamped, not raise an error."""
    adata = _make_preprocessed_adata(n_obs=30, n_vars=50)
    select_hvg(adata, n_top_genes=20)  # only 20 HVGs
    run_pca(adata, n_comps=9_999)  # absurdly large
    assert adata.obsm["X_pca"].shape[1] <= min(30, 20) - 1


def test_run_pca_invalidates_umap():
    """Re-running PCA must clear any stale UMAP embedding."""
    adata = _make_preprocessed_adata()
    run_pca(adata, n_comps=10)
    run_umap(adata, n_pcs=10, n_neighbors=10)
    assert "X_umap" in adata.obsm
    # Now re-run PCA — UMAP should be gone
    run_pca(adata, n_comps=8)
    assert "X_umap" not in adata.obsm
    assert "neighbors" not in adata.uns


# --- UMAP ---

def test_run_umap_adds_obsm_key():
    adata = _make_preprocessed_adata()
    run_pca(adata, n_comps=10)
    run_umap(adata, n_pcs=10, n_neighbors=10)
    assert "X_umap" in adata.obsm


def test_run_umap_embedding_is_2d():
    adata = _make_preprocessed_adata()
    run_pca(adata, n_comps=10)
    run_umap(adata, n_pcs=10, n_neighbors=10)
    assert adata.obsm["X_umap"].shape == (adata.n_obs, 2)


def test_run_umap_clamps_n_pcs():
    """n_pcs larger than available PCA components should not raise an error."""
    adata = _make_preprocessed_adata()
    run_pca(adata, n_comps=10)
    run_umap(adata, n_pcs=9_999, n_neighbors=10)  # will be clamped to 10
    assert adata.obsm["X_umap"].shape[1] == 2


def test_run_umap_is_reproducible():
    """Same random_state must produce identical embeddings."""
    adata = _make_preprocessed_adata()
    run_pca(adata, n_comps=10)

    run_umap(adata, n_pcs=10, n_neighbors=10, random_state=42)
    coords_a = adata.obsm["X_umap"].copy()

    run_umap(adata, n_pcs=10, n_neighbors=10, random_state=42)
    coords_b = adata.obsm["X_umap"].copy()

    assert np.allclose(coords_a, coords_b)


# --- Harmony ---

def _make_adata_with_batch(n_obs: int = 80, n_vars: int = 200, seed: int = 0) -> ad.AnnData:
    """Return a preprocessed AnnData with a 'batch' column in obs."""
    rng = np.random.default_rng(seed)
    X = rng.poisson(5, size=(n_obs, n_vars)).astype(float)
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(
            {"batch": ["A"] * (n_obs // 2) + ["B"] * (n_obs // 2)},
            index=[f"Cell{i}" for i in range(n_obs)],
        ),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(n_vars)]),
    )
    normalize_cpm(adata)
    log_transform(adata)
    select_hvg(adata, n_top_genes=50)
    run_pca(adata, n_comps=20)
    return adata


def test_run_harmony_adds_corrected_embedding():
    adata = _make_adata_with_batch()
    run_harmony(adata, batch_key="batch")
    assert "X_pca_harmony" in adata.obsm


def test_run_harmony_shape_matches_pca():
    adata = _make_adata_with_batch()
    run_harmony(adata, batch_key="batch")
    assert adata.obsm["X_pca_harmony"].shape == adata.obsm["X_pca"].shape


def test_run_harmony_stores_batch_key():
    adata = _make_adata_with_batch()
    run_harmony(adata, batch_key="batch")
    assert adata.uns.get("harmony_batch_key") == "batch"


def test_run_harmony_invalidates_umap():
    adata = _make_adata_with_batch()
    run_pca(adata, n_comps=10)
    run_umap(adata, n_pcs=10, n_neighbors=10)
    assert "X_umap" in adata.obsm
    run_harmony(adata, batch_key="batch")
    assert "X_umap" not in adata.obsm


def test_run_pca_invalidates_harmony():
    adata = _make_adata_with_batch()
    run_harmony(adata, batch_key="batch")
    assert "X_pca_harmony" in adata.obsm
    run_pca(adata, n_comps=10)  # re-run PCA
    assert "X_pca_harmony" not in adata.obsm


def test_run_umap_uses_harmony_when_available():
    adata = _make_adata_with_batch()
    run_harmony(adata, batch_key="batch")
    run_umap(adata, n_neighbors=10)
    assert "X_umap" in adata.obsm
    assert adata.obsm["X_umap"].shape == (adata.n_obs, 2)
