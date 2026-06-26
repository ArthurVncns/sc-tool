import numpy as np
import pandas as pd
import anndata as ad
import pytest

from analysis.preprocessing import (
    log_transform,
    normalize_cpm,
    normalize_pearson,
    reset_preprocessing,
    select_hvg,
    select_hvg_pearson,
)


def _make_adata(n_obs: int = 20, n_vars: int = 50, seed: int = 0) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    X = rng.poisson(10, size=(n_obs, n_vars)).astype(float)
    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"Cell{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(n_vars)]),
    )


# ── Workflow 1: normalize_cpm ─────────────────────────────────────────────────

def test_normalize_cpm_saves_counts_layer():
    adata = _make_adata()
    normalize_cpm(adata)
    assert "counts" in adata.layers


def test_normalize_cpm_counts_layer_is_independent_copy():
    adata = _make_adata()
    original_sum = float(adata.X[0].sum())
    normalize_cpm(adata, target_sum=1_000.0)
    assert abs(float(adata.layers["counts"][0].sum()) - original_sum) < 1e-6


def test_normalize_cpm_scales_to_target_sum():
    adata = _make_adata()
    normalize_cpm(adata, target_sum=1_000.0)
    row_sums = np.asarray(adata.X).sum(axis=1)
    assert np.allclose(row_sums, 1_000.0, rtol=1e-5)


def test_normalize_cpm_sets_workflow_flag():
    adata = _make_adata()
    normalize_cpm(adata)
    assert adata.uns.get("sc_tool_workflow") == "standard"


# ── Workflow 1: log_transform ─────────────────────────────────────────────────

def test_log_transform_sets_log1p_flag():
    adata = _make_adata()
    normalize_cpm(adata)
    log_transform(adata)
    assert "log1p" in adata.uns


def test_log_transform_values_are_log1p():
    adata = _make_adata()
    normalize_cpm(adata)
    original = adata.X.copy()
    log_transform(adata)
    assert np.allclose(adata.X, np.log1p(original), rtol=1e-5)


def test_log_transform_sets_raw():
    adata = _make_adata()
    normalize_cpm(adata)
    assert adata.raw is None
    log_transform(adata)
    assert adata.raw is not None


def test_log_transform_after_pearson_raises():
    adata = _make_adata()
    normalize_pearson(adata)
    with pytest.raises(ValueError, match="Pearson Residuals workflow"):
        log_transform(adata)


# ── Workflow 1: select_hvg ────────────────────────────────────────────────────

def test_select_hvg_seurat_adds_column():
    adata = _make_adata()
    normalize_cpm(adata)
    log_transform(adata)
    select_hvg(adata, n_top_genes=10)
    assert "highly_variable" in adata.var.columns
    assert "dispersions_norm" in adata.var.columns


def test_select_hvg_seurat_respects_n_top_genes():
    adata = _make_adata(n_vars=100)
    normalize_cpm(adata)
    log_transform(adata)
    select_hvg(adata, n_top_genes=20, flavor="seurat")
    assert adata.var["highly_variable"].sum() == 20


def test_select_hvg_seurat_v3_adds_variances_norm():
    adata = _make_adata(n_vars=100)
    normalize_cpm(adata)
    log_transform(adata)
    select_hvg(adata, n_top_genes=20, flavor="seurat_v3")
    assert "variances_norm" in adata.var.columns


def test_select_hvg_clamps_to_n_vars():
    adata = _make_adata(n_vars=30)
    normalize_cpm(adata)
    log_transform(adata)
    select_hvg(adata, n_top_genes=9_999)
    assert adata.var["highly_variable"].sum() <= 30


def test_select_hvg_seurat_v3_without_counts_layer_raises():
    adata = _make_adata()
    with pytest.raises(ValueError, match="counts"):
        select_hvg(adata, flavor="seurat_v3")


def test_select_hvg_modifies_in_place():
    """select_hvg should add highly_variable to the same adata object."""
    adata = _make_adata()
    normalize_cpm(adata)
    log_transform(adata)
    adata_ref = adata  # same object, not a copy
    select_hvg(adata, n_top_genes=10)
    assert "highly_variable" in adata_ref.var.columns


# ── Workflow 2: normalize_pearson ─────────────────────────────────────────────

def test_normalize_pearson_saves_counts_layer():
    adata = _make_adata()
    normalize_pearson(adata)
    assert "counts" in adata.layers


def test_normalize_pearson_sets_workflow_flag():
    adata = _make_adata()
    normalize_pearson(adata)
    assert adata.uns.get("sc_tool_workflow") == "pearson"


def test_normalize_pearson_sets_raw():
    adata = _make_adata()
    normalize_pearson(adata)
    assert adata.raw is not None


def test_normalize_pearson_removes_zero_count_genes():
    """Zero-count genes must be filtered to avoid division by zero in residuals."""
    X = np.zeros((20, 20))
    rng = np.random.default_rng(0)
    X[:, :10] = rng.poisson(5, size=(20, 10))  # first 10 genes expressed
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"Cell{i}" for i in range(20)]),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(20)]),
    )
    normalize_pearson(adata)
    assert adata.n_vars == 10
    assert not np.any(np.isnan(adata.X))


# ── Workflow 2: select_hvg_pearson ────────────────────────────────────────────

def test_select_hvg_pearson_adds_column():
    adata = _make_adata(n_vars=100)
    normalize_pearson(adata)
    select_hvg_pearson(adata, n_top_genes=20)
    assert "highly_variable" in adata.var.columns
    assert "residual_variances" in adata.var.columns


def test_select_hvg_pearson_respects_n_top_genes():
    adata = _make_adata(n_vars=100)
    normalize_pearson(adata)
    select_hvg_pearson(adata, n_top_genes=20)
    assert adata.var["highly_variable"].sum() == 20


def test_select_hvg_pearson_without_normalization_raises():
    adata = _make_adata()
    with pytest.raises(ValueError, match="normalize_pearson"):
        select_hvg_pearson(adata)


def test_select_hvg_pearson_after_standard_workflow_raises():
    adata = _make_adata()
    normalize_cpm(adata)
    with pytest.raises(ValueError, match="normalize_pearson"):
        select_hvg_pearson(adata)


# ── reset_preprocessing ───────────────────────────────────────────────────────

def test_reset_restores_raw_counts():
    adata = _make_adata()
    original_X = adata.X.copy()
    normalize_cpm(adata, target_sum=1_000.0)
    reset_preprocessing(adata)
    assert np.allclose(adata.X, original_X)


def test_reset_clears_workflow_flag():
    adata = _make_adata()
    normalize_cpm(adata)
    reset_preprocessing(adata)
    assert "sc_tool_workflow" not in adata.uns


def test_reset_clears_log1p():
    adata = _make_adata()
    normalize_cpm(adata)
    log_transform(adata)
    reset_preprocessing(adata)
    assert "log1p" not in adata.uns


def test_reset_clears_hvg_columns():
    adata = _make_adata(n_vars=100)
    normalize_cpm(adata)
    log_transform(adata)
    select_hvg(adata, n_top_genes=20)
    assert "highly_variable" in adata.var.columns
    reset_preprocessing(adata)
    assert "highly_variable" not in adata.var.columns
    assert "dispersions_norm" not in adata.var.columns


def test_reset_clears_raw():
    adata = _make_adata()
    normalize_cpm(adata)
    log_transform(adata)
    assert adata.raw is not None
    reset_preprocessing(adata)
    assert adata.raw is None


def test_reset_allows_workflow_switch():
    """After reset, a different workflow can be started."""
    adata = _make_adata(n_vars=100)
    normalize_cpm(adata)
    reset_preprocessing(adata)
    # Should now be able to start the Pearson workflow
    normalize_pearson(adata)
    assert adata.uns.get("sc_tool_workflow") == "pearson"


def test_reset_without_counts_layer_raises():
    adata = _make_adata()
    with pytest.raises(ValueError, match="adata.layers"):
        reset_preprocessing(adata)
