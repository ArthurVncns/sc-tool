import io

import numpy as np
import pandas as pd
import anndata as ad
import pytest

from analysis.annotation import run_clustering
from analysis.de import add_obs_column, get_de_results, load_metadata_file, run_de
from analysis.preprocessing import log_transform, normalize_cpm, select_hvg
from analysis.reduction import run_pca, run_umap


def _make_pipeline_adata(n_obs: int = 80, n_vars: int = 200, seed: int = 0) -> ad.AnnData:
    """Return an AnnData that has gone through the full pipeline up to clustering."""
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
    run_clustering(adata, resolution=0.5)
    return adata


# ── run_de ────────────────────────────────────────────────────────────────────

def test_run_de_adds_sc_tool_de_key():
    adata = _make_pipeline_adata()
    clusters = adata.obs["leiden"].unique()
    group = str(clusters[0])
    run_de(adata, groupby="leiden", group=group, reference="rest")
    assert "sc_tool_de" in adata.uns


def test_run_de_stores_correct_metadata():
    adata = _make_pipeline_adata()
    group = str(adata.obs["leiden"].unique()[0])
    run_de(adata, groupby="leiden", group=group, reference="rest", method="wilcoxon")
    de = adata.uns["sc_tool_de"]
    assert de["groupby"] == "leiden"
    assert de["group"] == group
    assert de["reference"] == "rest"
    assert de["method"] == "wilcoxon"


def test_run_de_two_groups():
    adata = _make_pipeline_adata()
    clusters = sorted(adata.obs["leiden"].unique().tolist())
    if len(clusters) < 2:
        pytest.skip("Need at least 2 clusters")
    run_de(adata, groupby="leiden", group=clusters[0], reference=clusters[1])
    assert "sc_tool_de" in adata.uns


# ── get_de_results ────────────────────────────────────────────────────────────

def test_get_de_results_returns_dataframe():
    adata = _make_pipeline_adata()
    group = str(adata.obs["leiden"].unique()[0])
    run_de(adata, groupby="leiden", group=group)
    df = get_de_results(adata)
    assert isinstance(df, pd.DataFrame)


def test_get_de_results_has_expected_columns():
    adata = _make_pipeline_adata()
    group = str(adata.obs["leiden"].unique()[0])
    run_de(adata, groupby="leiden", group=group, n_genes=50)
    df = get_de_results(adata)
    for col in ("names", "scores", "logfoldchanges", "pvals", "pvals_adj"):
        assert col in df.columns, f"Missing column: {col}"


def test_get_de_results_length():
    adata = _make_pipeline_adata()
    group = str(adata.obs["leiden"].unique()[0])
    run_de(adata, groupby="leiden", group=group, n_genes=30)
    df = get_de_results(adata)
    assert len(df) == 30


# ── add_obs_column ────────────────────────────────────────────────────────────

def test_add_obs_column_adds_to_obs():
    adata = _make_pipeline_adata()
    meta = pd.DataFrame({
        "barcode": adata.obs_names[:40].tolist(),
        "condition": ["healthy"] * 20 + ["disease"] * 20,
    })
    add_obs_column(adata, meta, id_column="barcode", label_column="condition", obs_key="cond")
    assert "cond" in adata.obs.columns


def test_add_obs_column_aligns_on_obs_names():
    adata = _make_pipeline_adata()
    # Only annotate first 20 cells; the rest get NaN
    meta = pd.DataFrame({
        "barcode": adata.obs_names[:20].tolist(),
        "group": ["A"] * 20,
    })
    add_obs_column(adata, meta, "barcode", "group", "grp")
    assert adata.obs["grp"].notna().sum() == 20
    assert adata.obs["grp"].isna().sum() == adata.n_obs - 20


# ── load_metadata_file ────────────────────────────────────────────────────────

def test_load_metadata_file_csv():
    csv_content = b"barcode,condition\nCell1,healthy\nCell2,disease\n"
    df = load_metadata_file(io.BytesIO(csv_content))
    assert df.shape == (2, 2)
    assert "condition" in df.columns


def test_load_metadata_file_tsv():
    tsv_content = b"barcode\tcondition\nCell1\thealthy\nCell2\tdisease\n"
    df = load_metadata_file(io.BytesIO(tsv_content))
    assert df.shape == (2, 2)
    assert "condition" in df.columns
