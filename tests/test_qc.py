import numpy as np
import pandas as pd
import anndata as ad

from analysis.qc import QCFilters, compute_qc_metrics, count_cells_passing_filters, filter_cells


def _make_adata(
    n_obs: int = 10,
    n_vars: int = 20,
    include_mt: bool = False,
    include_ribo: bool = False,
    include_hb: bool = False,
    seed: int = 42,
) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    X = rng.poisson(5, size=(n_obs, n_vars)).astype(float)
    var_names = [f"Gene{i}" for i in range(n_vars)]
    if include_mt:
        var_names[-3:] = ["MT-ND1", "MT-ND2", "MT-CO1"]
    if include_ribo:
        var_names[-6:-3] = ["RPS1", "RPL2", "RPS3"]
    if include_hb:
        var_names[-9:-6] = ["HBA1", "HBB", "HBD"]
    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"Cell{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=var_names),
    )


def test_compute_qc_metrics_adds_expected_columns():
    adata = compute_qc_metrics(_make_adata())
    assert "n_genes_by_counts" in adata.obs.columns
    assert "total_counts" in adata.obs.columns


def test_compute_qc_metrics_detects_mt_genes():
    adata = compute_qc_metrics(_make_adata(include_mt=True))
    assert "pct_counts_mt" in adata.obs.columns


def test_compute_qc_metrics_no_mt_genes():
    adata = compute_qc_metrics(_make_adata(include_mt=False))
    assert "pct_counts_mt" not in adata.obs.columns


def test_compute_qc_metrics_detects_ribo_genes():
    adata = compute_qc_metrics(_make_adata(include_ribo=True))
    assert "pct_counts_ribo" in adata.obs.columns


def test_compute_qc_metrics_detects_hb_genes():
    adata = compute_qc_metrics(_make_adata(include_hb=True))
    assert "pct_counts_hb" in adata.obs.columns


def test_compute_qc_metrics_no_ribo_or_hb():
    adata = compute_qc_metrics(_make_adata())
    assert "pct_counts_ribo" not in adata.obs.columns
    assert "pct_counts_hb" not in adata.obs.columns


def test_filter_cells_removes_cells():
    """filter_cells should keep only cells within the specified thresholds."""
    # Build an adata with known QC values directly so the test is explicit.
    # Cells 0-2 have low gene counts; cells 3-5 have high gene counts.
    adata = ad.AnnData(
        np.ones((6, 5)),
        obs=pd.DataFrame(
            {
                "n_genes_by_counts": [100, 200, 300, 5_000, 6_000, 7_000],
                "total_counts":      [500, 1_000, 1_500, 50_000, 60_000, 70_000],
            },
            index=[f"Cell{i}" for i in range(6)],
        ),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(5)]),
    )
    filters = QCFilters(min_genes=200, max_genes=4_000, min_counts=0, max_counts=99_999)
    filtered = filter_cells(adata, filters)

    assert filtered.n_obs == 2
    assert list(filtered.obs_names) == ["Cell1", "Cell2"]


def test_filter_cells_applies_pct_mt_threshold():
    adata = compute_qc_metrics(_make_adata(include_mt=True, n_obs=10))
    # max_pct_mt=0 should filter out any cell with any mitochondrial expression
    filters = QCFilters(min_genes=0, max_genes=9_999, min_counts=0, max_counts=9_999_999, max_pct_mt=0.0)
    filtered = filter_cells(adata, filters)
    assert all(filtered.obs["pct_counts_mt"] == 0.0)


def test_filter_cells_returns_copy():
    adata = compute_qc_metrics(_make_adata())
    filtered = filter_cells(adata, QCFilters())
    filtered.obs["test_col"] = 1
    assert "test_col" not in adata.obs.columns


def test_count_cells_passing_filters_matches_filter_cells():
    adata = compute_qc_metrics(_make_adata(n_obs=30, seed=1))
    filters = QCFilters(min_genes=5, max_genes=20, min_counts=10, max_counts=500)
    assert count_cells_passing_filters(adata, filters) == filter_cells(adata, filters).n_obs
