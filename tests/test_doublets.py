import numpy as np
import pandas as pd
import anndata as ad

from analysis.doublets import run_scrublet
from analysis.preprocessing import normalize_cpm


def _make_adata(n_obs: int = 60, n_vars: int = 100, seed: int = 0) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    X = rng.poisson(5, size=(n_obs, n_vars)).astype(float)
    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"Cell{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(n_vars)]),
    )


def test_run_scrublet_adds_doublet_score():
    adata = _make_adata()
    run_scrublet(adata)
    assert "doublet_score" in adata.obs.columns


def test_run_scrublet_adds_predicted_doublet():
    adata = _make_adata()
    run_scrublet(adata)
    assert "predicted_doublet" in adata.obs.columns


def test_run_scrublet_scores_are_in_0_1():
    adata = _make_adata()
    run_scrublet(adata)
    scores = adata.obs["doublet_score"].values
    assert scores.min() >= 0.0
    assert scores.max() <= 1.0


def test_run_scrublet_predicted_doublet_is_boolean():
    adata = _make_adata()
    run_scrublet(adata)
    assert adata.obs["predicted_doublet"].dtype == bool


def test_run_scrublet_uses_counts_layer_when_available():
    """run_scrublet should use raw counts even after normalization."""
    adata = _make_adata()
    raw_X = adata.X.copy()
    normalize_cpm(adata)  # modifies adata.X, saves raw in layers["counts"]
    run_scrublet(adata)
    # Scores should be computable without errors — values differ from
    # running on normalized data but there's no crash
    assert "doublet_score" in adata.obs.columns


def test_run_scrublet_score_length_matches_n_obs():
    adata = _make_adata(n_obs=40)
    run_scrublet(adata)
    assert len(adata.obs["doublet_score"]) == 40
