import pandas as pd
import numpy as np
import anndata as ad

from analysis.io import load_h5ad


def test_load_h5ad_roundtrip(tmp_path):
    """load_h5ad should return an AnnData with the same shape as the original."""
    original = ad.AnnData(np.ones((5, 10)))
    h5ad_path = tmp_path / "test.h5ad"
    original.write_h5ad(h5ad_path)

    with open(h5ad_path, "rb") as f:
        loaded = load_h5ad(f)

    assert loaded.shape == (5, 10)


def test_load_h5ad_preserves_obs_columns(tmp_path):
    """load_h5ad should preserve cell metadata (obs) columns."""
    original = ad.AnnData(
        np.ones((3, 4)),
        obs=pd.DataFrame({"cell_type": ["T cell", "B cell", "NK cell"]}),
    )
    h5ad_path = tmp_path / "test_obs.h5ad"
    original.write_h5ad(h5ad_path)

    with open(h5ad_path, "rb") as f:
        loaded = load_h5ad(f)

    assert "cell_type" in loaded.obs.columns
    assert list(loaded.obs["cell_type"]) == ["T cell", "B cell", "NK cell"]
