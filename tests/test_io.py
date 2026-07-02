import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad

from analysis.io import load_h5ad, load_10x_mtx


# ── h5ad ──────────────────────────────────────────────────────────────────────

def test_load_h5ad_roundtrip(tmp_path):
    original = ad.AnnData(np.ones((5, 10)))
    h5ad_path = tmp_path / "test.h5ad"
    original.write_h5ad(h5ad_path)

    with open(h5ad_path, "rb") as f:
        loaded = load_h5ad(f)

    assert loaded.shape == (5, 10)


def test_load_h5ad_preserves_obs_columns(tmp_path):
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


# ── 10x MTX ───────────────────────────────────────────────────────────────────

def _write_10x_mtx(folder: Path, n_cells: int = 4, n_genes: int = 6) -> None:
    """Write a minimal 10x MTX folder (gzip-compressed) for testing."""
    import gzip
    import io as _io

    rng = np.random.default_rng(0)
    X = sp.random(n_genes, n_cells, density=0.5, format="coo", random_state=0)
    X.data = rng.integers(1, 10, size=X.nnz).astype(float)

    # scipy.io.mmwrite needs a file path; compress afterwards
    raw_mtx = folder / "matrix.mtx"
    scipy.io.mmwrite(str(raw_mtx), X)
    with open(raw_mtx, "rb") as src, gzip.open(folder / "matrix.mtx.gz", "wb") as dst:
        dst.write(src.read())
    raw_mtx.unlink()

    with gzip.open(folder / "barcodes.tsv.gz", "wt") as f:
        f.writelines(f"Cell{i}-1\n" for i in range(n_cells))

    with gzip.open(folder / "features.tsv.gz", "wt") as f:
        f.writelines(
            f"ENSG{i:06d}\tGene{i}\tGene Expression\n" for i in range(n_genes)
        )


def test_load_10x_mtx_shape(tmp_path):
    n_cells, n_genes = 4, 6
    _write_10x_mtx(tmp_path, n_cells=n_cells, n_genes=n_genes)
    adata = load_10x_mtx(tmp_path)
    assert adata.shape == (n_cells, n_genes)


def test_load_10x_mtx_var_names_are_gene_symbols(tmp_path):
    _write_10x_mtx(tmp_path)
    adata = load_10x_mtx(tmp_path)
    # Gene symbols come from the second column of features.tsv
    assert all(name.startswith("Gene") for name in adata.var_names)


def test_load_10x_mtx_missing_folder_raises(tmp_path):
    import pytest
    with pytest.raises(FileNotFoundError):
        load_10x_mtx(tmp_path / "nonexistent")


def test_load_10x_mtx_missing_files_raises(tmp_path):
    import pytest
    # Empty folder — no required files
    with pytest.raises(ValueError, match="Required file"):
        load_10x_mtx(tmp_path)


# ── Seurat RDS ─────────────────────────────────────────────────────────────────
# Not tested automatically — requires R + Seurat + rpy2.
# The ImportError path is covered by the function's docstring contract.
