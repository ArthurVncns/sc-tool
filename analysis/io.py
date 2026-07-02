import tempfile
from pathlib import Path
from typing import BinaryIO

import anndata as ad
import scanpy as sc


def load_h5ad(source: BinaryIO) -> ad.AnnData:
    """Load an AnnData object from an .h5ad file-like source.

    HDF5 requires seekable random-access IO, so we write to a temporary file,
    read from it, then delete it immediately.
    """
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp.write(source.read())
        tmp_path = Path(tmp.name)
    try:
        return ad.read_h5ad(tmp_path)
    finally:
        tmp_path.unlink(missing_ok=True)


def load_h5(source: BinaryIO) -> ad.AnnData:
    """Load a 10x Genomics HDF5 (.h5) file produced by Cell Ranger.

    Args:
        source: Binary-readable object pointing to the .h5 file.

    Returns:
        AnnData with raw counts in X, genes as var, barcodes as obs.
    """
    with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as tmp:
        tmp.write(source.read())
        tmp_path = Path(tmp.name)
    try:
        return sc.read_10x_h5(tmp_path)
    finally:
        tmp_path.unlink(missing_ok=True)


def load_10x_mtx(folder: Path) -> ad.AnnData:
    """Load 10x Genomics sparse matrix format from a local folder.

    The folder must contain:
        matrix.mtx.gz   (or matrix.mtx)
        barcodes.tsv.gz (or barcodes.tsv)
        features.tsv.gz (or genes.tsv.gz for CellRanger 2.x)

    Args:
        folder: Path to the directory containing the three files.

    Returns:
        AnnData with raw counts in X, gene symbols as var_names.

    Raises:
        FileNotFoundError: If the folder does not exist.
        ValueError: If required files are missing from the folder.
    """
    folder = Path(folder)
    if not folder.exists():
        raise FileNotFoundError(f"Folder not found: {folder}")

    # Validate that the minimum required files are present (compressed or not)
    required = ["matrix.mtx", "barcodes.tsv"]
    feature_candidates = ["features.tsv", "genes.tsv"]  # CellRanger 3 vs 2

    for name in required:
        if not (folder / name).exists() and not (folder / f"{name}.gz").exists():
            raise ValueError(
                f"Required file '{name}' (or '{name}.gz') not found in {folder}."
            )

    has_features = any(
        (folder / name).exists() or (folder / f"{name}.gz").exists()
        for name in feature_candidates
    )
    if not has_features:
        raise ValueError(
            f"Feature file (features.tsv or genes.tsv) not found in {folder}."
        )

    return sc.read_10x_mtx(folder, var_names="gene_symbols", cache=False)


def load_seurat_rds(source: BinaryIO) -> ad.AnnData:
    """Convert a Seurat RDS file to AnnData via rpy2.

    Extracts raw counts, cell metadata, gene names, and any stored
    PCA / UMAP embeddings from the Seurat object.

    Args:
        source: Binary-readable object pointing to the .rds file.

    Returns:
        AnnData with raw counts in X (cells × genes).

    Raises:
        ImportError: If rpy2 is not installed.
        ValueError: If the RDS file does not contain a Seurat object,
            or if the Seurat or Matrix R packages are unavailable.
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
    except ImportError as exc:
        raise ImportError(
            "Loading Seurat RDS files requires rpy2 and an R installation.\n"
            "Install rpy2 with:  pip install rpy2\n"
            "Then install Seurat in R:  install.packages('Seurat')"
        ) from exc

    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    with tempfile.NamedTemporaryFile(suffix=".rds", delete=False) as tmp:
        tmp.write(source.read())
        tmp_path = Path(tmp.name)

    try:
        # Load the RDS and validate the object type
        ro.r(f'seurat_obj <- readRDS("{tmp_path}")')

        is_seurat = bool(ro.r('inherits(seurat_obj, "Seurat")')[0])
        if not is_seurat:
            raise ValueError(
                "The RDS file does not contain a Seurat object. "
                "sc_tool only supports Seurat RDS conversion."
            )

        try:
            ro.r("suppressPackageStartupMessages(library(Seurat))")
            ro.r("suppressPackageStartupMessages(library(Matrix))")
        except Exception as exc:
            raise ValueError(
                "The Seurat and Matrix R packages must be installed.\n"
                "In R, run:  install.packages(c('Seurat', 'Matrix'))"
            ) from exc

        # ── Cell and gene names ──────────────────────────────────────────
        cell_names = [str(x) for x in ro.r("colnames(seurat_obj)")]
        gene_names = [str(x) for x in ro.r("rownames(seurat_obj)")]
        default_assay = str(ro.r("DefaultAssay(seurat_obj)")[0])

        # ── Count matrix (genes × cells → cells × genes) ────────────────
        # Prefer raw counts; fall back to normalized data if counts slot is empty.
        ro.r(
            f'mat <- GetAssayData(seurat_obj, assay="{default_assay}", slot="counts")'
        )
        is_empty = bool(ro.r("prod(dim(mat)) == 0 || nnzero(mat) == 0")[0])
        if is_empty:
            ro.r(
                f'mat <- GetAssayData(seurat_obj, assay="{default_assay}", slot="data")'
            )

        nrow = int(ro.r("nrow(mat)")[0])  # genes
        ncol = int(ro.r("ncol(mat)")[0])  # cells

        # Extract the CSC sparse matrix components from R's dgCMatrix
        i_arr = np.asarray(ro.r("mat@i"), dtype=np.int32)
        p_arr = np.asarray(ro.r("mat@p"), dtype=np.int32)
        x_arr = np.asarray(ro.r("mat@x"), dtype=np.float32)

        counts = sp.csc_matrix((x_arr, i_arr, p_arr), shape=(nrow, ncol))
        counts = counts.T.tocsr()  # cells × genes

        # ── Cell metadata ────────────────────────────────────────────────
        meta_r = ro.r("seurat_obj@meta.data")
        with (ro.default_converter + pandas2ri.converter).context():
            meta_df = ro.conversion.rpy2py(meta_r)
        meta_df.index = cell_names

        # ── Build AnnData ────────────────────────────────────────────────
        adata = ad.AnnData(
            X=counts,
            obs=meta_df,
            var=pd.DataFrame(index=gene_names),
        )
        adata.obs_names = cell_names
        adata.var_names = gene_names

        # ── Embeddings (optional) ────────────────────────────────────────
        reductions = [str(x) for x in ro.r("names(seurat_obj@reductions)")]

        if "pca" in reductions:
            pca = np.asarray(ro.r("seurat_obj@reductions$pca@cell.embeddings"))
            adata.obsm["X_pca"] = pca

        if "umap" in reductions:
            umap = np.asarray(ro.r("seurat_obj@reductions$umap@cell.embeddings"))
            adata.obsm["X_umap"] = umap

        return adata

    finally:
        tmp_path.unlink(missing_ok=True)
