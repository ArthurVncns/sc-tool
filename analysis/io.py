import tempfile
from pathlib import Path
from typing import BinaryIO

import anndata as ad


def load_h5ad(source: BinaryIO) -> ad.AnnData:
    """Load an AnnData object from a binary file-like source.

    Args:
        source: Any binary-readable object (e.g. open file handle, Streamlit UploadedFile).

    Returns:
        Loaded AnnData object.

    Note:
        HDF5 (the format underlying .h5ad) requires seekable random-access IO.
        anndata's reader expects a file path, not a file-like object, so we
        write to a temporary file, read from it, then delete it immediately.
    """
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp.write(source.read())
        tmp_path = Path(tmp.name)

    try:
        return ad.read_h5ad(tmp_path)
    finally:
        tmp_path.unlink(missing_ok=True)
