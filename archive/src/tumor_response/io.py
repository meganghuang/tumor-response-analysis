"""Loading the Sade-Feldman et al. 2018 (Cell) melanoma scRNA-seq data, GEO GSE120575."""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

EXPRESSION_FILE = "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
METADATA_FILE = "GSE120575_patient_ID_single_cells.txt.gz"


def read_metadata(path: Path) -> pd.DataFrame:
    """Per-cell metadata indexed by cell title (e.g. "A10_P3_M11").

    The GEO file is a submission template: ~19 comment lines, then one row per cell, then
    trailing template sections. Sample IDs look like "Pre_P1", "Post_P1", or "Post_P1_2"
    (a second lesion from the same patient).
    """
    raw = pd.read_csv(path, sep="\t", skiprows=19, encoding="latin1", dtype=str)
    raw = raw[raw["Sample name"].str.match(r"Sample \d+$", na=False)]
    pid_col = next(c for c in raw.columns if c.startswith("characteristics: patinet ID"))
    meta = pd.DataFrame(
        {
            "sample": raw[pid_col].str.strip().values,
            "response": raw["characteristics: response"].str.strip().values,
            "therapy": raw["characteristics: therapy"].str.strip().values,
        },
        index=raw["title"].str.strip().values,
    )
    meta["timepoint"] = meta["sample"].str.split("_").str[0]
    meta["patient"] = meta["sample"].str.split("_").str[1]
    return meta


def read_expression(path: Path, chunksize: int = 2000) -> ad.AnnData:
    """Read the genes x cells log2(TPM+1) matrix into a sparse cells x genes AnnData.

    The first two header lines hold cell titles and sample IDs. Values are converted from
    log2(TPM+1) to ln(TPM+1), the scale scanpy assumes for log-normalized data.
    """
    with pd.io.common.get_handle(path, "r", compression="infer") as h:
        cells = h.handle.readline().rstrip("\n").split("\t")[1:]
    cells = [c.strip() for c in cells]
    n_cells = len(cells)
    while n_cells and not cells[n_cells - 1]:
        n_cells -= 1  # tolerate trailing tabs

    blocks, genes = [], []
    reader = pd.read_csv(path, sep="\t", skiprows=2, header=None, index_col=0, chunksize=chunksize)
    for chunk in reader:
        values = chunk.iloc[:, :n_cells].to_numpy(dtype=np.float32)
        np.nan_to_num(values, copy=False)
        blocks.append(sparse.csr_matrix(values.T * np.float32(np.log(2))))
        genes.extend(chunk.index.astype(str))

    adata = ad.AnnData(
        X=sparse.hstack(blocks, format="csr"),
        obs=pd.DataFrame(index=cells[:n_cells]),
        var=pd.DataFrame(index=genes),
    )
    adata.var_names_make_unique()
    return adata


def load_sade_feldman(raw_dir: Path) -> ad.AnnData:
    """Expression matrix with sample, patient, timepoint, response and therapy in .obs."""
    adata = read_expression(raw_dir / EXPRESSION_FILE)
    meta = read_metadata(raw_dir / METADATA_FILE)
    # Metadata is aligned by label, so duplicate cell titles would misalign silently.
    duplicated = adata.obs_names[adata.obs_names.duplicated()]
    if len(duplicated):
        raise ValueError(f"{len(duplicated)} duplicate cell titles, e.g. {list(duplicated[:3])}")
    if meta.index.duplicated().any():
        raise ValueError("duplicate cell titles in the metadata file")
    missing = adata.obs_names.difference(meta.index)
    if len(missing):
        raise ValueError(f"{len(missing)} cells have no metadata, e.g. {list(missing[:3])}")
    adata.obs = meta.loc[adata.obs_names].astype("category")
    return adata
