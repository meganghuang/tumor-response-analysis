"""Download and parse GSE120575 into an annotated AnnData object.

The GEO matrix is a 55,737 x 16,291 tab-separated TPM table with a two-row
header (cell barcode, then sample label) and a trailing empty field on every
data row. It is parsed in gene chunks and compressed to CSR on the fly, which
keeps peak memory near 1 GB instead of the ~3.6 GB a dense read would need.
"""
from __future__ import annotations

import gzip
import pathlib
import urllib.request

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

CHUNK_GENES = 4000


def download(cfg) -> tuple[pathlib.Path, pathlib.Path]:
    """Fetch the TPM matrix and the clinical metadata table from GEO."""
    raw = cfg.path("raw")
    out = []
    for key in ("tpm_file", "metadata_file"):
        name = cfg["dataset"][key]
        dst = raw / name
        if not dst.exists() or dst.stat().st_size < 1000:
            urllib.request.urlretrieve(cfg["dataset"]["base_url"] + name, dst)
        out.append(dst)
    return out[0], out[1]


def read_cell_header(tpm_path: pathlib.Path) -> tuple[list[str], list[str]]:
    """Return (cell ids, sample labels) from the two header rows."""
    with gzip.open(tpm_path, "rt") as fh:
        cells = fh.readline().rstrip("\n").split("\t")[1:]
        samples = fh.readline().rstrip("\n").split("\t")[1:]
    cells = [c for c in cells if c != ""]
    samples = samples[: len(cells)]
    return cells, samples


def parse_tpm(tpm_path: pathlib.Path, chunk_genes: int = CHUNK_GENES) -> ad.AnnData:
    """Parse the TPM text matrix into a cells x genes sparse AnnData."""
    cells, samples = read_cell_header(tpm_path)
    n_cells = len(cells)
    names = ["gene"] + cells + ["_trailing"]

    blocks, genes = [], []
    reader = pd.read_csv(
        tpm_path,
        sep="\t",
        skiprows=2,
        header=None,
        names=names,
        usecols=["gene"] + cells,
        index_col=0,
        dtype={c: np.float32 for c in cells},
        chunksize=chunk_genes,
        compression="gzip",
    )
    for chunk in reader:
        genes.extend(chunk.index.astype(str).tolist())
        blocks.append(sp.csr_matrix(chunk.to_numpy(dtype=np.float32)))

    matrix = sp.vstack(blocks, format="csr").T.tocsr()   # -> cells x genes
    adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame({"sample_label": samples}, index=pd.Index(cells, name="cell")),
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    adata.var_names_make_unique()
    return adata


def parse_metadata(meta_path: pathlib.Path) -> pd.DataFrame:
    """Parse the GEO SAMPLES block into one row per cell."""
    # The GEO metadata template is latin-1 encoded (it contains a micro sign).
    with gzip.open(meta_path, "rt", encoding="latin-1") as fh:
        lines = fh.read().splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("Sample name\t"))
    rows = []
    for line in lines[start + 1:]:
        if not line.startswith("Sample "):
            break
        parts = line.split("\t")
        rows.append(parts[:7])
    meta = pd.DataFrame(
        rows,
        columns=["geo_sample", "cell", "source", "organism", "sample_label", "response", "therapy"],
    ).set_index("cell")
    meta["timepoint"] = meta["sample_label"].str.split("_").str[0].str.capitalize()
    meta["patient"] = meta["sample_label"].str.extract(r"_(P\d+)")[0]
    # Two biopsies in the series carry an enrichment suffix; keep it in the sample id
    # so the biopsy stays the unit, but strip it from the patient id.
    meta["sample"] = meta["sample_label"]
    meta["response"] = meta["response"].str.strip()
    meta["therapy"] = meta["therapy"].str.strip()
    return meta[["sample", "patient", "timepoint", "response", "therapy", "geo_sample"]]


def build(cfg, force: bool = False) -> ad.AnnData:
    """Raw-TPM AnnData with clinical metadata attached.

    The text-matrix parse is the only slow step (~18 min) and is cached on its
    own, so a change to the metadata join never costs a re-parse.
    """
    tpm_path, meta_path = download(cfg)

    # Metadata first: it is cheap, and parsing it before the matrix means a
    # malformed metadata file fails in seconds rather than after the parse.
    meta = parse_metadata(meta_path)

    cache = cfg.path("processed", "tpm_matrix.h5ad")
    if cache.exists() and not force:
        adata = ad.read_h5ad(cache)
    else:
        adata = parse_tpm(tpm_path)
        adata.write_h5ad(cache, compression="gzip")

    missing = adata.obs_names.difference(meta.index)
    if len(missing):
        raise ValueError(f"{len(missing)} cells in the matrix have no metadata row")
    adata.obs = adata.obs.join(meta, how="left")
    adata.obs["responder"] = (adata.obs["response"] == "Responder").astype(int)
    for col in ["sample", "patient", "timepoint", "response", "therapy"]:
        adata.obs[col] = adata.obs[col].astype("category")
    return adata
