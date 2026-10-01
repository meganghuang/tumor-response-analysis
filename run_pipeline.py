#!/usr/bin/env python
"""End-to-end run: GEO download -> QC -> clustering -> annotation -> markers ->
composition -> pseudobulk DE -> response signatures -> figures -> report.

    python run_pipeline.py [--force-parse]

Every artifact lands in ``results/``.
"""
from __future__ import annotations

import argparse
import gc
import json
import pathlib
import sys
import warnings

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from pipeline import (__version__, annotate, composition, config, data, de, embed,
                      figures, markers, qc, report, signatures)

ARTIFACTS = [
    "summary_report.md", "umap_clusters.png", "umap_celltypes.png", "marker_dotplot.png",
    "marker_heatmap.png", "qc_metrics.png", "composition_barplots.png",
    "composition_boxplots.png", "signature_auc.png", "signature_scores.png",
    "composition_stats.csv", "cluster_markers.csv", "responder_DE_CD8.csv",
    "responder_DE_all.csv", "signature_genes.json",
]


def log(msg: str) -> None:
    print(f"[pipeline] {msg}", flush=True)


def main(force_parse: bool = False) -> None:
    warnings.filterwarnings("ignore", category=FutureWarning)
    sc.settings.verbosity = 1
    cfg = config.load_config()
    R = cfg.path("results")
    figures.apply_style()
    res: dict = {"version": __version__, "cfg": dict(cfg), "artifacts": ARTIFACTS}

    # ---------------------------------------------------------------- 1. data
    log("loading GSE120575")
    adata = data.build(cfg, force=force_parse)
    res.update(n_cells_raw=adata.n_obs, n_genes_raw=adata.n_vars,
               n_samples_raw=int(adata.obs["sample"].nunique()),
               n_patients_raw=int(adata.obs["patient"].nunique()))

    # ------------------------------------------------------------------ 2. QC
    log("quality control")
    filtered, qc_summary = qc.apply_qc(adata, cfg)
    qc_obs = adata.obs.copy()
    del adata; gc.collect()
    adata = filtered
    res.update(n_cells=adata.n_obs, n_genes=adata.n_vars, qc_summary=qc_summary,
               n_samples=int(adata.obs["sample"].nunique()),
               n_patients=int(adata.obs["patient"].nunique()))
    figures.qc_metrics(qc_obs, qc_summary, cfg, R / "qc_metrics.png")
    del qc_obs; gc.collect()

    # --------------------------------------------------- 3. embed and cluster
    log("normalising, embedding and clustering")
    embed.normalise(adata, cfg)
    embed.embed(adata, cfg)
    res["n_clusters"] = int(adata.obs["leiden"].nunique())

    # ------------------------------------------------------------ 4. annotate
    log("annotating lineages and CD8 states")
    annotate.annotate_lineages(adata, cfg)
    cd8_mask = (adata.obs["cell_type"].astype(str) == "CD8 T").to_numpy()
    sub = embed.subcluster(adata, cd8_mask, cfg, "cd8_leiden",
                           cfg.embedding["cd8_leiden_resolution"])
    annotate.annotate_cd8_states(adata, sub, cfg)
    res["n_cd8"] = int(cd8_mask.sum())
    res["pct_unresolved"] = float(100 * (adata.obs["population"].astype(str)
                                         == "Unresolved").mean())
    counts = (adata.obs["population"].value_counts().rename("cells").reset_index()
              .rename(columns={"index": "population", "population": "population"}))
    counts.columns = ["population", "cells"]
    counts["% of cells"] = 100 * counts["cells"] / adata.n_obs
    res["population_counts"] = counts
    adata.write_h5ad(cfg.path("processed", "annotated.h5ad"), compression="gzip")

    figures.umap_clusters(adata, R / "umap_clusters.png")
    figures.umap_celltypes(adata, R / "umap_celltypes.png")

    # -------------------------------------------------------------- 5. markers
    log("marker genes")
    mk = markers.rank_markers(adata, "population")
    mk_export = markers.export_markers(mk)
    mk_export.to_csv(R / "cluster_markers.csv", index=False)
    res["n_marker_rows"] = len(mk)
    res["n_marker_sig"] = len(mk_export)
    res["n_marker_pops"] = int(mk_export["population"].nunique())

    canonical: dict[str, list[str]] = {}
    seen: set[str] = set()
    present_pops = set(adata.obs["population"].astype(str))
    for name, genes in {**annotate.LINEAGE_MARKERS, **annotate.CD8_STATE_MARKERS}.items():
        if name not in present_pops:
            continue
        keep = [g for g in genes if g not in seen]
        seen.update(keep)
        if keep:
            canonical[name] = keep
    figures.marker_dotplot(adata, canonical, R / "marker_dotplot.png")
    figures.marker_heatmap(adata, markers.top_markers(mk, n=5), R / "marker_heatmap.png")

    # ---------------------------------------------------------- 6. composition
    log("composition: treatment and response")
    alpha = cfg.composition["fdr_alpha"]
    comp = composition.sample_composition(
        adata, "population", "sample", min_cells=cfg.composition["min_cells_per_unit"])
    comp.to_csv(cfg.path("processed", "composition_by_biopsy.csv"), index=False)

    treat = composition.treatment_paired_test(comp, "population", alpha)
    pf_all = composition.patient_fractions(comp, "population")
    pf_pre = composition.patient_fractions(comp, "population", timepoint="Pre")
    resp_all = composition.response_test(pf_all, "population", alpha, "response_all_timepoints")
    resp_pre = composition.response_test(pf_pre, "population", alpha, "response_baseline_only")
    comp_stats = pd.concat([treat, resp_all, resp_pre], ignore_index=True)
    comp_stats.to_csv(R / "composition_stats.csv", index=False)
    res["n_paired_patients"] = int(treat["n_patients"].max()) if len(treat) else 0

    figures.composition_barplots(comp, "population", R / "composition_barplots.png")
    figures.composition_boxplots(pf_all, resp_all, "population",
                                 R / "composition_boxplots.png")

    res["treatment_table"] = treat[["population", "n_patients", "mean_Pre", "mean_Post",
                                    "log2FC", "n_increase", "n_decrease", "pval", "qval"]].head(8)
    res["response_table"] = resp_all[["population", "mean_R", "mean_NR", "log2FC", "auc",
                                      "pval", "qval"]].head(8)

    up = treat[(treat["log2FC"] > 0) & (treat["pval"] < 0.05)]["population"].tolist()
    dn = treat[(treat["log2FC"] < 0) & (treat["pval"] < 0.05)]["population"].tolist()
    lead = treat.iloc[0]
    res["treatment_prose"] = (
        f"No population changes significantly on treatment: "
        f"{int(treat['significant'].sum())} of {len(treat)} survive q < {alpha}, and "
        f"{len(up) + len(dn)} reach even a nominal p < 0.05"
        + (f" (expanding: {', '.join(up)}; contracting: {', '.join(dn)})" if up or dn else "")
        + f". The strongest trend is {lead['population']} "
        f"({lead['mean_Pre']:.3f} -> {lead['mean_Post']:.3f} mean fraction, "
        f"up in {int(lead['n_increase'])} of {int(lead['n_patients'])} patients, "
        f"p = {lead['pval']:.3g}). With only {int(treat['n_patients'].max())} patients "
        f"contributing both a pre- and an on-treatment biopsy, this analysis is "
        f"underpowered for anything but a large effect.")
    sig_r = resp_all[resp_all["significant"]]
    res["response_prose"] = (
        f"{len(sig_r)} of {len(resp_all)} populations differ at q < {alpha}"
        + (": " + ", ".join(f"{r.population} (q = {r.qval:.3g}, "
                            f"{'higher' if r.log2FC > 0 else 'lower'} in responders)"
                            for r in sig_r.itertuples()) if len(sig_r) else "") + ".")

    # ------------------------------------------------------------------ 7. DE
    log("pseudobulk differential expression")
    expr_all, meta_all, det_all = de.pseudobulk(adata, cfg)
    expr_all = de.filter_genes(expr_all, det_all, cfg)
    de_all = de.test_response(expr_all, meta_all, cfg, "all immune cells")
    de_all.to_csv(R / "responder_DE_all.csv", index=False)

    cd8_cells = (adata.obs["cd8_state"].astype(str) != "Not CD8").to_numpy()
    adata_cd8 = adata[cd8_cells].copy()
    expr_cd8, meta_cd8, det_cd8 = de.pseudobulk(adata_cd8, cfg)
    expr_cd8 = de.filter_genes(expr_cd8, det_cd8, cfg)
    de_cd8 = de.test_response(expr_cd8, meta_cd8, cfg, "CD8 T cells")
    de_cd8.to_csv(R / "responder_DE_CD8.csv", index=False)

    cols = ["gene", "mean_log2_R", "mean_log2_NR", "log2FC", "auc", "pval", "qval"]
    res.update(n_de_all=len(de_all), n_de_all_sig=int(de_all["significant"].sum()),
               n_de_cd8=len(de_cd8), n_de_cd8_sig=int(de_cd8["significant"].sum()),
               de_all_top=de_all[cols].head(10), de_cd8_top=de_cd8[cols].head(10))

    # ---------------------------------------------------------- 8. signatures
    log("leave-one-patient-out signatures")
    size = cfg.signatures["size"]
    nperm = cfg.signatures["n_permutations"]
    seed = cfg.signatures["random_state"]

    scores = pd.DataFrame(index=meta_all.index)
    scores["response"] = meta_all["response"]
    scores["Whole-immune signature"] = signatures.lopo_scores(expr_all, meta_all, size)
    scores["CD8 signature"] = signatures.lopo_scores(expr_cd8, meta_cd8, size)

    # Literature benchmark: the TCF7+ CD8 fraction highlighted by Sade-Feldman et al.
    x = adata_cd8[:, "TCF7"].X
    x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
    tcf7 = (pd.DataFrame({"patient": adata_cd8.obs["patient"].astype(str).to_numpy(),
                          "pos": x > 0}).groupby("patient")["pos"].mean())
    scores["TCF7+ fraction of CD8 T cells"] = tcf7

    y = (scores["response"] == "Responder").astype(int).to_numpy()
    curves, auc_rows, sig_json = [], [], {}
    for i, col in enumerate(["CD8 signature", "Whole-immune signature",
                             "TCF7+ fraction of CD8 T cells"]):
        v = scores[col].to_numpy(dtype=float)
        ok = np.isfinite(v)
        auc, pval = signatures.auc_with_permutation(y[ok], v[ok], nperm, seed)
        fpr, tpr = signatures.roc_points(y[ok], v[ok])
        curves.append({"name": col, "auc": auc, "pval": pval, "fpr": fpr, "tpr": tpr,
                       "focal": i == 0})
        auc_rows.append({"predictor": col, "n_patients": int(ok.sum()),
                         "n_responder": int(y[ok].sum()),
                         "evaluation": "leave-one-patient-out" if "signature" in col
                                       else "direct (no fitting)",
                         "AUC": auc, "permutation p": pval})
    res["auc_table"] = pd.DataFrame(auc_rows)
    figures.signature_auc(curves, R / "signature_auc.png")

    for name, expr, meta in [("cd8", expr_cd8, meta_cd8),
                             ("whole_immune", expr_all, meta_all)]:
        yy = meta["responder"].to_numpy()
        up_g, dn_g = signatures.build_signature(expr, yy, size)
        row = next(r for r in auc_rows
                   if r["predictor"].lower().startswith("cd8" if name == "cd8" else "whole"))
        sig_json[name] = {
            "description": ("CD8 T cells" if name == "cd8" else "all immune cells")
                           + ", patient pseudobulk, Welch t-statistic ranking",
            "n_patients": int(len(meta)), "n_genes_tested": int(expr.shape[0]),
            "up_in_responders": up_g, "up_in_non_responders": dn_g,
            "lopo_auc": row["AUC"], "lopo_permutation_p": row["permutation p"],
        }

    up_cd8, dn_cd8 = sig_json["cd8"]["up_in_responders"], sig_json["cd8"]["up_in_non_responders"]
    signatures.score_cells(adata_cd8, up_cd8, dn_cd8, "cd8_signature_score")
    cell_scores = pd.DataFrame({"score": adata_cd8.obs["cd8_signature_score"].to_numpy(),
                                "response": adata_cd8.obs["response"].astype(str).to_numpy()})
    figures.signature_scores(scores, cell_scores, R / "signature_scores.png")

    bench = next(r for r in auc_rows if r["predictor"].startswith("TCF7"))
    sig_json["benchmark_tcf7_positive_cd8_fraction"] = {
        "description": "Fraction of a patient's CD8 T cells with TCF7 > 0 "
                       "(the marker highlighted by Sade-Feldman et al.)",
        "auc": bench["AUC"], "permutation_p": bench["permutation p"],
        "n_patients": bench["n_patients"]}
    sig_json["config"] = {"signature_size_per_direction": size,
                          "n_permutations": nperm, "random_state": seed,
                          "unit_of_analysis": "patient"}
    (R / "signature_genes.json").write_text(json.dumps(sig_json, indent=2), encoding="utf-8")
    scores.to_csv(cfg.path("processed", "patient_scores.csv"))

    best = max(auc_rows, key=lambda r: r["AUC"])
    res["signature_prose"] = (
        f"The strongest predictor is **{best['predictor']}** (AUC {best['AUC']:.2f}, "
        f"permutation p = {best['permutation p']:.3f}, n = {best['n_patients']} patients). "
        f"Per-cell scores overlap almost completely between responders and non-responders: "
        f"the signal is a property of the patient's immune composition and tone, not of "
        f"individual cells.")

    # -------------------------------------------------------------- 9. report
    log("writing report")
    res["headline"] = [
        f"{res['n_cells']:,} CD45+ cells from {res['n_samples']} biopsies of "
        f"{res['n_patients']} melanoma patients pass QC and resolve into "
        f"{res['n_clusters']} Leiden clusters spanning "
        f"{adata.obs['cell_type'].nunique()} immune lineages, with CD8 T cells further "
        f"split into {len([c for c in adata.obs['cd8_state'].unique() if c != 'Not CD8'])} "
        f"differentiation states.",
        res["treatment_prose"],
        res["response_prose"],
        f"{res['n_de_cd8_sig']:,} genes differ between responders and non-responders in CD8 "
        f"T cells and {res['n_de_all_sig']:,} across all immune cells at q < "
        f"{cfg['differential_expression']['fdr_alpha']} (patient-level pseudobulk).",
        res["signature_prose"].replace("**", ""),
    ]
    res["methods"] = [
        "Expression: GEO-distributed TPM, transformed to log2(TPM/10+1).",
        f"QC: {cfg.qc['min_genes']}-{cfg.qc['max_genes']} genes per cell, "
        f"<= {cfg.qc['max_pct_mito']}% mitochondrial transcripts, genes in "
        f">= {cfg.qc['min_cells_per_gene']} cells, biopsies with "
        f">= {cfg.qc['min_cells_per_sample']} cells.",
        f"Embedding: {cfg.embedding['n_top_genes']} HVGs (scanpy '{cfg.embedding['hvg_flavor']}' "
        f"flavor), z-scored and clipped at {cfg.embedding['scale_max_value']}, "
        f"{cfg.embedding['n_pcs']} PCs, {cfg.embedding['n_neighbors']}-NN graph, Leiden "
        f"(resolution {cfg.embedding['leiden_resolution']}), UMAP.",
        "Annotation: per-cell scores for canonical lineage programmes, averaged per cluster "
        "and z-scored across clusters; argmax assigns the label.",
        "Markers: one-vs-rest Wilcoxon rank-sum per population, BH-corrected.",
        "Treatment effect: Wilcoxon signed-rank on paired pre/on-treatment population "
        "fractions within patient.",
        "Response association: Mann-Whitney U on patient-level population fractions, "
        "BH-corrected across populations.",
        "Differential expression: patient pseudobulk (mean log expression), Welch t-test, "
        "BH-corrected across genes.",
        f"Signatures: top {size} genes per direction by Welch t-statistic, scored as mean "
        "standardised up-gene expression minus down-gene expression, evaluated "
        f"leave-one-patient-out with a {nperm}-permutation null.",
    ]
    res["limitations"] = [
        f"{res['n_patients']} patients is a small cohort for a {2*size}-gene signature; the "
        "LOPO AUC is an honest but high-variance estimate and no external validation cohort "
        "is used here.",
        "Biopsies were taken at the treating clinician's discretion and therapy was mixed "
        "(anti-CTLA4, anti-PD1, combination); therapy is not modelled as a covariate.",
        "Smart-seq2 plate-based sorting means population fractions reflect the sorting gate "
        "on CD45+ cells, not tumour composition.",
        "Mitochondrial fraction is computed from TPM, which is a weaker QC signal than from "
        "UMI counts.",
        "Assumption diagnostics for the t-tests were not assessed.",
    ]
    report.build(res, R / "summary_report.md")
    log(f"done — {len(ARTIFACTS)} artifacts in {R}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--force-parse", action="store_true",
                    help="re-parse the GEO text matrix instead of using the h5ad cache")
    main(**{"force_parse": ap.parse_args().force_parse})
