"""Goal 3: signatures that distinguish responders from non-responders.

1. Cell-type composition (all samples and baseline-only)
2. TCF7+ fraction of CD8 T cells, the key predictor reported by Sade-Feldman et al.
3. Sample-level differential expression, in all immune cells and within CD8 T cells
4. A learned gene signature, evaluated with leave-one-patient-out cross-validation
"""
import pandas as pd
import scanpy as sc

from tumor_response.annotation import cd8_mask
from tumor_response.composition import compare_fractions, composition_by_sample
from tumor_response.config import load_config, project_path
from tumor_response.plotting import fraction_boxplots, save, score_by_group, set_style, volcano
from tumor_response.signatures import (
    auc, fraction_positive, leave_one_patient_out_scores, pseudobulk, pseudobulk_de,
)

NR, R = "Non-responder", "Responder"


def main():
    cfg = load_config()
    d, a = cfg["data"], cfg["analysis"]
    processed_dir = project_path(d["processed_dir"])
    figures_dir = project_path(cfg["output"]["figures_dir"])
    tables_dir = project_path(cfg["output"]["tables_dir"])
    set_style()

    adata = sc.read_h5ad(processed_dir / "02_annotated.h5ad")
    sample_col, resp_col = d["sample_column"], d["response_column"]
    sample_info = adata.obs.groupby(sample_col, observed=True)[
        [resp_col, d["patient_column"], d["timepoint_column"]]].first().astype(str)
    labels, patients = sample_info[resp_col], sample_info[d["patient_column"]]
    summary = []

    # 1. Composition
    comp = composition_by_sample(adata, sample_col)
    cell_types = adata.obs["cell_type"].cat.categories.astype(str).tolist()
    for name, sub in [("all_samples", comp), ("baseline", comp[comp[d["timepoint_column"]] == "Pre"])]:
        table = compare_fractions(sub, cell_types, resp_col, NR, R)
        table.to_csv(tables_dir / f"response_composition_{name}.csv")
        n = table.attrs["n"]
        save(fraction_boxplots(sub, cell_types, resp_col, [NR, R], table,
                               f"Cell-type fractions by response, {name.replace('_', ' ')} "
                               f"(NR n={n[NR]}, R n={n[R]})"),
             figures_dir / f"response_composition_{name}.png")

    # 2. TCF7+ CD8 T cells
    cd8 = cd8_mask(adata)
    print(f"CD8 T cells: {cd8.sum()}")
    tcf7 = fraction_positive(adata, "TCF7", cd8, sample_col, a["min_cells_per_sample"])
    tcf7_auc = auc(tcf7, labels, R)
    summary.append({"signature": "TCF7+ fraction of CD8 T cells (from the paper)",
                    "n_samples": len(tcf7), "auc": tcf7_auc, "cross_validated": "not needed (no fitting)"})
    tcf7.rename("tcf7_pos_fraction_of_cd8").to_frame().join(sample_info).to_csv(
        tables_dir / "tcf7_cd8_fraction_by_sample.csv")
    save(score_by_group(tcf7, labels, [NR, R], "TCF7+ fraction of CD8 T cells",
                        f"TCF7+ CD8 T cells (AUC={tcf7_auc:.2f})"),
         figures_dir / "response_tcf7_cd8.png")

    # 3 & 4. Differential expression and learned signatures
    for name, mask in [("all_immune", None), ("cd8", cd8)]:
        pb = pseudobulk(adata, sample_col, mask, min_cells=a["min_cells_per_sample"])
        de = pseudobulk_de(pb, labels, NR, R)
        de.to_csv(tables_dir / f"response_de_{name}.csv")
        print(f"{name}: {len(pb)} samples, {pb.shape[1]} genes, "
              f"{(de['qvalue'] < 0.1).sum()} genes at q < 0.1")
        title = "CD8 T cells" if name == "cd8" else "All immune cells"
        save(volcano(de, f"Responder vs non-responder, {title} (per-sample means)", NR, R),
             figures_dir / f"response_de_{name}.png")

        scores = leave_one_patient_out_scores(pb, labels, patients, NR, R, a["signature_size"])
        cv_auc = auc(scores.dropna(), labels, R)
        summary.append({"signature": f"Learned {a['signature_size']}+{a['signature_size']} gene "
                                     f"signature, {title}",
                        "n_samples": int(scores.notna().sum()), "auc": cv_auc,
                        "cross_validated": "leave-one-patient-out"})
        scores.rename("signature_score").to_frame().join(sample_info).to_csv(
            tables_dir / f"response_signature_scores_{name}.csv")
        save(score_by_group(scores.dropna(), labels, [NR, R], "Signature score (held-out)",
                            f"Learned signature, {title}\n(cross-validated AUC={cv_auc:.2f})"),
             figures_dir / f"response_signature_{name}.png")

    summary = pd.DataFrame(summary)
    summary.to_csv(tables_dir / "response_signature_summary.csv", index=False)
    print(summary.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
