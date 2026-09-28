"""Goal 3: signatures that distinguish responders from non-responders.

1. Cell-type composition, per patient (primary) and per biopsy (secondary, non-independent)
2. TCF7+ fraction of CD8 T cells, the key predictor reported by Sade-Feldman et al.
3. Sample-level differential expression, in all immune cells and within CD8 T cells
4. A learned gene signature, evaluated with leave-one-patient-out cross-validation

Response is annotated per biopsy in GEO, and 13 of the 32 patients gave more than one
biopsy, so biopsies are not independent observations. Every comparison is therefore reported
twice: once with the patient as the unit, and once with the biopsy as the unit. Four patients
(P1, P4, P5, P28) have biopsies whose labels disagree; they are excluded from the
patient-level comparison as `Conflicting` rather than being assigned to a group.
"""
import pandas as pd
import scanpy as sc

from tumor_response.annotation import cd8_mask
from tumor_response.composition import (
    CONFLICTING, _single_label, aggregate_to_patient, compare_fractions,
    composition_by_sample, conflicting_patients,
)
from tumor_response.config import load_config, project_path
from tumor_response.plotting import fraction_boxplots, save, score_by_group, set_style, volcano
from tumor_response.signatures import (
    auc, by_patient, fraction_positive, leave_one_patient_out_scores, pseudobulk,
    pseudobulk_matrices, pseudobulk_de,
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
    sample_col, resp_col, pat_col = d["sample_column"], d["response_column"], d["patient_column"]
    sample_info = adata.obs.groupby(sample_col, observed=True)[
        [resp_col, pat_col, d["timepoint_column"]]].first().astype(str)
    labels, patients = sample_info[resp_col], sample_info[pat_col]
    patient_labels = sample_info.groupby(pat_col)[resp_col].agg(_single_label)
    summary = []

    # 1. Composition
    comp = composition_by_sample(adata, sample_col)
    cell_types = adata.obs["cell_type"].cat.categories.astype(str).tolist()
    conflicted = conflicting_patients(comp, pat_col, resp_col)
    print(f"Patients with conflicting response labels across biopsies: {conflicted}")

    patient_comp = aggregate_to_patient(comp, cell_types, pat_col)
    n_conflict = int((patient_comp[resp_col] == CONFLICTING).sum())
    print(f"Patient-level composition: {len(patient_comp)} patients, "
          f"{n_conflict} excluded as {CONFLICTING}")

    groups = [
        ("per_patient", patient_comp, "one value per patient"),
        ("all_samples", comp, "one value per biopsy, patients repeat"),
        ("baseline", comp[comp[d["timepoint_column"]] == "Pre"], "baseline biopsies only"),
    ]
    for name, sub, note in groups:
        table = compare_fractions(sub, cell_types, resp_col, NR, R)
        table.to_csv(tables_dir / f"response_composition_{name}.csv")
        n = table.attrs["n"]
        save(fraction_boxplots(sub, cell_types, resp_col, [NR, R], table,
                               f"Cell-type fractions by response, {note} "
                               f"(NR n={n[NR]}, R n={n[R]})"),
             figures_dir / f"response_composition_{name}.png")
        print(f"{name}: NR n={n[NR]}, R n={n[R]}, "
              f"{(table['qvalue'] < 0.1).sum()} cell types at q < 0.1")

    # 2. TCF7+ CD8 T cells, per biopsy and per patient
    cd8 = cd8_mask(adata)
    print(f"CD8 T cells: {cd8.sum()}")
    tcf7 = fraction_positive(adata, "TCF7", cd8, sample_col, a["min_cells_per_sample"])
    tcf7.rename("tcf7_pos_fraction_of_cd8").to_frame().join(sample_info).to_csv(
        tables_dir / "tcf7_cd8_fraction_by_sample.csv")

    tcf7_pt = by_patient(tcf7, patients)
    keep_pt = patient_labels.reindex(tcf7_pt.index) != CONFLICTING
    tcf7_pt = tcf7_pt[keep_pt.fillna(False).to_numpy()]
    for unit, values, lab in [("patient", tcf7_pt, patient_labels), ("biopsy", tcf7, labels)]:
        value_auc = auc(values, lab, R)
        summary.append({"signature": "TCF7+ fraction of CD8 T cells (from the paper)",
                        "unit": unit, "n": len(values), "auc": value_auc,
                        "cross_validated": "not needed (no fitting)"})
        save(score_by_group(values, lab, [NR, R], "TCF7+ fraction of CD8 T cells",
                            f"TCF7+ CD8 T cells, per {unit} (AUC={value_auc:.2f})"),
             figures_dir / f"response_tcf7_cd8_per_{unit}.png")

    # 3 & 4. Differential expression and learned signatures
    for name, mask in [("all_immune", None), ("cd8", cd8)]:
        title = "CD8 T cells" if name == "cd8" else "All immune cells"
        pb = pseudobulk(adata, sample_col, mask, min_cells=a["min_cells_per_sample"],
                        min_frac_expressed=a["min_frac_expressed"])
        de = pseudobulk_de(pb, labels, NR, R)
        de.to_csv(tables_dir / f"response_de_{name}.csv")
        print(f"{name}: {len(pb)} samples, {pb.shape[1]} genes, "
              f"{(de['qvalue'] < 0.1).sum()} genes at q < 0.1")
        save(volcano(de, f"Responder vs non-responder, {title} (per-biopsy means)", NR, R),
             figures_dir / f"response_de_{name}.png")

        # Unfiltered matrices: the detection filter is refitted inside each fold, so the
        # held-out patient's cells never influence which genes are available.
        means, detect, n_cells = pseudobulk_matrices(
            adata, sample_col, mask, min_cells=a["min_cells_per_sample"])
        folds = [("all biopsies", means.index)]
        if conflicted:
            kept = means.index[~patients.loc[means.index].isin(conflicted).to_numpy()]
            folds.append((f"excluding {len(conflicted)} label-conflicted patients", kept))
        for note, idx in folds:
            scores = leave_one_patient_out_scores(
                means.loc[idx], labels.loc[idx], patients.loc[idx], NR, R, a["signature_size"],
                detect=detect, n_cells=n_cells, min_frac_expressed=a["min_frac_expressed"])
            cv_auc = auc(scores.dropna(), labels, R)
            summary.append({"signature": f"Learned {a['signature_size']}+{a['signature_size']} "
                                         f"gene signature, {title}",
                            "unit": f"biopsy, {note}", "n": int(scores.notna().sum()),
                            "auc": cv_auc, "cross_validated": "leave-one-patient-out"})
            tag = "all" if idx is means.index else "no_conflict"
            scores.rename("signature_score").to_frame().join(sample_info).to_csv(
                tables_dir / f"response_signature_scores_{name}_{tag}.csv")
            save(score_by_group(scores.dropna(), labels, [NR, R], "Signature score (held-out)",
                                f"Learned signature, {title}, {note}\n"
                                f"(cross-validated AUC={cv_auc:.2f})"),
                 figures_dir / f"response_signature_{name}_{tag}.png")

    summary = pd.DataFrame(summary)
    summary.to_csv(tables_dir / "response_signature_summary.csv", index=False)
    print(summary.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
