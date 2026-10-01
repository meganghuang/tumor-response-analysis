"""Goal 1: cell populations that increase or decrease with treatment (Pre vs Post)."""
import scanpy as sc

from tumor_response.composition import (
    compare_fractions, composition_by_sample, paired_fractions, single_timepoint_samples,
    therapy_switchers,
)
from tumor_response.config import load_config, project_path
from tumor_response.plotting import fraction_boxplots, paired_slopes, save, set_style


def main():
    cfg = load_config()
    d = cfg["data"]
    processed_dir = project_path(d["processed_dir"])
    figures_dir = project_path(cfg["output"]["figures_dir"])
    tables_dir = project_path(cfg["output"]["tables_dir"])
    set_style()

    adata = sc.read_h5ad(processed_dir / "02_annotated.h5ad")
    comp = composition_by_sample(adata, d["sample_column"])
    comp.to_csv(tables_dir / "composition_by_sample.csv")
    cell_types = adata.obs["cell_type"].cat.categories.astype(str).tolist()
    time = d["timepoint_column"]

    # Primary unpaired test: only patients biopsied at a single timepoint, so the two groups
    # share no patients. Mann-Whitney U compares independent groups, and 11 patients appear
    # at both timepoints; including them would make the groups overlap.
    disjoint = single_timepoint_samples(comp, d["patient_column"], time)
    independent = compare_fractions(disjoint, cell_types, time, "Pre", "Post")
    independent.to_csv(tables_dir / "treatment_composition_independent_patients.csv")
    n = independent.attrs["n"]
    save(fraction_boxplots(disjoint, cell_types, time, ["Pre", "Post"], independent,
                           f"Cell-type fractions, patients biopsied once: baseline (n={n['Pre']}) "
                           f"vs on treatment (n={n['Post']})"),
         figures_dir / "treatment_composition_independent_patients.png")

    # Kept for comparison with the original run: every sample, ignoring that 11 patients
    # contribute to both groups. Reported as a non-independent secondary result.
    unpaired = compare_fractions(comp, cell_types, time, "Pre", "Post")
    unpaired.to_csv(tables_dir / "treatment_composition_all_samples.csv")
    n_all = unpaired.attrs["n"]
    save(fraction_boxplots(comp, cell_types, time, ["Pre", "Post"], unpaired,
                           f"Cell-type fractions, all samples, groups share 11 patients: "
                           f"baseline (n={n_all['Pre']}) vs on treatment (n={n_all['Post']})"),
         figures_dir / "treatment_composition_all_samples.png")

    # Same comparison within responders and within non-responders: a population can move in
    # opposite directions depending on outcome.
    for resp, sub in comp.groupby(d["response_column"]):
        table = compare_fractions(sub, cell_types, time, "Pre", "Post")
        table.to_csv(tables_dir / f"treatment_composition_{resp.lower()}s.csv")

    # Patients sampled both before and during treatment: the most direct test.
    paired, long = paired_fractions(comp, cell_types, d["patient_column"], time)
    paired.to_csv(tables_dir / "treatment_composition_paired.csv")
    save(paired_slopes(long, cell_types, paired), figures_dir / "treatment_composition_paired.png")

    print(f"Paired patients: {paired.attrs['n_patients']}")
    switchers = therapy_switchers(comp, d["patient_column"])
    if switchers:
        print("Patients whose biopsies span more than one therapy, so their Pre/Post pair is "
              f"not one treatment course: {switchers}")
    print(paired.round(3).to_string())


if __name__ == "__main__":
    main()
