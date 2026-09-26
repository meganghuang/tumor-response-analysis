"""Goal 1: cell populations that increase or decrease with treatment (Pre vs Post)."""
import scanpy as sc

from tumor_response.composition import compare_fractions, composition_by_sample, paired_fractions
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

    # All samples: baseline lesions vs on-treatment lesions (mostly different patients).
    unpaired = compare_fractions(comp, cell_types, time, "Pre", "Post")
    unpaired.to_csv(tables_dir / "treatment_composition_all_samples.csv")
    n = unpaired.attrs["n"]
    save(fraction_boxplots(comp, cell_types, time, ["Pre", "Post"], unpaired,
                           f"Cell-type fractions, baseline (n={n['Pre']}) vs on treatment (n={n['Post']})"),
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
    print(paired.round(3).to_string())


if __name__ == "__main__":
    main()
