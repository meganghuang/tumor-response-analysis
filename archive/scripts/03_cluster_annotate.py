"""Cluster cells, assign cell types from known markers, and find marker genes (goal 2)."""
import scanpy as sc

from tumor_response.annotation import available_markers, cluster, find_markers, score_cell_types
from tumor_response.config import load_config, project_path
from tumor_response.plotting import set_style


def main():
    cfg = load_config()
    processed_dir = project_path(cfg["data"]["processed_dir"])
    tables_dir = project_path(cfg["output"]["tables_dir"])
    sc.settings.figdir = project_path(cfg["output"]["figures_dir"])
    set_style()

    adata = sc.read_h5ad(processed_dir / "01_preprocessed.h5ad")
    adata = cluster(adata, cfg["clustering"]["resolution"])
    adata = score_cell_types(adata)
    print(adata.obs.groupby("leiden", observed=True)["cell_type"].first().to_string())
    print(adata.obs["cell_type"].value_counts().to_string())

    d = cfg["data"]
    sc.pl.umap(adata, color=["cell_type", "leiden"], legend_loc="on data", legend_fontsize=7,
               frameon=False, save="_cell_types.png", show=False)
    sc.pl.umap(adata, color=[d["response_column"], d["timepoint_column"], d["patient_column"]],
               frameon=False, save="_samples.png", show=False)

    # Known markers: confirms the cluster labels.
    markers = available_markers(adata)
    sc.pl.dotplot(adata, markers, groupby="leiden", standard_scale="var",
                  save="_known_markers_by_cluster.png", show=False)

    # Data-driven markers: genes that distinguish each cluster / cell type from the rest.
    find_markers(adata, "leiden").to_csv(tables_dir / "markers_by_cluster.csv", index=False)
    by_type = find_markers(adata, "cell_type")
    by_type.to_csv(tables_dir / "markers_by_cell_type.csv", index=False)
    sc.pl.rank_genes_groups_dotplot(adata, key="rank_genes_cell_type", n_genes=5,
                                    standard_scale="var", save="_top_markers_by_cell_type.png",
                                    show=False)

    adata.write_h5ad(processed_dir / "02_annotated.h5ad")


if __name__ == "__main__":
    main()
