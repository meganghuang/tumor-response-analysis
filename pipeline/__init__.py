"""Reproducible pipeline for the Sade-Feldman melanoma checkpoint-immunotherapy atlas.

Stages
------
1. ``data``         download + parse GSE120575 into an AnnData with clinical metadata
2. ``qc``           per-cell quality control
3. ``embed``        normalisation, HVGs, PCA, kNN graph, Leiden clustering, UMAP
4. ``annotate``     lineage annotation and CD8 T-cell state assignment
5. ``markers``      per-population marker genes
6. ``composition``  which populations expand/contract with treatment and with response
7. ``de``           patient-level pseudobulk differential expression responder vs non-responder
8. ``signatures``   leave-one-patient-out signature learning and ROC evaluation
9. ``figures``      all publication figures
10. ``report``      the assembled summary_report.md
"""

__version__ = "1.0.0"
