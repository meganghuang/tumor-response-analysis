# Audit of the scanpy pipeline

Read-through of `src/tumor_response/` and `scripts/` at commit `58b7ed9`, with the GEO
files checked directly to confirm what the code actually receives. Findings are split into
defects that change reported numbers (A), methodological weaknesses that should be stated
rather than silently carried (B), and robustness gaps (C). Section D lists what was checked
and found correct, so the audit is not mistaken for a blanket criticism.

## A. Defects that change reported numbers

### A1. Mitochondrial fraction is computed from log-scale values

`qc.annotate_qc_metrics` calls `sc.pp.calculate_qc_metrics` on `adata.X`, which by that
point holds `ln(TPM+1)` (converted in `io.read_expression`). `calculate_qc_metrics` sums
`X`, so `total_counts` is a sum of logarithms and

    pct_counts_mt = 100 * sum(ln(TPM_mito + 1)) / sum(ln(TPM_all + 1))

which is not a fraction of transcripts. The logarithm compresses large values, so the
handful of very highly expressed mitochondrial genes contribute far less to the numerator
than their share of transcripts warrants: the reported percentage is biased downward for
exactly the high-mitochondrial cells the `max_pct_mito: 20` filter exists to remove. The
filter is therefore more permissive than intended.

**Fix.** Compute the fraction on the back-transformed TPM scale, `expm1(X)`, before any
scaling. Also confirm at runtime that `MT-` prefixed genes exist in this annotation at all —
if they do not, the filter is a silent no-op and should be reported as such rather than
appearing to have been applied.

### A2. Forty-eight samples from thirty-two patients are tested as independent

The response analyses in `05_response_signatures.py` use the sample as the unit: 19 patients
contribute one sample, 10 contribute two, and 3 contribute three. The effective number of
independent units is 32, not 48.

This affects `composition.compare_fractions` (the `q = 0.002` for memory/naive T cells and
`q = 0.01` for B cells), `signatures.pseudobulk_de`, and `signatures.auc` — including the
headline TCF7+ CD8 figure of AUC 0.72, which is computed over samples. Mann-Whitney U and
the Benjamini-Hochberg q-values derived from it assume independent observations, so the
p-values are anticonservative and the AUC's effective sample size is overstated.

**Fix.** Make a patient-level analysis (one value per patient) the primary result and report
the sample-level numbers alongside it, so the difference is visible rather than hidden.

### A3. Four patients carry contradictory response labels

Response is annotated per biopsy in GEO, and for four patients the lesions disagree:

| patient | sample | timepoint | response | therapy |
|---|---|---|---|---|
| P1 | Pre_P1 | Pre | Responder | anti-CTLA4 |
| P1 | Post_P1 | Post | Responder | anti-PD1 |
| P1 | Post_P1_2 | Post | Non-responder | anti-PD1 |
| P4 | Pre_P4 | Pre | Non-responder | anti-CTLA4+PD1 |
| P4 | Post_P4 | Post | Responder | anti-CTLA4+PD1 |
| P5 | Post_P5 | Post | Non-responder | anti-PD1 |
| P5 | Post_P5_2 | Post | Responder | anti-PD1 |
| P28 | Pre_P28 | Pre | Responder | anti-CTLA4+PD1 |
| P28 | Post_P28 | Post | Non-responder | anti-CTLA4+PD1 |
| P28 | Post_P28_2 | Post | Non-responder | anti-CTLA4+PD1 |

P1 and P5 have two lesions at the *same* timepoint with opposite labels. The consequence is
that these patients appear in both arms of every responder-versus-non-responder comparison,
and in the cross-validation a held-out patient can be simultaneously its own positive and
negative case. Nothing in the code detects this.

P1 also switches therapy between its baseline and on-treatment biopsy (anti-CTLA4 to
anti-PD1), so its "paired" Pre/Post pair does not describe one treatment course.

**Fix.** Choose a rule and state it in the output rather than letting the conflict pass
silently: the patient-level primary analysis excludes the four conflicted patients, and the
sample-level (lesion-level) analysis keeps them but is reported as lesion-level with the
conflict named.

### A4. The unpaired Pre-versus-Post test mixes in eleven paired patients

`04_treatment_changes.py` runs `compare_fractions(comp, ..., "Pre", "Post")` across all
samples, but 11 patients (P1, P2, P3, P4, P6, P7, P8, P12, P15, P20, P28) contribute to
both groups. Mann-Whitney U compares two independent groups; here the groups share
patients. The README notes the design but the test is not adjusted.

**Fix.** Restrict the unpaired test to the 8 Pre-only and 13 Post-only patients, and treat
the paired Wilcoxon on the 11 both-timepoint patients as the primary timepoint result.

### A5. The expression gene filter sits outside the cross-validation fold

`05_response_signatures.py` calls `pseudobulk(..., min_frac_expressed=0.1)` once on all
samples, then passes the resulting matrix to `leave_one_patient_out_scores`. Gene selection
and z-scoring are correctly restricted to training samples inside the loop, but the *gene
universe* was chosen with the held-out patient's cells included.

The leakage is unsupervised — no response labels enter the filter — so the inflation of the
0.85 and 0.89 cross-validated AUCs should be small. It is also trivial to remove, and a
cross-validation claim is worth stating without an asterisk.

**Fix.** Apply the detection filter inside each fold, on training samples only.

## B. Methodological weaknesses worth stating

### B1. Cross-validated scores are pooled across folds

Each fold produces a different gene set and a different z-score reference, so scores from
different folds are not on a common scale; pooling them into a single AUC assumes they are
comparable. This is common practice, but a per-fold AUC summary or within-fold rank
normalisation is a cheap sensitivity check.

### B2. Compositional data is tested one cell type at a time

Cell-type fractions sum to one, so an increase in one population mechanically depresses
every other, and each per-cell-type test ignores that constraint. `Unresolved` sits in the
denominator, so a change in how many cells fail to be labelled shifts every reported
fraction. A centred log-ratio transform, or a Dirichlet-multinomial model, would respect the
constraint.

### B3. Marker p-values come from a cell-level test

`annotation.find_markers` runs a Wilcoxon test over thousands of individual cells, which
contradicts the per-sample principle the rest of the repo follows and produces p-values
driven by cell count rather than by the number of independent patients. The gene *ranking*
is still useful; the p-values should be labelled descriptive.

### B4. An upper gene-count filter is not a doublet filter here

`max_genes: 10000` is a droplet-data idiom. These are plate-based Smart-seq2 libraries with
one cell per well, so a high detected-gene count is not doublet evidence. How many cells the
ceiling removes should be reported before it is kept.

### B5. Cell-type labels mix lineage and state

`Cycling` and `Memory-naive T` are states that span CD4 and CD8 lineages, so the composition
comparisons mix "which lineages are present" with "what state the T cells are in". This is a
defensible choice for replicating the paper's populations, but it makes the fractions harder
to interpret as lineage abundance and should be said out loud.

## C. Robustness gaps

- **C1.** `annotation.cd8_mask` indexes `adata.raw[:, ["CD8A", "CD8B"]]` and raises
  `KeyError` if either gene is absent from the annotation.
- **C2.** `io.load_sade_feldman` aligns metadata with `meta.loc[adata.obs_names]`; duplicate
  cell titles would misalign silently. Uniqueness should be asserted.
- **C3.** `signatures.top_genes` can return fewer than `n` genes without saying so.

## D. Checked and correct

These were verified against the GEO files rather than assumed:

- **Unit conversion.** `log2(TPM+1) * ln(2) = ln(TPM+1)`. Correct, and it matters, because
  scanpy's downstream defaults assume natural-log data.
- **Column alignment.** The expression header has 16,292 fields (a leading empty cell plus
  16,291 cell titles) while data rows have 16,293 fields, the last being empty from a
  trailing tab. Truncating to `n_cells = 16291` is correct, and a one-column shift — which
  would have silently corrupted every downstream result — does not occur. All 16,291 cells
  have matching metadata.
- **Metadata offset.** `skiprows=19` lands exactly on the `Sample name` header row, and the
  `Sample \d+$` filter correctly strips the template's trailing sections.
- **AUC direction.** `mannwhitneyu(positive, negative).statistic / (n_pos * n_neg)` is
  `P(positive > negative)`, as documented, in both `compare_fractions` and `signatures.auc`.
- **Normalisation order.** `adata.raw` is set after normalisation and before scaling, and
  `score_genes`, `rank_genes_groups`, `pseudobulk` and `cd8_mask` all read `use_raw=True`,
  so no analysis is run on z-scored values.
- **Cross-validation z-reference.** `signature_score(..., ref=pb[train])` uses training
  samples only. Gene selection is likewise inside the fold. Only the detection filter (A5)
  escapes.
- **Multiple testing.** Benjamini-Hochberg via `scipy.stats.false_discovery_control`.
- **Cross-implementation agreement on the load and QC.** An independent R implementation
  reading the same archives arrives at the same 16,291 cells, the same 37 mitochondrial
  genes, the same 13.94% mean mitochondrial fraction, and the same 1,786 cells removed,
  leaving the same 14,505 x 35,491 matrix. The corrected mitochondrial calculation was also
  checked against the TPM constraint directly: per-cell sums of `expm1(X)` have a median of
  994,726, so the values really are transcripts per million and a 12.7% median mitochondrial
  content is the true quantity.

## E. Found while running, not in the read-through

- **`sc.pp.neighbors` cannot use its default backend here.** pynndescent builds a joblib
  thread pool, which needs a named pipe the sandbox denies, so the step died after QC.
  `preprocess.normalize_and_reduce` now passes `transformer="sklearn"`, an exact
  nearest-neighbour search that is also fast at this scale and removes a source of
  run-to-run drift.
- **`min_genes`/`max_genes` are inert on this dataset.** Detected genes run from 1,093 to
  9,315, so neither bound removes a cell. Predicted in B4 and confirmed on the data.
- **The `Unresolved` label is doing more work than intended.** 45.9% of cells fall in it,
  and the parallel Seurat implementation resolves 3,322 of them (23% of all cells) as CD8
  exhausted. The margin rule is discarding the paper's central population, which is an
  analysis-design problem rather than a coding defect, but it dominates the interpretation
  of every composition table.
