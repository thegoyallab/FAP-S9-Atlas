# Reproducibility notes

## Statistical unit
Donors/patients are the biological replicates. Cells and Visium spots are intermediate observations and are not treated as independent biological replicates.

## Single-cell analysis
The primary fibroblast analysis uses exact CELLxGENE `tissue`, collection-scoped donor IDs, and donor collapse across assay/suspension strata before inference. Final inference files correspond to the donor-collapsed `31d/35c/33c/33d/38b` analysis lineage.

## Spatial analysis
The primary spatial cohort is GSE210616. Repeated sections are collapsed within patient. FAP is removed from the fibroblast-affinity signature before FAP-specific localization tests to avoid circularity.

## Evolution
Final phylogenetic, ASR and codon-model results are frozen. The corrected DPP6/DPP10 ancestral reconstruction favors Ser at the vertebrate ancestral node (posterior probability 0.80255); the precise intermediate substitution path is unresolved.

## Paths
The frozen scripts preserve some project-local paths for audit fidelity. For public reuse, set a local project root and update input paths to the downloaded public datasets. The release contains source tables sufficient to verify the reported main-figure numerical summaries without redistributing raw third-party data.

## Targetability extension (v2)
The consolidated Nature Communications v8 analysis adds organ-concordant tumor-normal FAP enrichment, cross-cancer paired FAP-only/DPP4-only state contrasts, same-study malignant-versus-benign lung remodeling, abundance-versus-state coupling, cellular selectivity and an explicitly non-composite integrated targetability map. Final scripts are under `scripts/targetability/`; compact figure-level source data are under `source_data/targetability/`. See `docs/TARGETABILITY_ANALYSIS.md`.
