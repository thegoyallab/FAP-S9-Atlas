# FAP / human S9 peptidase atlas

Reproducibility release for the manuscript:

**Evolutionary specialization and single-cell mapping reveal conserved FAP-dominant stromal remodeling across human cancers**

## What is included
- frozen scripts for the complete 1,086-dataset donor-aware human S9 atlas;
- donor-collapsed FAP/DPP4 fibroblast-state analyses;
- patient-level spatial validation;
- expanded FAP/DPP4 and S9-family evolutionary analyses;
- final targetability analyses covering organ-concordant tumor-normal enrichment, conserved cross-cancer state programs, malignant-versus-benign remodeling, abundance-versus-state coupling and cellular selectivity;
- compact source-data tables underlying main Figures 1-8 and target-biology Supplementary Figures S8-S9;
- environment/version captures and analysis-boundary documentation.

## Analysis design
1. **Family atlas:** nine human S9-family genes across all 1,086 eligible CELLxGENE Census datasets.
2. **Evolution:** family phylogeny, FAP-DPP4 topology testing, synteny, ancestral reconstruction, FAP-stem codon models and sequence-expression coupling.
3. **Fibroblast landscape:** exact-tissue donor-level FAP/DPP4 states across 114 replicated contexts.
4. **Breast and spatial validation:** within-donor breast-state contrasts plus 22-patient TNBC spatial validation.
5. **Cross-cancer program:** paired FAP-only versus DPP4-only analysis across nine organ-concordant cancer groups / 270 donors.
6. **Targetability:** matched tumor-normal enrichment, nonmalignant remodeling, cellular selectivity and an explicitly non-composite target-biology map.

## Important boundaries
Raw third-party CELLxGENE, GEO, UniProt and reference-genome datasets are not redistributed. Some frozen analysis scripts retain local project paths for audit fidelity; adjust those paths for a new environment. Transcript abundance is not treated as equivalent to surface-protein density, PET uptake, ligand binding or absorbed radioligand dose. The integrated targetability map is not a clinical ranking.

Read `docs/REPRODUCIBILITY.md`, `docs/ANALYSIS_BOUNDARIES.md`, and `docs/TARGETABILITY_ANALYSIS.md` before reuse.

## Repository
https://github.com/thegoyallab/FAP-S9-Atlas

## License status
No code/data license has been assigned in this frozen release. See `LICENSE_DECISION_REQUIRED.md`.
