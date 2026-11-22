[Reading 33 lines from start (total: 33 lines, 0 remaining)]

# FAP / human S9 peptidase atlas

Reproducibility release for the manuscript:

**Single-cell, spatial and evolutionary profiling of human S9 peptidases reveals FAP-dominant fibroblast remodeling in breast cancer**

## What is included
- frozen analysis scripts used for the final donor-collapsed single-cell analysis;
- frozen patient-level spatial validation scripts;
- final evolutionary figure-generation scripts;
- source-data tables supporting the principal numerical claims and main figures;
- vector previews of Figures 1-5;
- environment/version captures and analysis-boundary documentation.

## What is not included
Raw CELLxGENE, GEO, UniProt, or reference-genome data are not redistributed. Obtain them from their original public repositories using the accessions listed in `docs/DATA_ACCESSIONS.md`.

## Main analysis design
1. **Family atlas:** nine human S9-family genes across the complete 1,086-dataset CELLxGENE Census extraction.
2. **Fibroblast analysis:** exact-tissue donor-level FAP/DPP4 states after assay/suspension collapse; 114 replicated contexts.
3. **Breast cancer:** within-donor FAP-only vs DPP4-only module contrasts, replicated across three collections.
4. **Spatial validation:** 22 TNBC patients / 43 Visium sections; patient is the inferential unit.
5. **Evolution:** phylogeny, topology testing, synteny, ancestral reconstruction, codon models and sequence-expression coupling.

## Reproducibility cautions
Read `docs/REPRODUCIBILITY.md` and `docs/ANALYSIS_BOUNDARIES.md` before reuse. Frozen scripts retain some local project paths for audit fidelity; update those paths for your environment.

## Release status
This repository is the public reproducibility release accompanying the manuscript. The final manuscript author list and code/data license remain to be confirmed by the authors before archival DOI minting.


## Repository
https://github.com/thegoyallab/FAP-S9-Atlas

[executed on device: bioline (c855b460-d935-487c-b56c-083cb9e61e47)]