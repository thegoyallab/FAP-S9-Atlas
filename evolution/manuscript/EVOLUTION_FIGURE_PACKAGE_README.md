# Final evolution figure package — complete 1086-dataset Census analysis

## Main panel mapping
A. `Deep_S9_curated_ML_publication_FINAL` — frozen 36-sequence, 1,192-aa curated deep-S9 ML tree; IQ-TREE model Q.PFAM+R4; source lnL -52433.7444. Support labels are displayed only at SH-aLRT >=90 and UFBoot >=95.
B. `Panel_B_S9_duplication_history_FINAL` — validated duplication-history schematic and topology-test summary.
C. `Panel_C_FAP_DPP4_synteny_FINAL` — validated sarcopterygian FAP-DPP4 locus architecture and intergenic distances.
D. `Panel_D_FAP_A657_ASR_FINAL` — FAP A657-equivalent ancestral-state posterior probabilities.
E. `Panel_E_FAP_branchsite_structure_FINAL` — PyMOL FAP catalytic-domain render combined with best multi-start branch-site/BEB statistics. This panel contains a raster-backed structural image embedded within vector output.
F. `Panel_F_DPP6_DPP10_ASR_FINAL` — corrected DPP6/DPP10 catalytic-site ancestral-state posterior probabilities. It intentionally does not use the superseded S->A->D/G interpretation.
G. `S9_sister_pair_expression_divergence_final1086` — final 1,086-dataset evolution-expression integration panel.

## Final expression statistics
- 1,431 robust tissue x cell-type strata across all nine genes and 44 tissues.
- Global sequence-distance vs expression-distance Spearman rho = -0.0849421; 9,999-permutation P = 0.7076.
- FAP-DPP4 rho = 0.135404 (95% bootstrap interval 0.074012-0.210391).
- DPP8-DPP9 rho = 0.630923 (0.570331-0.695833).
- DPP6-DPP10 rho = 0.784677 (0.610035-0.873450).
- PREP-PREPL rho = 0.578557 (0.451548-0.678253).
- FAP-DPP4 minus mean(other three sister pairs) = -0.529315; 95% interval -0.596699 to -0.421627; 5,000/5,000 bootstrap contrasts below zero.

## Interpretation safeguards
- Do not infer causal coupling between sequence divergence and expression divergence.
- Do not claim the teleost FAP/DPP4 duplication specifically derives from 3R from current synteny evidence.
- Distinguish sequence-likelihood topology evidence from duplication-loss parsimony.
- FAP A657 is a functional/ASR specificity determinant and is not a BEB-selected codon.
- DPP6/DPP10: ancestral Ser is favored at the vertebrate node; crown DPP6/DPP10 are D/G; the exact intermediate path is unresolved.

Checkpoint1000 files remain outside this package in their preserved project directories and were not overwritten.

## Final composite figure
The journal-layout composite is in `composite/` and uses a 180 x 240 mm page:
- `Evolution_main_figure_FINAL_1086.svg`
- `Evolution_main_figure_FINAL_1086.pdf`
- `Evolution_main_figure_FINAL_1086_600dpi.png` (4252 x 5669 pixels)
- `Evolution_main_figure_FINAL_1086_600dpi.tif` (LZW RGB, 600 dpi)
- `Evolution_main_figure_FINAL_1086_preview.jpg`

The composite is raster-backed for reliable cross-platform assembly, while the individual A-G panel SVG/PDF sources remain in `main_panels/` for vector-based journal production and later layout refinement.