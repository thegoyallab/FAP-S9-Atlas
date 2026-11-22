# Targetability analysis boundaries

The targetability extension was designed to separate four properties that should not be conflated: cancer-versus-normal target enrichment, fibroblast biological state, cellular source/selectivity and nonmalignant remodeling background.

- Cancer-versus-normal analyses use organ-concordant disease-site rules; discordant metastases and hematologic malignancies are excluded from the primary solid-tumor inference.
- The strongest tumor-normal tier uses within-donor exact-tissue pairs with at least five pairs.
- Cross-cancer state analyses compare FAP-only with DPP4-only fibroblasts within the same donor; anatomical variants of the same organ are donor-collapsed before meta-analysis.
- The benign-remodeling comparison is anchored to one lung collection containing healthy, fibrotic, inflammatory/infectious and cancer donors in the same exact tissue.
- Cellular-selectivity analyses compare the curated primary-fibroblast compartment with all other annotated cells without imposing a new broad-lineage ontology; exact Cell Ontology labels are retained for off-compartment signals.
- No composite targetability score or cancer ranking is calculated.
- RNA abundance/detection does not establish surface protein, enzymatic activity, PET SUV, radioligand binding or dosimetry.
