# AgingExpress

An exploratory benchmark for testing **Genomic Intelligence age-conditioned expression predictions** against measured human age-associated expression effects from [HCATA](http://hcata-xiaodonglab.org:3304/).

The first dataset, [`hcata_muscle_v1`](data/hcata_muscle_v1/), reuses a frozen 150-gene panel: 50 HAGR age-up, 50 age-down and 50 expression-matched background genes. It extracts HCATA's skeletal-muscle effects from PMID **36516485 / GSE167186**, with **17 donors aged 19–90**. There are **2,426 gene/context effects across 20 contexts**, with **124 primary gene targets** in the combined **All Cells, Skeletal Muscle** context. [`coverage.tsv`](data/hcata_muscle_v1/coverage.tsv) explicitly retains the 26 genes lacking that primary reference.

The dataset contains reference age-effect tables, donor/GEO metadata, exact GI context requests, an empty prediction template, raw source snapshots with SHA-256 hashes, extraction scripts and a scorer. **GI predictions have not been run for this dataset.** Counts and completed checks are recorded in [`validation_summary.json`](data/hcata_muscle_v1/validation_summary.json).

Start with the [GI test protocol](docs/gi_protocol.md) and [data dictionary](docs/data_dictionary.md). There are 300 age-20/age-75 requests for the full frozen panel; **248 requests have a primary reference**. Another 270 requests provide optional baseline/wording controls. The primary score compares the ranking of predicted age changes with HCATA's reported annual Log2FC scaled over 55 years. Source adjusted p-values are decoded from HCATA's −log10 display scale; effects are never selected by significance.

**Reference strength is limited:** none of the 124 primary effects passes HCATA's reported adjusted-p < 0.05 threshold, and 16 coefficients round to zero. A correlation with these point estimates is descriptive; it does not establish accuracy against confidently detected age effects. The fixed panel is retained rather than selecting stronger-looking genes after inspecting outcomes.

```bash
python -m pip install -r requirements.txt
python scripts/build_hcata.py
python -m unittest discover -s tests -v
# After filling the prediction template with actual GI outputs:
python scripts/evaluate_gi.py /path/to/predictions.tsv --output results/run_001
```

Python 3.12 was used for validation. Building from the committed snapshots requires only the standard library and no live HCATA connection. To retrieve current public responses, run `python scripts/download_hcata.py` and `python scripts/download_geo.py`, then rebuild. Downloads reuse existing cache files; use a separate copy of the dataset directory with its `source/` directory removed for a new snapshot. Do not overwrite a reference release silently.

The main limitations are the small donor cohort, a large age gap (six donors aged 19–22 and eleven aged 70–90), uncertain cell annotations, and unverified health/sex covariates. HCATA labels every donor healthy, but the original study includes sarcopenia/frailty; that label is retained as a source claim. Reference effects come from single-nucleus counts and are **not absolute TPM measurements**. The gene panel is selected and previous pilot model predictions had been inspected. Model-training overlap is unverified. These data support an exploratory age-response test, not a claim of independent or causal ageing validation.

Sources: Bartz et al., *Communications Biology* (2025), [HCATA paper](https://doi.org/10.1038/s42003-025-08845-8); *Aging* (2022), [muscle source study](https://doi.org/10.18632/aging.204435); [GEO GSE167186](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167186). Source terms are recorded in [SOURCE_TERMS.md](SOURCE_TERMS.md).

Report: [GI age-expression 150-gene report](reports/gi_age_expression_150_gene_report.docx).
