# GI expression test protocol

This is an exploratory follow-up to an earlier 150-gene pilot, using a new HCATA reference. Model predictions from that pilot had already been inspected. The panel is fixed before extraction of HCATA outcomes, but this is not a preregistered or blinded validation.

## Frozen design

- Model: `g0-expression`, as exposed by Genomic Intelligence on 2026-10-03. Record the actual model/version returned during execution. Changing the model defines a separate experiment.
- Organism/tissue: human skeletal muscle. No sex or health status is asserted in the prediction descriptions.
- Ages: 20 and 75, both represented in the HCATA cohort. There are no donors aged 23–69; the fitted annual effect does not establish a continuous trajectory through that gap.
- Primary target: HCATA's **All Cells, Skeletal Muscle** effect in PMID 36516485. This is a combined single-nucleus context, not a bulk RNA-seq TPM measurement.
- Panel: 50 HAGR up, 50 HAGR down and 50 expression-matched background genes. HAGR labels are selection metadata, not the evaluation target. Preserve all available reference genes, including nonsignificant and rounded-zero effects. Report coverage instead of replacing missing genes.
- Primary comparison: Spearman correlation across genes between `log2(TPM75 + 1) - log2(TPM20 + 1)` and `55 × HCATA reported annual Log2FC`. This tests ranking of age changes; it does not test absolute TPM calibration. The scaled reference is a fitted contrast, not a measured paired 20-to-75 change.
- Secondary comparisons: direction accuracy on nonzero reported effects, ratio sensitivity without a pseudocount for positive predictions, and the source-adjusted-p < 0.05 subset (explicitly exploratory).

## Run GI manually or with the connector

1. Read `data/hcata_muscle_v1/gi_requests.tsv`. Start with `variant=primary` and `primary_target_available=1`: 248 requests for 124 genes. The full frozen panel has 300 primary requests; 26 genes lack a primary reference and are not scored. The remaining rows provide 150 age-unspecified baselines and two alternative wordings on 30 frozen control genes (120 requests).
2. For each unique gene, call `fetch_gene_for_expression(gene=<ensembl_gene_id>, species="human")`. Confirm the resolved Ensembl ID matches the panel. Use the exact TSS-centred **9,198 bp** window prepared by this tool. Do not use a whole gene body as the expression input.
3. Save sequence provenance: reference assembly/release, transcript/TSS choice, genomic coordinates and strand, length, and SHA-256 of the uppercase DNA string. Reuse exactly the same sequence for every condition and wording of that gene. GI sequence handles are temporary; they are not durable sequence identities.
4. Call `predict_expression(sequence_ref=<handle>, description=<exact description column>, model="g0-expression")`. Keep the entire returned response and its metadata. Do not include HAGR labels, HCATA effects, or gene-age expectations in the context text. Read the documented TPM output; do not assume another output unit is TPM.
5. Copy `predictions_template.tsv` to a run-specific file. Populate the TPM, sequence hash, resolved Ensembl ID, model and a single experiment `run_id`. Leave unrun control rows blank. Log sequence-fetch failures and model errors separately; never enter zero to represent a failed request.
6. Score the predictions:

```bash
python -m pip install -r requirements.txt
python scripts/evaluate_gi.py /path/to/predictions.tsv --output results/run_001
```

The evaluator rejects duplicate/unknown request IDs, negative/nonfinite TPM, Ensembl/model mismatches, mixed run IDs, changed sequences between conditions, and missing primary predictions. For documented execution failures, use `--allow-incomplete`; the output lists missing requests and the number of complete gene pairs. Controls are optional. Alternative-wording comparisons use the same genes in both wordings.

## Interpretation

Report the primary result first, with gene coverage and rounded-zero reference count. **None of the 124 primary targets passes the source adjusted-p < 0.05 threshold; 16 coefficients round to zero.** These are uncertain point estimates, not confidently established gene directions. A correlation is descriptive and should not be called validated ageing accuracy. A constant prediction or reference produces an undefined correlation, recorded as `null`, rather than a misleading zero. The 2,000-iteration bootstrap resamples genes; its interval describes this panel's cross-gene statistic and does not estimate biological uncertainty across donors. The permutation check also concerns genes in this selected panel. Neither controls for source-cohort health, sex or cell-composition confounding.

Report the 30-gene wording sensitivity separately, including the mean absolute change in predicted age contrast. Baseline requests are useful for examining how adding age changes predictions; the reference has no matched TPM measurements for validating baseline levels.

HCATA reports annual coefficients to three decimal places through the API. A displayed zero can hide a small effect; its direction is unresolved. The source-adjusted p is retained independently of this rounded coefficient. Treat ties accordingly. Do not optimize descriptions, pseudocounts, gene filters, or model settings against this dataset and then describe the resulting score as held-out performance.

GI's current model catalog lists ENCODE bulk, GTEx/recount3 and cellxgene pseudobulk among training sources. Overlap with this study has not been checked. This benchmark cannot establish independence from training data. Record model metadata and resolve overlap before claiming external validation.
