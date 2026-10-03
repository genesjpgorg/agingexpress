# Data dictionary

All curated data are UTF-8, tab-separated files with a header. Empty fields mean missing/unrun; they do not mean zero. The immutable source snapshot and checksums are in `data/hcata_muscle_v1/source/`.

| File | Row unit | Purpose |
|---|---|---|
| `gene_panel.tsv` | Frozen gene | Ensembl/symbol/Entrez identifiers, HAGR selection group, and 30-gene wording-control flag |
| `age_effects.tsv` | Gene × HCATA muscle context | All returned muscle effects for the panel; no significance filter |
| `primary_targets.tsv` | Gene | The combined All Cells, Skeletal Muscle context only |
| `contexts.tsv` | Context | Primary and two alternative HCATA annotations; numbered clusters and Unknown labels are retained |
| `donors.tsv` | Donor | HCATA sample ID, age/sex/health labels, original GEO accession/group and verification status |
| `coverage.tsv` | Frozen gene | Returned context count, primary availability, and explicit missing status |
| `gi_requests.tsv` | Model request | Exact gene, age, model and context text; 300 primary plus 270 optional control requests |
| `predictions_template.tsv` | Model request | Empty fields for actual GI outputs and execution provenance |
| `panel_provenance.json` | Panel snapshot | Earlier selection method, seed, control IDs, source hashes and original panel hash |
| `validation_summary.json` | Build | Counts and integrity checks; no model evaluation results |

## Effect columns

| Column | Meaning |
|---|---|
| `ensembl_gene_id`, `gene_symbol` | Frozen mapping; HCATA numeric gene IDs are left-padded into stable ENSG identifiers |
| `study_pmid`, `geo_series` | Source cohort: PMID 36516485 / GSE167186 |
| `context_id` | Deterministic hash of the three HCATA context labels; not a biological cell ID |
| `cell_type_hcata`, `cell_type_alternative_2/3` | Source annotation strings, outer whitespace trimmed; alternative labels expose annotation uncertainty |
| `primary_tissue_context` | 1 for the combined All Cells, Skeletal Muscle context; otherwise 0 |
| `n_donors_in_study`, `age_min_years`, `age_max_years` | Study-level coverage, not verified donor counts for each cluster |
| `lfc_log2_per_year_reported` | HCATA `lfc`, described by the paper as annual Log2FC; API precision retained |
| `slope_reported`, `intercept_reported` | Source fit/display values, retained for auditing; intercept is not treated as TPM |
| `neg_log10_adjusted_p` | Raw HCATA field named `p_value`; web-client scale and threshold identify it as −log10(adjusted p) |
| `adjusted_p_derived` | `10 ** (-neg_log10_adjusted_p)`; not a new multiple-testing correction |
| `adjusted_p_below_0_05` | Source adjusted p < 0.05 flag; no additional correction within this selected panel |
| `scaled_log2_change_age20_to75` | `55 × lfc`; model-scaled contrast, not raw measured TPM or paired change |
| `direction_reported` | up/down/rounded_zero based on the reported coefficient; zero is unresolved at source precision |
| `hcata_plotting_id` | Public source plot identifier, useful for tracing the record in HCATA |
| `source_path` | Archived API response containing this record |

HCATA's `sex=unknown` and `disease_status=healthy` are copied as source claims, not repaired or verified. GEO confirms donor ages and Young/Old grouping but does not establish donor sex or rule out sarcopenia/frailty. No within-donor expression matrix is included. API scatterplot values lack donor identifiers and are not used as donor-level measurements. Raw effect snapshots include other tissues returned for the same genes, but the curated benchmark uses only the muscle cohort.
