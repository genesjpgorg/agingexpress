#!/usr/bin/env python3
"""Build a deterministic GI benchmark from the archived public HCATA responses."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import json
import math
from pathlib import Path
import re

PMID = 36516485
AGES = {"young": 20, "old": 75}


def read_tsv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    fields = fields or list(rows[0])
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def raw_json(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def verify_sources(root):
    entries = []
    for filename in ("downloads.json", "geo_downloads.json"):
        entries.extend(json.loads((root / "source" / filename).read_text()))
    for row in entries:
        data = (root / row["path"]).read_bytes()
        if hashlib.sha256(data).hexdigest() != row["sha256"]:
            raise ValueError(f"Checksum mismatch: {row['path']}")
        if hashlib.sha256(gzip.decompress(data)).hexdigest() != row["uncompressed_sha256"]:
            raise ValueError(f"Uncompressed checksum mismatch: {row['path']}")
    return len(entries)


def geo_metadata(source):
    records = {}
    for path in sorted(source.glob("GSM*.soft.gz")):
        text = gzip.decompress(path.read_bytes()).decode()
        donor = re.search(r"^!Sample_description = (HM\d+)\s*$", text, re.M).group(1)
        age = re.search(r"^!Sample_characteristics_ch1 = age: (\d+)", text, re.M).group(1)
        group = re.search(r"^!Sample_characteristics_ch1 = group: (.+)$", text, re.M).group(1).strip()
        records[donor] = {"geo_accession": path.name.split(".")[0], "geo_age_years": age,
                          "geo_group": group}
    return records


def build(root):
    source_count = verify_sources(root)
    panel = read_tsv(root / "gene_panel.tsv")
    if len(panel) != 150 or len({r["ensembl_gene_id"] for r in panel}) != 150:
        raise ValueError("Expected 150 distinct frozen genes")
    for row in panel:
        if not re.fullmatch(r"ENSG\d{11}", row["ensembl_gene_id"]):
            raise ValueError(f"Invalid Ensembl identifier: {row}")
    by_number = {int(r["ensembl_gene_id"][4:]): r for r in panel}
    all_effects = []
    for path in sorted((root / "source").glob("effects_batch_*.json.gz")):
        for row in raw_json(path):
            if int(row["gene"]) not in by_number:
                raise ValueError("HCATA returned a gene outside the frozen panel")
            all_effects.append((row, str(path.relative_to(root))))

    effects, seen = [], set()
    for row, source_path in all_effects:
        if int(row["pmid"]) != PMID:
            continue
        gene = by_number[int(row["gene"])]
        labels = tuple(str(row[k]).strip() for k in ("cell_type", "cell_type2", "cell_type3"))
        context_id = "ctx_" + hashlib.sha256("|".join(labels).encode()).hexdigest()[:12]
        key = (gene["ensembl_gene_id"], context_id)
        if key in seen:
            raise ValueError(f"Duplicate effect: {key}")
        seen.add(key)
        lfc = Decimal(str(row["lfc"]))
        neglogp = float(row["p_value"])
        if not math.isfinite(float(lfc)) or not math.isfinite(neglogp) or neglogp < 0:
            raise ValueError(f"Invalid source effect: {row}")
        p_adjusted = 10 ** (-neglogp)
        is_primary = " ".join(labels[0].lower().split()) == "all cells, skeletal muscle"
        effects.append({
            "ensembl_gene_id": gene["ensembl_gene_id"], "gene_symbol": gene["gene_symbol"],
            "study_pmid": PMID, "geo_series": "GSE167186", "context_id": context_id,
            "cell_type_hcata": labels[0], "cell_type_alternative_2": labels[1],
            "cell_type_alternative_3": labels[2], "primary_tissue_context": int(is_primary),
            "n_donors_in_study": 17, "age_min_years": 19, "age_max_years": 90,
            "lfc_log2_per_year_reported": str(row["lfc"]),
            "slope_reported": str(row["slope"]), "intercept_reported": str(row["inter"]),
            "neg_log10_adjusted_p": str(row["p_value"]),
            "adjusted_p_derived": format(p_adjusted, ".17g"),
            "adjusted_p_below_0_05": int(neglogp > -math.log10(0.05)),
            "scaled_log2_change_age20_to75": str(lfc * Decimal(55)),
            "direction_reported": "up" if lfc > 0 else "down" if lfc < 0 else "rounded_zero",
            "hcata_plotting_id": row["plotting_id"], "source_path": source_path,
        })
    if not effects:
        raise ValueError("No muscle effects retrieved")
    effects.sort(key=lambda r: (int(next(p["panel_order"] for p in panel if p["ensembl_gene_id"] == r["ensembl_gene_id"])), r["context_id"]))
    write_tsv(root / "age_effects.tsv", effects)
    primary = [r for r in effects if r["primary_tissue_context"]]
    if len(primary) != len({r["ensembl_gene_id"] for r in primary}):
        raise ValueError("More than one primary effect per gene")
    write_tsv(root / "primary_targets.tsv", primary, list(effects[0]))
    by_gene = {r["ensembl_gene_id"]: r for r in primary}
    coverage = []
    for row in panel:
        gene_id = row["ensembl_gene_id"]
        coverage.append({"ensembl_gene_id": gene_id, "gene_symbol": row["gene_symbol"],
                         "n_muscle_contexts": sum(r["ensembl_gene_id"] == gene_id for r in effects),
                         "primary_available": int(gene_id in by_gene),
                         "status": "available" if gene_id in by_gene else "missing_hcata_all_cells_effect"})
    write_tsv(root / "coverage.tsv", coverage)

    contexts = {r["context_id"]: {k: r[k] for k in ("context_id", "cell_type_hcata", "cell_type_alternative_2", "cell_type_alternative_3", "primary_tissue_context")} for r in effects}
    write_tsv(root / "contexts.tsv", [contexts[k] for k in sorted(contexts)])

    geo = geo_metadata(root / "source")
    samples = [s for s in raw_json(root / "source" / "hcata_samples.json.gz") if int(s["study_id"]) == PMID]
    donors = []
    for sample in sorted(samples, key=lambda s: int(s["sample_id"])):
        donor = re.search(r"HM\d+", sample["notes"]).group(0)
        original = geo[donor]
        if float(original["geo_age_years"]) != float(sample["age"]):
            raise ValueError(f"HCATA/GEO age mismatch: {donor}")
        donors.append({"study_pmid": PMID, "hcata_sample_id": sample["sample_id"],
                       "donor_id": donor, "age_years": sample["age"],
                       "sex_hcata": sample["sex"], "disease_status_hcata": sample["disease_status"],
                       "tissue_hcata": sample["tissue"],
                       **original, "health_status_verified": "no",
                       "covariate_warning": "Portal healthy label is not independently confirmed; source study includes sarcopenia/frailty."})
    if len(donors) != 17 or len(geo) != 17 or len({r['donor_id'] for r in donors}) != 17:
        raise ValueError("Expected 17 donor metadata records")
    write_tsv(root / "donors.tsv", donors)

    requests = []
    for gene in panel:
        variants = ["primary", "baseline"] + (["alt1", "alt2"] if gene["wording_audit"] == "1" else [])
        for variant in variants:
            conditions = ["baseline"] if variant == "baseline" else ["young", "old"]
            for condition in conditions:
                age = AGES.get(condition, "")
                if variant == "baseline":
                    description = "Human skeletal muscle tissue from an adult donor."
                elif variant == "primary":
                    description = f"Human skeletal muscle tissue from a {age}-year-old adult donor."
                elif variant == "alt1":
                    description = f"Human skeletal muscle tissue from an adult donor aged {age} years."
                else:
                    description = f"Skeletal muscle tissue from a human adult donor who is {age} years old."
                requests.append({"request_id": f"{gene['ensembl_gene_id']}__{variant}_{condition}",
                                 "ensembl_gene_id": gene["ensembl_gene_id"], "gene_symbol": gene["gene_symbol"],
                                 "model": "g0-expression", "variant": variant, "condition": condition,
                                 "age_years": age, "description": description,
                                 "primary_target_available": int(gene["ensembl_gene_id"] in by_gene)})
    write_tsv(root / "gi_requests.tsv", requests)
    write_tsv(root / "predictions_template.tsv", [{"request_id": r["request_id"], "tpm": "", "sequence_sha256": "", "resolved_ensembl_gene_id": "", "model": "", "run_id": ""} for r in requests])
    counts = {"dataset_id": "hcata_muscle_v1", "snapshot_date_utc": "2026-10-03",
              "study_pmid": PMID, "geo_series": "GSE167186", "frozen_genes": len(panel),
              "primary_target_genes": len(primary), "muscle_effect_rows": len(effects),
              "muscle_contexts": len(contexts), "source_effect_rows_all_tissues": len(all_effects),
              "donors": len(donors), "age_min": min(int(r["age_years"]) for r in donors),
              "age_max": max(int(r["age_years"]) for r in donors),
              "donors_age_19_to_22": sum(19 <= int(r["age_years"]) <= 22 for r in donors),
              "donors_age_70_to_90": sum(70 <= int(r["age_years"]) <= 90 for r in donors),
              "primary_adjusted_p_below_0_05": sum(r["adjusted_p_below_0_05"] for r in primary),
              "primary_reported_rounded_zero": sum(r["direction_reported"] == "rounded_zero" for r in primary),
              "primary_gi_requests": sum(r["variant"] == "primary" for r in requests),
              "total_gi_requests_including_controls": len(requests),
              "verified_archived_source_files": source_count,
              "health_labels_independently_verified": False,
              "gi_predictions_run": False,
              "checks_passed": ["source checksums", "unique frozen identifiers", "unique gene/context effects", "finite effect/probability fields", "one primary effect per gene", "HCATA and GEO donor ages agree", "17 distinct GEO donor records"]}
    (root / "validation_summary.json").write_text(json.dumps(counts, indent=2) + "\n")
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, default=Path("data/hcata_muscle_v1"))
    build(parser.parse_args().dataset)
