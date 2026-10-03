#!/usr/bin/env python3
"""Score GI predictions against the frozen HCATA muscle reference.

Only model outputs supplied by the user are scored. This script does not invoke
GI and does not treat blank template rows as predictions.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np
from scipy.stats import spearmanr


def read_tsv(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def spearman(x, y):
    if len(x) < 3 or len(set(x)) < 2 or len(set(y)) < 2:
        return None
    value = float(spearmanr(x, y).statistic)
    return value if math.isfinite(value) else None


def summarize(pairs, seed, iterations):
    x = np.array([r["predicted_log2_change_plus1"] for r in pairs])
    y = np.array([r["reference_scaled_log2_change"] for r in pairs])
    rho = spearman(x.tolist(), y.tolist())
    result = {"n_genes": len(pairs), "spearman_rho": rho,
              "bootstrap_unit": "gene, conditional on this frozen panel; not donor uncertainty"}
    rng = np.random.default_rng(seed)
    if rho is not None and iterations:
        boot = []
        for _ in range(iterations):
            idx = rng.integers(0, len(x), len(x))
            value = spearman(x[idx].tolist(), y[idx].tolist())
            if value is not None:
                boot.append(value)
        null = [spearman(x.tolist(), rng.permutation(y).tolist()) for _ in range(iterations)]
        null = [v for v in null if v is not None]
        result["gene_bootstrap_95pct_interval"] = np.quantile(boot, [0.025, 0.975]).tolist() if boot else None
        result["two_sided_gene_permutation_p"] = (1 + sum(abs(v) >= abs(rho) for v in null)) / (1 + len(null))
    nonzero = [r for r in pairs if r["reference_scaled_log2_change"] != 0]
    result["direction_n_nonzero_reference"] = len(nonzero)
    result["direction_accuracy_nonzero_reference"] = sum(np.sign(r["predicted_log2_change_plus1"]) == np.sign(r["reference_scaled_log2_change"]) for r in nonzero) / len(nonzero) if nonzero else None
    result["direction_up_reference_n"] = sum(r["reference_scaled_log2_change"] > 0 for r in nonzero)
    result["direction_down_reference_n"] = sum(r["reference_scaled_log2_change"] < 0 for r in nonzero)
    positive = [r for r in pairs if r["tpm_young"] > 0 and r["tpm_old"] > 0]
    result["no_pseudocount_sensitivity_n"] = len(positive)
    result["no_pseudocount_sensitivity_spearman"] = spearman([math.log2(r["tpm_old"]) - math.log2(r["tpm_young"]) for r in positive], [r["reference_scaled_log2_change"] for r in positive])
    significant = [r for r in pairs if r["adjusted_p_below_0_05"]]
    result["significant_reference_subset_n"] = len(significant)
    result["significant_reference_subset_spearman_exploratory"] = spearman([r["predicted_log2_change_plus1"] for r in significant], [r["reference_scaled_log2_change"] for r in significant])
    return result


def score(dataset, predictions, allow_incomplete=False, iterations=2000, seed=42):
    requests = {r["request_id"]: r for r in read_tsv(dataset / "gi_requests.tsv")}
    targets = {r["ensembl_gene_id"]: r for r in read_tsv(dataset / "primary_targets.tsv")}
    completed, seen, hashes, run_ids = {}, set(), {}, set()
    for row in predictions:
        request_id = row.get("request_id", "")
        if request_id not in requests or request_id in seen:
            raise ValueError(f"Unknown or duplicate request_id: {request_id}")
        seen.add(request_id)
        if not row.get("tpm", "").strip():
            continue
        tpm = float(row["tpm"])
        if not math.isfinite(tpm) or tpm < 0:
            raise ValueError(f"Invalid TPM: {request_id}")
        request = requests[request_id]
        gene = request["ensembl_gene_id"]
        if row.get("resolved_ensembl_gene_id") != gene:
            raise ValueError(f"Ensembl resolution mismatch: {request_id}")
        if row.get("model") != request["model"]:
            raise ValueError(f"Model mismatch: {request_id}")
        sequence_hash = row.get("sequence_sha256", "")
        if not re.fullmatch(r"[0-9a-f]{64}", sequence_hash):
            raise ValueError(f"Missing/invalid sequence hash: {request_id}")
        if gene in hashes and hashes[gene] != sequence_hash:
            raise ValueError(f"Different sequences used for one gene: {gene}")
        hashes[gene] = sequence_hash
        if not row.get("run_id", "").strip():
            raise ValueError(f"Missing run_id: {request_id}")
        run_ids.add(row["run_id"])
        completed[request_id] = tpm
    if len(run_ids) > 1:
        raise ValueError("Mixed run_id values; score each experiment separately")
    expected = [r for r in requests.values() if r["variant"] == "primary" and r["ensembl_gene_id"] in targets]
    missing = [r["request_id"] for r in expected if r["request_id"] not in completed]
    if missing and not allow_incomplete:
        raise ValueError(f"Missing {len(missing)} primary predictions. Complete the run or use --allow-incomplete and report coverage.")
    pairs = []
    for gene, target in targets.items():
        for variant in ("primary", "alt1", "alt2"):
            keys = [f"{gene}__{variant}_{c}" for c in ("young", "old")]
            if not all(k in completed for k in keys):
                continue
            young, old = [completed[k] for k in keys]
            pairs.append({"ensembl_gene_id": gene, "gene_symbol": target["gene_symbol"],
                          "variant": variant, "tpm_young": young, "tpm_old": old,
                          "predicted_log2_change_plus1": math.log2(old + 1) - math.log2(young + 1),
                          "reference_scaled_log2_change": float(target["scaled_log2_change_age20_to75"]),
                          "adjusted_p_below_0_05": int(target["adjusted_p_below_0_05"])})
    primary = [r for r in pairs if r["variant"] == "primary"]
    metrics = {"dataset_id": "hcata_muscle_v1", "run_id": next(iter(run_ids), None),
               "seed": seed, "iterations": iterations, "completed_request_count": len(completed),
               "reference_gene_count": len(targets), "primary_pair_count": len(primary),
               "reference_significant_gene_count": sum(int(t['adjusted_p_below_0_05']) for t in targets.values()),
               "reference_rounded_zero_gene_count": sum(float(t['scaled_log2_change_age20_to75']) == 0 for t in targets.values()),
               "missing_primary_request_ids": missing,
               "primary_all_available_genes": summarize(primary, seed, iterations),
               "warning": "Exploratory study-level benchmark with no source-significant primary effects in v1. HCATA normalized nucleus counts are not TPM; age effects and gene bootstrap intervals do not establish validated ageing accuracy, causal ageing or donor-level precision."}
    # Compare alternative wording on the SAME genes, rather than a 150-vs-30 comparison.
    for variant in ("alt1", "alt2"):
        alternative = {r["ensembl_gene_id"]: r for r in pairs if r["variant"] == variant}
        comparable = [r for r in primary if r["ensembl_gene_id"] in alternative]
        alt_pairs = [alternative[r["ensembl_gene_id"]] for r in comparable]
        metrics[f"wording_{variant}"] = {"primary_on_same_genes": summarize(comparable, seed, iterations),
                                        "alternative_on_same_genes": summarize(alt_pairs, seed, iterations),
                                        "mean_absolute_change_in_predicted_delta": float(np.mean([abs(a["predicted_log2_change_plus1"] - b["predicted_log2_change_plus1"]) for a,b in zip(comparable, alt_pairs)])) if comparable else None}
    return metrics, pairs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("predictions", type=Path)
    parser.add_argument("--dataset", type=Path, default=Path("data/hcata_muscle_v1"))
    parser.add_argument("--output", type=Path, default=Path("results/gi_hcata_muscle"))
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--iterations", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if args.iterations < 0:
        parser.error("--iterations must be nonnegative")
    metrics, pairs = score(args.dataset, read_tsv(args.predictions), args.allow_incomplete, args.iterations, args.seed)
    metrics["prediction_file_sha256"] = hashlib.sha256(args.predictions.read_bytes()).hexdigest()
    metrics["reference_file_sha256"] = hashlib.sha256((args.dataset / "primary_targets.tsv").read_bytes()).hexdigest()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=False) + "\n")
    if pairs:
        with (args.output / "paired_predictions.tsv").open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(pairs[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(pairs)
    print(json.dumps(metrics["primary_all_available_genes"], indent=2))


if __name__ == "__main__":
    main()
