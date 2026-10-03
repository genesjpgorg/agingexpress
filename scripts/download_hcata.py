#!/usr/bin/env python3
"""Download public HCATA responses for a frozen Ensembl gene panel.

The public Angular client converts ENSG identifiers to integers and requests
comma-separated batches from differentialExpression. No outcome filtering is
applied here. Each exact response is retained, compressed deterministically.
"""
import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import json
from pathlib import Path
import time
import urllib.request

BASE = "http://134.84.61.115:3304/api/"  # Public API URL in HCATA's web client.


def fetch(url, path, timeout):
    if path.exists():
        return json.loads(gzip.decompress(path.read_bytes()))
    last = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                data = response.read()
            parsed = json.loads(data)
            path.parent.mkdir(parents=True, exist_ok=True)
            temp = path.with_suffix(path.suffix + ".tmp")
            temp.write_bytes(gzip.compress(data, mtime=0))
            temp.replace(path)
            return parsed
        except Exception as exc:
            last = exc
            if attempt < 2:
                time.sleep(attempt + 1)
    raise RuntimeError(f"Could not retrieve {url}: {last}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, default=Path("data/hcata_muscle_v1"))
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=300)
    args = parser.parse_args()
    root = args.dataset
    with (root / "gene_panel.tsv").open() as stream:
        panel = list(csv.DictReader(stream, delimiter="\t"))
    ids = [str(int(row["ensembl_gene_id"][4:])) for row in panel]
    source = root / "source"
    metadata_base = "http://hcata-xiaodonglab.org:3304/api/"
    tasks = [(metadata_base + "sample", source / "hcata_samples.json.gz"),
             (metadata_base + "metadata", source / "hcata_studies.json.gz")]
    for start in range(0, len(ids), 25):
        tasks.append((BASE + "differentialExpression/" + ",".join(ids[start:start + 25]),
                      source / f"effects_batch_{start // 25:03d}.json.gz"))

    def run(task):
        url, path = task
        result = fetch(url, path, args.timeout)
        if not isinstance(result, list):
            raise ValueError(f"Unexpected response from {url}")
        print(f"{path.name}: {len(result)} records", flush=True)
        return {"path": str(path.relative_to(root)), "url": url,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "uncompressed_sha256": hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest(),
                "records": len(result)}

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        manifest = list(executor.map(run, tasks))
    (source / "downloads.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
