#!/usr/bin/env python3
"""Archive the public GEO metadata for HCATA's 17 muscle donors."""
import argparse
import concurrent.futures
import gzip
import hashlib
import json
from pathlib import Path
import urllib.request


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, default=Path("data/hcata_muscle_v1"))
    args = parser.parse_args()
    root = args.dataset
    accessions = ["GSE167186"] + [f"GSM{5098737 + i}" for i in range(17)]

    def fetch(accession):
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=self&form=text&view=full"
        name = "GSE167186_series" if accession.startswith("GSE") else accession
        path = root / "source" / (name + ".soft.gz")
        if not path.exists():
            with urllib.request.urlopen(url, timeout=60) as response:
                raw = response.read()
            if not raw.startswith(("^SERIES = " if accession.startswith("GSE") else "^SAMPLE = ").encode() + accession.encode()):
                raise ValueError(f"Unexpected GEO response: {accession}")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(gzip.compress(raw, mtime=0))
        data = path.read_bytes()
        return {"path": str(path.relative_to(root)), "url": url,
                "sha256": hashlib.sha256(data).hexdigest(),
                "uncompressed_sha256": hashlib.sha256(gzip.decompress(data)).hexdigest()}

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        manifest = list(executor.map(fetch, accessions))
    (root / "source" / "geo_downloads.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Verified {len(manifest)} GEO metadata snapshots")


if __name__ == "__main__":
    main()
