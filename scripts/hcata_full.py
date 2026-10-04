#!/usr/bin/env python3
"""Full-HCATA benchmark for Genomic Intelligence (GI) age-conditioned expression predictions.

Every gene HCATA reports, in every human study, scored with a GI expression model. Steps, each
resumable (cached downloads/responses are reused):

  genes      HCATA gene universe (topGenes/findAll)                   -> <root>/genes.tsv
  effects    HCATA differentialExpression for all genes, batched      -> <cache>/effects/*.json.gz
             + tables: every effect, and the primary (All Cells)      -> <root>/age_effects.tsv.gz,
               context per study scaled to that study's age contrast     <root>/primary_targets.tsv.gz
  sequences  Ensembl canonical TSS +/- flank, sense strand, GRCh38    -> <cache>/sequences.fa.gz,
                                                                         <root>/sequences.tsv.gz
  requests   one GI request per (gene, study, young/old)              -> <root>/gi_requests.tsv.gz
  predict    batched, rate-limited, concurrent GI calls               -> <cache>/gi/<model>/...,
                                                                         <root>/predictions_<model>.tsv.gz
  evaluate   Spearman of predicted vs reference age change per study  -> <root>/results_<model>/

Study tissue, the young/old ages and curation notes come from <root>/studies.tsv (checked against
PubMed; HCATA swaps two study labels). The reference for a gene in a study is HCATA's annual
log2 fold change times (age_old - age_young); the prediction is
log2(TPM_old + 1) - log2(TPM_young + 1). This tests the ranking of age changes, not TPM calibration.

GI auth: GI_API_KEY in the environment, or --key-file (default ~/.config/gi/api_key).
Only the Python standard library is needed, except `evaluate`, which uses numpy.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import gzip
import hashlib
import io
import json
import math
import os
import random
import re
import sys
import threading
import time
import urllib.error
import urllib.request
from decimal import Decimal
from pathlib import Path

HCATA = "http://134.84.61.115:3304/api/"
GI = os.environ.get("GI_BASE_URL", "https://api.genomicintelligence.ai")
ENSEMBL_FTP = "https://ftp.ensembl.org/pub/release-{r}/"
csv.field_size_limit(sys.maxsize)


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def read_tsv(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    fields = fields or (list(rows[0]) if rows else [])
    tmp = path.with_name(path.name + ".tmp")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(tmp, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def http_get(url: str, timeout: int = 600, tries: int = 4) -> bytes:
    last = None
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except Exception as exc:  # noqa: BLE001 - retried, then re-raised
            last = exc
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"GET {url} failed: {last}")


def download(url: str, path: Path) -> Path:
    if path.exists():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    log(f"downloading {url}")
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as fh:
        while chunk := r.read(1 << 22):
            fh.write(chunk)
    tmp.replace(path)
    return path


# ---------------------------------------------------------------- genes


def cmd_genes(a) -> None:
    if a.panel:  # fixed gene panel (e.g. hcata_muscle_v1/gene_panel.tsv) instead of all HCATA genes
        rows = [dict(ensembl_gene_id=r["ensembl_gene_id"], gene_symbol=r.get("gene_symbol", ""),
                     panel_group=r.get("hagr_group", r.get("panel_group", ""))) for r in read_tsv(a.panel)]
        write_tsv(a.root / "genes.tsv", rows)
        log(f"{len(rows)} panel genes from {a.panel} -> genes.tsv")
        for snap in a.import_effects:  # raw differentialExpression snapshots, e.g. source/effects_batch_*.json.gz
            recs = json.loads(gzip.decompress(Path(snap).read_bytes()))
            ids = sorted({f"ENSG{int(r['gene']):011d}" for r in recs})
            dest = a.cache / "effects" / f"batch_import_{hashlib.sha256(Path(snap).read_bytes()).hexdigest()[:12]}.json.gz"
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(gzip.compress(json.dumps({"genes": ids, "records": recs,
                                                       "imported_from": str(snap)}).encode(), mtime=0))
            log(f"imported {len(recs)} records for {len(ids)} genes from {snap}")
        return
    raw = a.cache / "hcata_findall.json.gz"
    if not raw.exists():
        raw.parent.mkdir(parents=True, exist_ok=True)
        raw.write_bytes(gzip.compress(http_get(HCATA + "topGenes/findAll"), mtime=0))
    genes = json.loads(gzip.decompress(raw.read_bytes()))
    rows = [dict(ensembl_gene_id=g["enID"], gene_symbol=g["gene"], chr=g["chr"], start=g["start"],
                 end=g["end"], hcata_total_sig=g["totalSig"]) for g in genes
            if re.fullmatch(r"ENSG\d{11}", str(g.get("enID", "")))]
    rows.sort(key=lambda r: r["ensembl_gene_id"])
    write_tsv(a.root / "genes.tsv", rows)
    log(f"{len(rows)} HCATA genes with Ensembl IDs ({len(genes) - len(rows)} without) -> genes.tsv")


# ---------------------------------------------------------------- effects


def cmd_effects(a) -> None:
    genes = [r["ensembl_gene_id"] for r in read_tsv(a.root / "genes.tsv")]
    if a.limit:
        genes = genes[: a.limit]
    out = a.cache / "effects"
    out.mkdir(parents=True, exist_ok=True)
    covered = set()
    for f in out.glob("batch_*.json.gz"):
        covered.update(json.loads(gzip.decompress(f.read_bytes()))["genes"])
    left = [g for g in genes if g not in covered]
    todo = list(enumerate([left[i : i + a.batch_size] for i in range(0, len(left), a.batch_size)]))
    log(f"effects: {len(covered)} genes cached, {len(left)} to fetch in {len(todo)} batches of "
        f"{a.batch_size}, {a.workers} workers")
    t0, done, failed = time.time(), 0, []

    def fetch(item):
        i, b = item
        ids = ",".join(str(int(g[4:])) for g in b)
        data = http_get(HCATA + "differentialExpression/" + ids, timeout=a.timeout)
        recs = json.loads(data)
        if not isinstance(recs, list):
            raise ValueError("unexpected response")
        path = out / f"batch_{b[0]}_{len(b)}.json.gz"  # named by content, not position
        tmp = path.with_name(path.name + ".tmp")
        tmp.write_bytes(gzip.compress(json.dumps({"genes": b, "records": recs}).encode(), mtime=0))
        tmp.replace(path)
        return len(recs)

    with cf.ThreadPoolExecutor(a.workers) as ex:
        futs = {ex.submit(fetch, t): t for t in todo}
        for f in cf.as_completed(futs):
            i, b = futs[f]
            try:
                f.result()
            except Exception as exc:  # noqa: BLE001
                # a whole batch can 500 on one bad gene: retry its genes one by one later
                failed.append(i)
                log(f"batch {i} failed: {exc}")
            done += 1
            if done % 20 == 0 or done == len(todo):
                rate = done / (time.time() - t0)
                log(f"effects {done}/{len(todo)} batches, {rate * 60:.1f}/min, "
                    f"ETA {(len(todo) - done) / max(rate, 1e-9) / 3600:.1f} h, {len(failed)} failed")
    if failed:
        log(f"{len(failed)} batches failed; rerun `effects` (e.g. --batch-size 5) to retry only "
            "their genes; cached genes are skipped")
    build_effect_tables(a)


def build_effect_tables(a) -> None:
    studies = {int(s["pmid"]): s for s in read_tsv(a.root / "studies.tsv")}
    wanted = {r["ensembl_gene_id"] for r in read_tsv(a.root / "genes.tsv")}
    effects, primary, seen = [], [], set()
    for path in sorted((a.cache / "effects").glob("batch_*.json.gz")):
        for r in json.loads(gzip.decompress(path.read_bytes()))["records"]:
            gene = f"ENSG{int(r['gene']):011d}"
            if gene not in wanted:
                continue
            labels = tuple(str(r.get(k, "")).strip() for k in ("cell_type", "cell_type2", "cell_type3"))
            key = (gene, r["pmid"], labels)
            if key in seen:
                continue
            seen.add(key)
            lfc, neglogp = Decimal(str(r["lfc"])), float(r["p_value"])
            row = dict(ensembl_gene_id=gene, pmid=r["pmid"], context=labels[0],
                       context_2=labels[1], context_3=labels[2],
                       lfc_log2_per_year_reported=str(r["lfc"]),
                       neg_log10_adjusted_p=str(r["p_value"]),
                       adjusted_p=repr(10 ** -neglogp) if math.isfinite(neglogp) else "")
            effects.append(row)
            s = studies.get(int(r["pmid"]))
            norm = lambda x: " ".join(x.split())  # noqa: E731 - labels differ in spacing
            if s and labels[0].startswith("All Cells") and norm(labels[0]) == norm(s["hcata_context_label"]):
                span = int(s["age_old"]) - int(s["age_young"])
                primary.append(dict(row, tissue=s["tissue"], age_young=s["age_young"],
                                    age_old=s["age_old"],
                                    reference_log2_change=str(lfc * span),
                                    direction=("up" if lfc > 0 else "down" if lfc < 0
                                               else "rounded_zero")))
    write_tsv(a.root / "age_effects.tsv.gz", effects)
    write_tsv(a.root / "primary_targets.tsv.gz", primary)
    by = {}
    for p in primary:
        by[p["pmid"]] = by.get(p["pmid"], 0) + 1
    log(f"{len(effects)} effects, {len(primary)} primary targets; per study: {by}")


# ---------------------------------------------------------------- sequences


def gi_models() -> dict:
    data = json.loads(http_get(GI + "/v1/tasks/expression/models", timeout=60))
    return {m["id"]: m for m in data["models"]}


def canonical_tss(gtf: Path, wanted: set[str]) -> dict[str, dict]:
    """Ensembl_canonical transcript TSS per gene (1-based genomic coordinate)."""
    out = {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "transcript" or 'tag "Ensembl_canonical"' not in f[8]:
                continue
            gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            if gid not in wanted:
                continue
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            tver = re.search(r'transcript_version "([^"]+)"', f[8])
            strand = f[6]
            out[gid] = dict(chrom=f[0], strand=strand, transcript=tid + (f".{tver.group(1)}" if tver else ""),
                            tss=int(f[3]) if strand == "+" else int(f[4]))
    return out


def iter_fasta(path: Path):
    name, buf = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.rstrip())
    if name:
        yield name, "".join(buf)


COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def cmd_sequences(a) -> None:
    flank = a.flank or gi_models()[a.model]["bio_spec"]["recommended_flank_bp"]
    base = ENSEMBL_FTP.format(r=a.ensembl_release)
    ens = a.cache / f"ensembl_{a.ensembl_release}"
    gtf = download(base + f"gtf/homo_sapiens/Homo_sapiens.GRCh38.{a.ensembl_release}.gtf.gz",
                   ens / "genes.gtf.gz")
    fa = download(base + "fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
                  ens / "genome.fa.gz")
    genes = {r["ensembl_gene_id"]: r for r in read_tsv(a.root / "genes.tsv")}
    tss = canonical_tss(gtf, set(genes))
    log(f"canonical TSS for {len(tss)}/{len(genes)} HCATA genes (Ensembl {a.ensembl_release}); "
        f"flank {flank} bp each side")
    by_chrom = {}
    for g, t in tss.items():
        by_chrom.setdefault(t["chrom"], []).append(g)
    rows, n_ok = [], 0
    fasta_out = a.cache / f"sequences_{a.root.name}_flank{flank}.fa.gz"
    tmp = fasta_out.with_name(fasta_out.name + ".tmp")
    with gzip.open(tmp, "wt") as out:
        for chrom, seq in iter_fasta(fa):
            for g in by_chrom.pop(chrom, []):
                t = tss[g]
                lo, hi = t["tss"] - 1 - flank, t["tss"] - 1 + flank + 1  # 0-based, TSS at index flank
                row = dict(ensembl_gene_id=g, gene_symbol=genes[g]["gene_symbol"], chrom=chrom,
                           strand=t["strand"], transcript=t["transcript"], tss_1based=t["tss"])
                if lo < 0 or hi > len(seq):
                    rows.append(dict(row, status="tss_too_close_to_contig_end"))
                    continue
                s = seq[lo:hi].upper()
                if t["strand"] == "-":
                    s = s.translate(COMP)[::-1]
                n_frac = s.count("N") / len(s)
                rows.append(dict(row, status="ok", length=len(s), tss_index=flank,
                                 n_fraction=f"{n_frac:.4f}",
                                 sequence_sha256=hashlib.sha256(s.encode()).hexdigest()))
                out.write(f">{g}\n{s}\n")
                n_ok += 1
    tmp.replace(fasta_out)
    for chrom, gs in by_chrom.items():  # contigs absent from the primary assembly
        rows += [dict(ensembl_gene_id=g, gene_symbol=genes[g]["gene_symbol"], chrom=chrom,
                      status="contig_not_in_primary_assembly") for g in gs]
    rows += [dict(ensembl_gene_id=g, gene_symbol=genes[g]["gene_symbol"],
                  status="no_canonical_transcript_in_release") for g in genes if g not in tss]
    write_tsv(a.root / "sequences.tsv.gz", rows,
              ["ensembl_gene_id", "gene_symbol", "status", "chrom", "strand", "transcript",
               "tss_1based", "length", "tss_index", "n_fraction", "sequence_sha256"])
    (a.root / "sequences_provenance.json").write_text(json.dumps(dict(
        assembly="GRCh38 primary assembly", ensembl_release=a.ensembl_release, gtf=base + "gtf/",
        tss="Ensembl_canonical transcript start, sense strand 5'->3'", flank_bp=flank,
        model_for_flank=a.model, fasta=str(fasta_out), n_ok=n_ok, n_genes=len(genes)), indent=1))
    log(f"{n_ok} sequences -> {fasta_out}")


# ---------------------------------------------------------------- requests


def description(study: dict, age: int) -> str:
    article = "an" if str(age).startswith("8") or age in (11, 18) else "a"
    return f"Human {study['description_tissue']} from {article} {age}-year-old donor."


def cmd_requests(a) -> None:
    studies = {s["pmid"]: s for s in read_tsv(a.root / "studies.tsv")}
    seqs = {r["ensembl_gene_id"]: r for r in read_tsv(a.root / "sequences.tsv.gz") if r["status"] == "ok"}
    rows, missing = [], 0
    for t in read_tsv(a.root / "primary_targets.tsv.gz"):
        if t["ensembl_gene_id"] not in seqs:
            missing += 1
            continue
        s = studies[t["pmid"]]
        for cond in ("young", "old"):
            age = int(s[f"age_{cond}"])
            rows.append(dict(request_id=f"{t['ensembl_gene_id']}__{t['pmid']}__{cond}",
                             ensembl_gene_id=t["ensembl_gene_id"], pmid=t["pmid"],
                             tissue=s["tissue"], condition=cond, age_years=age,
                             description=description(s, age)))
    write_tsv(a.root / "gi_requests.tsv.gz", rows)
    log(f"{len(rows)} GI requests ({len(rows) // 2} gene-study pairs); "
        f"{missing} primary targets without a sequence")


# ---------------------------------------------------------------- predict


class RateLimiter:
    """Token bucket that follows GI's RateLimit-* headers and backs off on 429."""

    def __init__(self, per_minute: float):
        self.rate = per_minute / 60
        self.tokens, self.t = 1.0, time.monotonic()
        self.lock = threading.Lock()
        self.pause_until = 0.0

    def acquire(self) -> None:
        while True:
            with self.lock:
                now = time.monotonic()
                if now >= self.pause_until:
                    self.tokens = min(max(1.0, self.rate * 10), self.tokens + (now - self.t) * self.rate)
                    self.t = now
                    if self.tokens >= 1:
                        self.tokens -= 1
                        return
                    wait = (1 - self.tokens) / self.rate
                else:
                    wait = self.pause_until - now
            time.sleep(min(wait, 5))

    def update(self, headers) -> None:
        limit = headers.get("RateLimit-Limit")
        if limit:
            try:  # burst = rpm / 6 (x-rate-limit-burst-divisor)
                rpm = float(str(limit).split(",")[0].split(";")[0]) * 6
                with self.lock:
                    self.rate = 0.9 * rpm / 60
            except ValueError:
                pass

    def back_off(self, seconds: float) -> None:
        with self.lock:
            self.pause_until = max(self.pause_until, time.monotonic() + seconds)


def load_fasta(path: Path, wanted: set[str]) -> dict[str, str]:
    return {n: s for n, s in iter_fasta(path) if n in wanted}


def cmd_predict(a) -> None:
    key = os.environ.get("GI_API_KEY") or (a.key_file.read_text().strip() if a.key_file.exists() else "")
    if not key:
        sys.exit(f"no GI key: set GI_API_KEY or write it to {a.key_file}")
    reqs = read_tsv(a.root / "gi_requests.tsv.gz")
    if a.pmid:
        reqs = [r for r in reqs if r["pmid"] in set(map(str, a.pmid))]
    if a.limit:
        genes = sorted({r["ensembl_gene_id"] for r in reqs})[: a.limit]
        reqs = [r for r in reqs if r["ensembl_gene_id"] in set(genes)]
    out = a.cache / "gi" / a.model
    path = lambda rid: out / rid[11:15] / f"{rid}.json"  # noqa: E731 - shard by gene id digits
    todo = [r for r in reqs if not path(r["request_id"]).exists()]
    prov = json.loads((a.root / "sequences_provenance.json").read_text())
    log(f"predict {a.model}: {len(reqs)} requests, {len(todo)} to run, {a.workers} workers, "
        f"batches of {a.batch_size}")
    seq_meta = {r["ensembl_gene_id"]: r for r in read_tsv(a.root / "sequences.tsv.gz") if r["status"] == "ok"}
    limiter = RateLimiter(a.rpm)
    errors_path = a.root / f"predict_errors_{a.model}.jsonl"
    t0, done, failed = time.time(), 0, 0

    def call(r, seq):
        body = json.dumps(dict(sequence=seq, sequence_name=r["ensembl_gene_id"], model=a.model,
                               tss_index=int(seq_meta[r["ensembl_gene_id"]]["tss_index"]),
                               options=dict(description=r["description"]))).encode()
        last = ""
        for attempt in range(a.retries):
            limiter.acquire()
            req = urllib.request.Request(GI + "/v1/tasks/expression/predict", data=body, method="POST",
                                         headers={"Authorization": f"Bearer {key}",
                                                  "Content-Type": "application/json"})
            try:
                with urllib.request.urlopen(req, timeout=a.timeout) as resp:
                    limiter.update(resp.headers)
                    data = json.loads(resp.read())
                p = path(r["request_id"])
                p.parent.mkdir(parents=True, exist_ok=True)
                tmp = p.with_name(p.name + ".tmp")
                tmp.write_text(json.dumps(dict(request=r, response=data,
                                               sequence_sha256=seq_meta[r["ensembl_gene_id"]]["sequence_sha256"])))
                tmp.replace(p)
                return None
            except urllib.error.HTTPError as e:
                limiter.update(e.headers)
                text = e.read().decode(errors="replace")[:500]
                last = f"{e.code} {text}"
                if e.code == 429 or e.code >= 500:
                    wait = float(e.headers.get("Retry-After") or 2 ** attempt)
                    limiter.back_off(wait + random.random())
                    continue
                return dict(request_id=r["request_id"], status=e.code, error=text)  # 4xx: final
            except Exception as exc:  # noqa: BLE001 - network: retry
                time.sleep(2 ** attempt)
                last = str(exc)
        return dict(request_id=r["request_id"], status="retries_exhausted", error=last)

    # batches of genes: load only their sequences, run all their requests concurrently
    genes = sorted({r["ensembl_gene_id"] for r in todo})
    gene_batches = [genes[i : i + a.batch_size] for i in range(0, len(genes), a.batch_size)]
    fasta = Path(prov["fasta"])
    with cf.ThreadPoolExecutor(a.workers) as ex:
        for bi, gb in enumerate(gene_batches, 1):
            seqs = load_fasta(fasta, set(gb))
            batch = [r for r in todo if r["ensembl_gene_id"] in seqs]
            for res in ex.map(lambda r: call(r, seqs[r["ensembl_gene_id"]]), batch):
                done += 1
                if res:
                    failed += 1
                    with open(errors_path, "a") as fh:
                        fh.write(json.dumps(res) + "\n")
            rate = done / (time.time() - t0)
            log(f"batch {bi}/{len(gene_batches)}: {done}/{len(todo)} requests, {failed} failed, "
                f"{rate * 60:.0f}/min, ETA {(len(todo) - done) / max(rate, 1e-9) / 3600:.1f} h")
    collect_predictions(a, reqs, path)


def collect_predictions(a, reqs, path) -> None:
    rows, models = [], set()
    for r in reqs:
        p = path(r["request_id"])
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        resp = d["response"]
        pred = resp["data"]["prediction"]
        models.add(resp["data"].get("model"))
        rows.append(dict(request_id=r["request_id"], ensembl_gene_id=r["ensembl_gene_id"],
                         pmid=r["pmid"], condition=r["condition"], age_years=r["age_years"],
                         tpm=pred.get("expression_tpm"), log_tpm=pred.get("expression_log_tpm"),
                         unit=pred.get("unit"), model=resp["data"].get("model"),
                         sequence_sha256=d["sequence_sha256"],
                         scored_window=json.dumps(resp["meta"].get("task_specific_counts", {}).get("scored_window")),
                         request_id_gi=resp["meta"].get("request_id")))
    write_tsv(a.root / f"predictions_{a.model}.tsv.gz", rows)
    log(f"{len(rows)}/{len(reqs)} predictions -> predictions_{a.model}.tsv.gz (models returned: {models})")


# ---------------------------------------------------------------- evaluate


def cmd_evaluate(a) -> None:
    import numpy as np

    def rank(x):
        order = np.argsort(x, kind="mergesort")
        r = np.empty(len(x))
        r[order] = np.arange(len(x))
        _, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)  # average ties
        sums = np.bincount(inv, weights=r)
        return (sums / cnt)[inv]

    def spearman(x, y):
        if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
            return None
        return float(np.corrcoef(rank(np.asarray(x)), rank(np.asarray(y)))[0, 1])

    preds = {}
    for p in read_tsv(a.root / f"predictions_{a.model}.tsv.gz"):
        preds[(p["ensembl_gene_id"], p["pmid"], p["condition"])] = float(p["tpm"])
    studies = {s["pmid"]: s for s in read_tsv(a.root / "studies.tsv")}
    pairs = {}
    for t in read_tsv(a.root / "primary_targets.tsv.gz"):
        k = (t["ensembl_gene_id"], t["pmid"])
        y, o = preds.get(k + ("young",)), preds.get(k + ("old",))
        if y is None or o is None:
            continue
        pairs.setdefault(t["pmid"], []).append(dict(
            ensembl_gene_id=k[0], pred=math.log2(o + 1) - math.log2(y + 1),
            ref=float(t["reference_log2_change"]), adj_p=float(t["adjusted_p"] or 1), tpm_young=y, tpm_old=o))
    rng = np.random.default_rng(a.seed)
    res = {"model": a.model, "studies": {}}
    out = a.root / f"results_{a.model}"
    out.mkdir(exist_ok=True)
    for pmid, rows in sorted(pairs.items()):
        x = np.array([r["pred"] for r in rows])
        y = np.array([r["ref"] for r in rows])
        nz = y != 0
        sig = np.array([r["adj_p"] < 0.05 for r in rows])
        boot = []
        for _ in range(a.bootstrap):
            i = rng.integers(0, len(x), len(x))
            v = spearman(x[i], y[i])
            if v is not None:
                boot.append(v)
        res["studies"][pmid] = dict(
            tissue=studies[pmid]["tissue"], age_young=int(studies[pmid]["age_young"]),
            age_old=int(studies[pmid]["age_old"]), n_genes=len(rows),
            spearman=spearman(x, y),
            spearman_95ci=(np.quantile(boot, [0.025, 0.975]).tolist() if boot else None),
            direction_accuracy_nonzero=(float(np.mean(np.sign(x[nz]) == np.sign(y[nz]))) if nz.any() else None),
            n_reference_significant=int(sig.sum()),
            spearman_significant_subset=spearman(x[sig], y[sig]),
            predicted_change_sd=float(np.std(x)), predicted_constant=bool(np.ptp(x) == 0))
        write_tsv(out / f"pairs_{pmid}.tsv.gz", rows)
    allx = [r["pred"] for rows in pairs.values() for r in rows]
    ally = [r["ref"] for rows in pairs.values() for r in rows]
    res["pooled_spearman_all_pairs"] = spearman(allx, ally)
    vals = [s["spearman"] for s in res["studies"].values() if s["spearman"] is not None]
    res["median_study_spearman"] = float(np.median(vals)) if vals else None
    res["note"] = ("Reference = HCATA annual log2FC x (age_old - age_young) for the study's All Cells "
                   "context; prediction = log2(TPM_old+1) - log2(TPM_young+1). HCATA values are "
                   "single-cell normalized counts, not TPM; bootstrap resamples genes, not donors.")
    (out / "summary.json").write_text(json.dumps(res, indent=1))
    for pmid, s in res["studies"].items():
        log(f"{pmid} {s['tissue']:<16} n={s['n_genes']:<6} rho={s['spearman']} "
            f"dir={s['direction_accuracy_nonzero']} sig={s['n_reference_significant']}")
    log(f"pooled rho {res['pooled_spearman_all_pairs']}, median study rho {res['median_study_spearman']}")


# ---------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=Path("data/hcata_full_v1"))
    ap.add_argument("--cache", type=Path, default=Path("/mnt/filesystem-c8/genes.jpg/agingexpress_cache"))
    ap.add_argument("--model", default="g0-expression-8192")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("genes")
    p.add_argument("--panel", type=Path, default=None, help="TSV with ensembl_gene_id: use this panel")
    p.add_argument("--import-effects", nargs="*", default=[],
                   help="existing raw HCATA differentialExpression snapshots (.json.gz) to reuse")
    p = sub.add_parser("effects")
    p.add_argument("--batch-size", type=int, default=25, help="genes per HCATA request (100 gave 500s)")
    p.add_argument("--workers", type=int, default=3, help="parallel HCATA requests (small server)")
    p.add_argument("--timeout", type=int, default=900)
    p.add_argument("--limit", type=int, default=0, help="first N genes only (testing)")
    sub.add_parser("effect-tables")
    p = sub.add_parser("sequences")
    p.add_argument("--ensembl-release", type=int, default=115)
    p.add_argument("--flank", type=int, default=0, help="bp each side of TSS; 0 = model's recommended")
    sub.add_parser("requests")
    p = sub.add_parser("predict")
    p.add_argument("--workers", type=int, default=16, help="concurrent GI requests")
    p.add_argument("--rpm", type=float, default=120, help="initial requests/min; adapts to RateLimit headers")
    p.add_argument("--batch-size", type=int, default=500, help="genes per batch (sequences loaded per batch)")
    p.add_argument("--retries", type=int, default=8)
    p.add_argument("--timeout", type=int, default=300)
    p.add_argument("--pmid", type=int, nargs="*", default=[], help="only these studies")
    p.add_argument("--limit", type=int, default=0, help="first N genes only (testing)")
    p.add_argument("--key-file", type=Path, default=Path.home() / ".config/gi/api_key")
    p = sub.add_parser("evaluate")
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    a.root.mkdir(parents=True, exist_ok=True)
    a.cache.mkdir(parents=True, exist_ok=True)
    {"genes": cmd_genes, "effects": cmd_effects, "effect-tables": build_effect_tables,
     "sequences": cmd_sequences, "requests": cmd_requests, "predict": cmd_predict,
     "evaluate": cmd_evaluate}[a.cmd](a)


if __name__ == "__main__":
    main()
