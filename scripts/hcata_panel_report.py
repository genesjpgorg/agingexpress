#!/usr/bin/env python3
"""Statistics, figures and DOCX report for the 150-gene GI study across 12 HCATA studies.

Reads data/hcata_panel150_v1 (built by scripts/hcata_full.py: panel, HCATA primary targets,
sequences, GI predictions for g0-expression and g0-expression-8192) and writes
reports/gi_hcata_12_study_150_gene_report.docx plus its figures and a stats JSON.
The existing 150-gene muscle report is used as the style template.

Requires numpy, matplotlib and python-docx.
"""

from __future__ import annotations

import argparse
import copy
import csv
import gzip
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

MODELS = ["g0-expression", "g0-expression-8192"]
LABEL = {"g0-expression": "g0-expression (9.2 kb)", "g0-expression-8192": "g0-expression-8192 (82 kb)"}
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
COLOR = {"g0-expression": "#2a78d6", "g0-expression-8192": "#eb6834"}
DIVERGING = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#e34948"])


def read_tsv(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def rank(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x))
    r[order] = np.arange(len(x))
    _, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
    return (np.bincount(inv, weights=r) / cnt)[inv]


def spearman(x, y) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan")
    return float(np.corrcoef(rank(x), rank(y))[0, 1])


def style(ax, title, xlabel, ylabel):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, loc="left", color=INK, fontsize=9.5)
    ax.set_xlabel(xlabel, color=INK2, fontsize=8.5)
    ax.set_ylabel(ylabel, color=INK2, fontsize=8.5)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=7.5)


# ---------------------------------------------------------------- analysis


def analyse(root: Path, raw_dir: Path, perms: int, seed: int) -> dict:
    studies = {s["pmid"]: s for s in read_tsv(root / "studies.tsv")}
    targets = read_tsv(root / "primary_targets.tsv.gz")
    ref = {(t["ensembl_gene_id"], t["pmid"]): t for t in targets}
    panel = {g["ensembl_gene_id"]: g for g in read_tsv(root / "genes.tsv")}
    preds = {}
    for m in MODELS:
        for p in read_tsv(root / f"predictions_{m}.tsv.gz"):
            preds[(m, p["ensembl_gene_id"], p["pmid"], p["condition"])] = float(p["tpm"])
    rng = np.random.default_rng(seed)
    out = {"studies": {}, "models": {}}
    order = sorted(studies, key=lambda p: (studies[p]["tissue"], p))
    pairs = {m: [] for m in MODELS}
    for pmid in order:
        s = studies[pmid]
        keys = sorted(g for (g, p) in ref if p == pmid)
        y = np.array([float(ref[(g, pmid)]["reference_log2_change"]) for g in keys])
        adj = np.array([float(ref[(g, pmid)]["adjusted_p"] or 1) for g in keys])
        entry = dict(tissue=s["tissue"], age_young=int(s["age_young"]), age_old=int(s["age_old"]),
                     n_genes=len(keys), n_significant=int((adj < 0.05).sum()),
                     n_rounded_zero=int((y == 0).sum()), note=s["cohort_disease_note"],
                     curation=s["curation_note"], median_abs_reference=float(np.median(np.abs(y))))
        for m in MODELS:
            young = np.array([preds[(m, g, pmid, "young")] for g in keys])
            old = np.array([preds[(m, g, pmid, "old")] for g in keys])
            x = np.log2(old + 1) - np.log2(young + 1)
            rho = spearman(x, y)
            boot = [spearman(x[i], y[i]) for i in (rng.integers(0, len(x), len(x)) for _ in range(perms))]
            null = np.array([spearman(x, rng.permutation(y)) for _ in range(perms)])
            nz = (y != 0) & (x != 0)
            entry[m] = dict(spearman=rho, ci95=np.nanquantile(boot, [0.025, 0.975]).tolist(),
                            perm_p_two_sided=float((np.sum(np.abs(null) >= abs(rho)) + 1) / (perms + 1)),
                            direction_agreement=float(np.mean(np.sign(x[nz]) == np.sign(y[nz]))),
                            n_direction=int(nz.sum()), median_abs_predicted=float(np.median(np.abs(x))),
                            n_pred_up=int((x > 0).sum()), n_pred_down=int((x < 0).sum()),
                            n_pred_zero=int((x == 0).sum()))
            pairs[m] += [(g, pmid, xi, yi, yo) for g, xi, yi, yo in zip(keys, x, y, young)]
        out["studies"][pmid] = entry
    for m in MODELS:
        rhos = [out["studies"][p][m]["spearman"] for p in order]
        allx = np.array([r[2] for r in pairs[m]])
        ally = np.array([r[3] for r in pairs[m]])
        out["models"][m] = dict(pooled_spearman=spearman(allx, ally), median_study_spearman=float(np.median(rhos)),
                                n_positive=int(sum(r > 0 for r in rhos)),
                                n_ci_above_zero=int(sum(out["studies"][p][m]["ci95"][0] > 0 for p in order)),
                                n_ci_below_zero=int(sum(out["studies"][p][m]["ci95"][1] < 0 for p in order)),
                                n_perm_p_below_005=int(sum(out["studies"][p][m]["perm_p_two_sided"] < 0.05 for p in order)),
                                median_abs_predicted=float(np.median(np.abs(allx))),
                                n_pairs=len(allx), n_requests=2 * len(allx))
    # agreement between the two models
    xa = np.array([r[2] for r in pairs[MODELS[0]]])
    xb = np.array([r[2] for r in pairs[MODELS[1]]])
    la = np.log2(np.array([preds[(MODELS[0], g, p, c)] for (g, p, *_ ) in pairs[MODELS[0]] for c in ("young", "old")]) + 1)
    lb = np.log2(np.array([preds[(MODELS[1], g, p, c)] for (g, p, *_ ) in pairs[MODELS[0]] for c in ("young", "old")]) + 1)
    out["model_agreement"] = dict(
        level_spearman=spearman(la, lb), change_spearman=spearman(xa, xb),
        change_sign_agreement=float(np.mean(np.sign(xa[(xa != 0) & (xb != 0)]) == np.sign(xb[(xa != 0) & (xb != 0)]))),
        study_rho_spearman=spearman([out["studies"][p][MODELS[0]]["spearman"] for p in order],
                                    [out["studies"][p][MODELS[1]]["spearman"] for p in order]))
    # per-gene consistency of predicted age change across studies
    for m in MODELS:
        by_gene = {}
        for g, p, x, *_ in pairs[m]:
            by_gene.setdefault(g, []).append(x)
        cons = [np.mean(np.sign(v) == np.sign(np.mean(v))) for v in by_gene.values() if len(v) >= 4]
        out["models"][m]["gene_sign_consistency_median"] = float(np.median(cons))
    # agreement between experimental references (HCATA effects of the same genes across studies)
    mat = np.full((len(order), len(order)), np.nan)
    for i, a in enumerate(order):
        for j, b in enumerate(order):
            shared = [g for (g, p) in ref if p == a and (g, b) in ref]
            if len(shared) >= 20:
                mat[i, j] = spearman([float(ref[(g, a)]["lfc_log2_per_year_reported"]) for g in shared],
                                     [float(ref[(g, b)]["lfc_log2_per_year_reported"]) for g in shared])
    same = {"lung": ["30554520", "32832599", "36108172"], "colon": ["32888429", "34450029"],
            "brain": ["31178122", "31316211"]}
    same_t = {}
    for t, ps in same.items():
        vals = [mat[order.index(a), order.index(b)] for i, a in enumerate(ps) for b in ps[i + 1:]]
        same_t[t] = [float(v) for v in vals]
    off = mat[~np.eye(len(order), dtype=bool)]
    out["reference_agreement"] = dict(order=order, matrix=mat.tolist(), same_tissue=same_t,
                                      median_all_pairs=float(np.nanmedian(off)))
    # baseline diagnostic: predicted young-age level vs HCATA fitted level at that age
    raw = {}
    for f in sorted(raw_dir.glob("effects_batch_*.json.gz")):
        for r in json.loads(gzip.decompress(f.read_bytes())):
            raw[(f"ENSG{int(r['gene']):011d}", str(r["pmid"]), " ".join(r["cell_type"].split()))] = r
    base = {}
    for pmid in order:
        lab = " ".join(studies[pmid]["hcata_context_label"].split())
        age = int(studies[pmid]["age_young"])
        rows = [(g, raw[(g, pmid, lab)]) for (g, p) in ref if p == pmid and (g, pmid, lab) in raw]
        fitted = np.array([float(r["inter"]) + float(r["slope"]) * age for _, r in rows])
        base[pmid] = {m: spearman([math.log2(preds[(m, g, pmid, "young")] + 1) for g, _ in rows], fitted)
                      for m in MODELS}
    out["baseline_level_spearman"] = base
    out["panel"] = dict(n_genes=len(panel), groups={k: sum(1 for g in panel.values() if g.get("panel_group") == k)
                                                     for k in sorted({g.get("panel_group") for g in panel.values()})})
    out["_pairs"] = {m: [(g, p, float(x), float(y)) for g, p, x, y, _ in pairs[m]] for m in MODELS}
    return out


# ---------------------------------------------------------------- figures


def figures(res: dict, out: Path) -> list[Path]:
    order = res["reference_agreement"]["order"]
    st = res["studies"]
    labels = [f"{st[p]['tissue']} ({p})" for p in order]
    paths = []

    # Figure 1: per-study Spearman with bootstrap CI, both models
    f, ax = plt.subplots(figsize=(7, 3.6), dpi=200)
    f.patch.set_facecolor(SURFACE)
    yy = np.arange(len(order))[::-1]
    for k, m in enumerate(MODELS):
        off = 0.17 if k == 0 else -0.17
        r = [st[p][m]["spearman"] for p in order]
        lo = [st[p][m]["ci95"][0] for p in order]
        hi = [st[p][m]["ci95"][1] for p in order]
        ax.hlines(yy + off, lo, hi, color=COLOR[m], linewidth=1.6)
        ax.plot(r, yy + off, "o", color=COLOR[m], markersize=4.5, markeredgecolor=SURFACE,
                markeredgewidth=1, label=LABEL[m])
    ax.axvline(0, color=INK2, linewidth=0.9)
    ax.set_yticks(yy)
    ax.set_yticklabels(labels)
    style(ax, "Predicted vs HCATA age change, per study", "Spearman ρ across genes (95% gene bootstrap)", "")
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc="upper center", bbox_to_anchor=(0.45, -0.16), ncol=2)
    f.tight_layout()
    paths.append(out / "fig1_study_correlations.png")
    f.savefig(paths[-1], facecolor=SURFACE)

    # Figure 2: reference agreement matrix (left), baseline level diagnostic (right)
    f, (a1, a2) = plt.subplots(1, 2, figsize=(7, 3.3), dpi=200, gridspec_kw=dict(width_ratios=[1.15, 1]))
    f.patch.set_facecolor(SURFACE)
    mat = np.array(res["reference_agreement"]["matrix"], dtype=float)
    np.fill_diagonal(mat, np.nan)  # a study against itself is trivially 1
    cmap = DIVERGING.copy()
    cmap.set_bad("#d9d8d4")
    im = a1.imshow(np.ma.masked_invalid(mat), cmap=cmap, vmin=-0.6, vmax=0.6)
    tshort = {"skeletal muscle": "muscle", "bone marrow": "marrow"}
    short = [f"{tshort.get(st[p]['tissue'], st[p]['tissue'])} {p[-3:]}" for p in order]
    a1.set_xticks(range(len(order)))
    a1.set_xticklabels(short, rotation=90, fontsize=6)
    a1.set_yticks(range(len(order)))
    a1.set_yticklabels(short, fontsize=6)
    a1.set_title("HCATA age slopes: agreement between studies", loc="left", color=INK, fontsize=9)
    a1.tick_params(colors=INK2, length=0)
    for s in a1.spines.values():
        s.set_visible(False)
    cb = f.colorbar(im, ax=a1, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=6.5, colors=INK2)
    cb.outline.set_visible(False)
    cb.set_label("Spearman ρ (shared genes)", fontsize=7, color=INK2)
    base = res["baseline_level_spearman"]
    yy = np.arange(len(order))[::-1]
    for k, m in enumerate(MODELS):
        off = 0.15 if k == 0 else -0.15
        a2.plot([base[p][m] for p in order], yy + off, "o", color=COLOR[m], markersize=4,
                markeredgecolor=SURFACE, markeredgewidth=1, label=m)
    a2.axvline(0, color=INK2, linewidth=0.9)
    a2.set_yticks(yy)
    a2.set_yticklabels(short, fontsize=6.5)
    style(a2, "Expression level at the young age", "Spearman ρ, GI level vs HCATA fit", "")
    a2.legend(frameon=False, fontsize=6.5, labelcolor=INK2, loc="upper left")
    f.tight_layout()
    paths.append(out / "fig3_reference_checks.png")
    f.savefig(paths[-1], facecolor=SURFACE)

    # Figure 3: agreement of predicted age changes between models (left), effect sizes (right)
    f, (a1, a2) = plt.subplots(1, 2, figsize=(7, 3.1), dpi=200)
    f.patch.set_facecolor(SURFACE)
    xa = np.array([r[2] for r in res["_pairs"][MODELS[0]]])
    xb = np.array([r[2] for r in res["_pairs"][MODELS[1]]])
    a1.scatter(xa, xb, s=5, color=INK2, alpha=0.45, linewidths=0)
    lim = float(np.quantile(np.abs(np.r_[xa, xb]), 0.995)) * 1.15
    n_out = int(np.sum((np.abs(xa) > lim) | (np.abs(xb) > lim)))
    a1.plot([-lim, lim], [-lim, lim], color=INK2, linewidth=0.8)
    if n_out:
        a1.annotate(f"{n_out} pairs outside the axes", (0.02, 0.96), xycoords="axes fraction",
                    fontsize=6.5, color=INK2, va="top")
    a1.axhline(0, color=GRID, linewidth=0.8)
    a1.axvline(0, color=GRID, linewidth=0.8)
    a1.set_xlim(-lim, lim)
    a1.set_ylim(-lim, lim)
    style(a1, f"Predicted age change, the two models (ρ {res['model_agreement']['change_spearman']:.2f})",
          "g0-expression, log2 change", "g0-expression-8192, log2 change")
    yy = np.arange(len(order))[::-1]
    a2.plot([st[p]["median_abs_reference"] for p in order], yy, "s", color=INK, markersize=4,
            label="HCATA reference")
    for k, m in enumerate(MODELS):
        a2.plot([st[p][m]["median_abs_predicted"] for p in order], yy + (0.18 if k == 0 else -0.18), "o",
                color=COLOR[m], markersize=4, markeredgecolor=SURFACE, markeredgewidth=1, label=m)
    a2.set_yticks(yy)
    a2.set_yticklabels(short, fontsize=6.5)
    style(a2, "Median absolute age change per study", "|log2 change|, young to old age", "")
    a2.legend(frameon=False, fontsize=6.5, labelcolor=INK2, loc="lower right")
    f.tight_layout()
    paths.append(out / "fig2_model_agreement.png")
    f.savefig(paths[-1], facecolor=SURFACE)
    return paths


# ---------------------------------------------------------------- docx


def build_docx(res: dict, figs: list[Path], template: Path, out: Path, meta: dict) -> None:
    import docx
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Inches, Pt

    d = docx.Document(template)
    body = d.element.body
    for el in list(body):
        if el.tag != qn("w:sectPr"):
            body.remove(el)
    P = lambda text, st="Normal": d.add_paragraph(text, style=st)  # noqa: E731

    def table(rows, widths=None, bold_first_col=False):
        t = d.add_table(rows=len(rows), cols=len(rows[0]))
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        for i, row in enumerate(rows):
            for j, val in enumerate(row):
                cell = t.cell(i, j)
                cell.text = ""
                run = cell.paragraphs[0].add_run(str(val))
                run.font.size = Pt(9)
                run.bold = i == 0 or (bold_first_col and j == 0)
                tcPr = cell._tc.get_or_add_tcPr()
                borders = OxmlElement("w:tcBorders")
                for side in ("top", "bottom"):
                    b = OxmlElement(f"w:{side}")
                    b.set(qn("w:val"), "single")
                    b.set(qn("w:sz"), "4")
                    b.set(qn("w:color"), "BFBFBF")
                    borders.append(b)
                tcPr.append(borders)
                if i == 0:
                    shd = OxmlElement("w:shd")
                    shd.set(qn("w:val"), "clear")
                    shd.set(qn("w:fill"), "EDEDEA")
                    tcPr.append(shd)
                if widths:
                    cell.width = Inches(widths[j])
        d.add_paragraph()
        return t

    def figure(path, caption):
        d.add_paragraph().add_run().add_picture(str(path), width=Inches(7.0))
        P(caption, "Caption")

    st, mo, ag, ra = res["studies"], res["models"], res["model_agreement"], res["reference_agreement"]
    order = ra["order"]
    a, b = MODELS
    f3 = lambda v: f"{v:.3f}"  # noqa: E731
    f2 = lambda v: f"{v:.2f}"  # noqa: E731
    fp = lambda v: "< 0.001" if v < 0.001 else f"{v:.3f}"  # noqa: E731

    P("GI Expression Predictions of Age-Related Changes Across 12 Human Studies (HCATA)", "Title")
    P(f"150-gene panel, two GI expression models  |  {meta['date']}", "Subtitle")
    P(f"Neither GI expression model reproduced HCATA's age-related expression changes. Across the "
      f"12 studies, the correlation of predicted with reported age changes was "
      f"{f3(mo[a]['pooled_spearman'])} for g0-expression and {f3(mo[b]['pooled_spearman'])} for "
      f"g0-expression-8192 (all gene-study pairs pooled). Median per-study correlations were "
      f"{f3(mo[a]['median_study_spearman'])} and {f3(mo[b]['median_study_spearman'])}.")
    P(f"Some studies had bootstrap intervals excluding zero, but in both directions: "
      f"{mo[a]['n_ci_above_zero']} above and {mo[a]['n_ci_below_zero']} below zero for g0-expression, "
      f"{mo[b]['n_ci_above_zero']} above and {mo[b]['n_ci_below_zero']} below for g0-expression-8192. "
      f"The two models agreed on expression levels (ρ {f2(ag['level_spearman'])}) but only weakly on "
      f"predicted age changes (ρ {f2(ag['change_spearman'])}), so individual study results often changed "
      f"sign between models. The HCATA references are themselves weak: 0–5 panel genes per study pass "
      f"adjusted p < 0.05 (21 in the pancreas study), and the same-tissue studies agree poorly with each "
      f"other.")
    P("This study extends the earlier 150-gene muscle pilot. It used the same frozen gene panel, but scored "
      "it against HCATA's single-cell age effects in 12 human studies covering 8 tissues, with two GI "
      "expression models. Only the donor age and tissue in the text description change between paired "
      "requests; the DNA input of each gene is fixed.")
    figure(figs[0], "Figure 1. Each row is an HCATA study (tissue and PubMed ID). The point is the Spearman "
           "correlation, across panel genes, between the GI-predicted log2(TPM+1) change from the study's "
           "young to old age and HCATA's annual log2 fold change multiplied by the same age gap. Lines are "
           "95% intervals from 2,000 gene bootstrap resamples.")
    P(f"All {mo[a]['n_requests']:,} requests per model ({mo[a]['n_pairs']:,} gene-study pairs) completed "
      "with no failed predictions. This study tests age-conditioned text prompts; it does not validate "
      "a biological age clock, an intervention or a causal mechanism.", "Caption")

    P("How the experiment was performed", "Heading 1")
    P("Experimental reference data", "Heading 2")
    P("The Human Cell Aging Transcriptome Atlas (HCATA) [1] fits, for every gene and cell-type context of "
      "each included single-cell study, a regression of normalised expression on donor age. Its public "
      "API reports the annual log2 fold change (lfc), the fitted slope and intercept, and a "
      "−log10 adjusted p value. The repository already held HCATA's full responses for the 150 panel "
      "genes, retrieved for the muscle benchmark (data/hcata_muscle_v1/source). These responses include "
      "every human study that reports these genes: 12 studies and 28,818 gene-context effects.")
    P("For each study, the primary target was its combined 'All Cells' context. The reference age change "
      "for a gene was lfc × (old age − young age). Effects were never filtered by significance, and "
      "rounded-zero coefficients were kept. Two HCATA context labels are swapped relative to the "
      "studies: PMID 36108172 is labelled 'Liver' but is a COPD lung study, and PMID 35021063 is "
      "labelled 'Lung2' but is a hepatic macrophage study. The PubMed titles and HCATA's own sample "
      "table agree on lung and liver respectively, and these were used. PMID 30554520 has an empty "
      "tissue label (lung).")
    rows = [["PMID", "Tissue", "Ages used", "Donor ages", "Genes", "Adj. p < 0.05", "Cohort note"]]
    for p in order:
        s = st[p]
        rows.append([p, s["tissue"], f"{s['age_young']} / {s['age_old']}",
                     meta["donor_ranges"][p], s["n_genes"], s["n_significant"], s["note"]])
    table(rows, widths=[0.75, 0.95, 0.7, 0.7, 0.5, 0.65, 2.25])
    P("The young and old ages are each study's 10th and 90th percentile donor ages in HCATA's sample "
      "table, so both fall inside the observed cohort. Several cohorts include disease samples (COPD, "
      "interstitial lung disease, colorectal cancer, multiple sclerosis, obesity/type 2 diabetes, head "
      "and neck cancer). HCATA's regressions pool them, and the GI descriptions do not mention disease.")
    P("Gene panel and DNA inputs", "Heading 2")
    P(f"The panel is the frozen 150-gene panel of the muscle benchmark: 50 HAGR age-up, 50 HAGR age-down "
      f"and 50 expression-matched background genes. Of these, 148 have HCATA records; per study, "
      f"122–144 genes have an 'All Cells' effect. No genes were added or removed after seeing results.")
    P("Sequences came from Ensembl release 115 (GRCh38 primary assembly). For each gene, the TSS of the "
      "Ensembl canonical transcript was taken, and 40,960 bases on each side were extracted on the "
      "gene's sense strand (81,921 bp; at most 1.2% N). Each request sent this sequence with tss_index "
      "40,960, and the API cut each model's own window around it: 9,198 bp for g0-expression and up to "
      "82 kb for g0-expression-8192. The same sequence was used for every condition of a gene. Its "
      "SHA-256 is recorded in sequences.tsv.gz.")
    P("Exact descriptions", "Heading 2")
    P("Template: Human <tissue phrase> from a(n) <age>-year-old donor.")
    P("Example (skeletal muscle): Human skeletal muscle tissue from a 20-year-old donor. / … from an "
      "82-year-old donor.")
    P("The tissue phrases were lung tissue, bone marrow, liver tissue, skeletal muscle tissue, colon "
      "tissue, pancreatic tissue, brain tissue and blood. The word 'adult' was dropped because the "
      "youngest pancreas age is 14. Each gene-study pair received one young and one old request per "
      "model. No baseline or alternative-wording controls were run in this study.")
    P("Models and execution", "Heading 2")
    P(f"Predictions used the Genomic Intelligence REST API [2] (POST /v1/tasks/expression/predict). "
      f"g0-expression-8192 ran first ({meta['runs'][b]}), then g0-expression ({meta['runs'][a]}), each with "
      f"3,104 requests and 0 failures. The model id returned in every response matched the requested "
      f"model. Requests ran 16 at a time under a rate limiter that follows the API's RateLimit headers. "
      f"The GI response field expression_tpm was used as TPM.")

    P("What the quantitative results mean", "Heading 1")
    P("The primary statistic is the Spearman correlation across genes, within a study, between the "
      "predicted change log2(TPM_old + 1) − log2(TPM_young + 1) and the reference change. It asks "
      "whether genes predicted to rise most with age are those HCATA reports rising most. It does not "
      "test absolute TPM calibration: HCATA values are normalised single-nucleus or single-cell counts, "
      "not TPM. The intervals resample genes, not donors. The permutation P value shuffles the "
      "reference across genes within a study (2,000 permutations, two-sided). All P values are "
      "exploratory and unadjusted for the 24 tests.")
    rows = [["PMID", "Tissue", "n", "ρ g0-expression [95% CI]", "P", "ρ g0-expression-8192 [95% CI]", "P"]]
    for p in order:
        s = st[p]
        rows.append([p, s["tissue"], s["n_genes"],
                     f"{f2(s[a]['spearman'])} [{f2(s[a]['ci95'][0])}, {f2(s[a]['ci95'][1])}]",
                     fp(s[a]['perm_p_two_sided']),
                     f"{f2(s[b]['spearman'])} [{f2(s[b]['ci95'][0])}, {f2(s[b]['ci95'][1])}]",
                     fp(s[b]['perm_p_two_sided'])])
    rows.append(["All", "pooled", mo[a]["n_pairs"], f2(mo[a]["pooled_spearman"]), "",
                 f2(mo[b]["pooled_spearman"]), ""])
    rows.append(["", "median study", "", f2(mo[a]["median_study_spearman"]), "",
                 f2(mo[b]["median_study_spearman"]), ""])
    table(rows, widths=[0.75, 1.0, 0.45, 1.55, 0.5, 1.75, 0.5])
    def ex(m):  # studies whose interval excludes zero, as "tissue PMID (rho)"
        return ", ".join(f"{st[p]['tissue']} {p} ({f2(st[p][m]['spearman'])})" for p in order
                         if st[p][m]["ci95"][0] > 0 or st[p][m]["ci95"][1] < 0)
    P(f"{mo[a]['n_positive']} of 12 study correlations were positive for g0-expression and "
      f"{mo[b]['n_positive']} for g0-expression-8192. Intervals excluding zero: g0-expression {ex(a)}; "
      f"g0-expression-8192 {ex(b)}. Chance alone would give about 0.6 such studies per model, so "
      f"g0-expression shows more than expected. However, these results go in both directions and "
      f"mostly do not replicate in the other model. "
      f"The only results pointing the same way in both models are negative. The multiple sclerosis "
      f"brain study ({f2(st['31316211'][a]['spearman'])} / {f2(st['31316211'][b]['spearman'])}) excludes zero "
      f"in both. The bone-marrow study ({f2(st['30518681'][a]['spearman'])} / "
      f"{f2(st['30518681'][b]['spearman'])}) excludes zero only for g0-expression. In both, predictions run "
      f"opposite to HCATA.")
    rows = [["Additional measure", "g0-expression", "g0-expression-8192"],
            ["Direction agreement, median over studies",
             f"{np.median([st[p][a]['direction_agreement'] for p in order]):.1%}",
             f"{np.median([st[p][b]['direction_agreement'] for p in order]):.1%}"],
            ["Studies with permutation P < 0.05", mo[a]["n_perm_p_below_005"], mo[b]["n_perm_p_below_005"]],
            ["Median |predicted age change|, log2", f3(mo[a]["median_abs_predicted"]), f3(mo[b]["median_abs_predicted"])],
            ["Per-gene sign consistency across studies (median)", f"{mo[a]['gene_sign_consistency_median']:.0%}",
             f"{mo[b]['gene_sign_consistency_median']:.0%}"]]
    table(rows, widths=[3.3, 1.7, 1.7])
    P(f"Direction agreement counts genes where both the predicted and reported changes are nonzero. "
      f"The median |reference change| across studies was "
      f"{np.median([st[p]['median_abs_reference'] for p in order]):.3f} log2 units, comparable in size to "
      f"the predicted changes. For a given gene, the predicted direction of change is the same in "
      f"{mo[a]['gene_sign_consistency_median']:.0%} (g0-expression) and "
      f"{mo[b]['gene_sign_consistency_median']:.0%} (g0-expression-8192) of studies at the median. The models therefore give each gene a partly consistent 'age response', "
      f"but it does not line up with HCATA's.")

    P("Agreement between the two models", "Heading 2")
    P(f"The two models agree closely on expression levels: Spearman {f2(ag['level_spearman'])} over all "
      f"6,208 predicted values. They agree only modestly on predicted age changes: Spearman "
      f"{f2(ag['change_spearman'])}, with the same sign for {ag['change_sign_agreement']:.0%} of pairs. Their "
      f"per-study correlations with HCATA are only weakly related (Spearman "
      f"{f2(ag['study_rho_spearman'])} across the 12 studies). For example, lung study 32832599 gives "
      f"{f2(st['32832599'][b]['spearman'])} with g0-expression-8192 and {f2(st['32832599'][a]['spearman'])} "
      f"with g0-expression, while liver changes from {f2(st['35021063'][b]['spearman'])} to "
      f"{f2(st['35021063'][a]['spearman'])}.")
    figure(figs[2], "Figure 2. Left: predicted age changes for all 1,552 gene-study pairs from the two models; "
           "the diagonal marks identical predictions. Right: per study, the median absolute predicted age "
           "change for each model and the median absolute HCATA reference change.")

    P("Checking the experimental references", "Heading 1")
    P("Agreement between HCATA studies", "Heading 2")
    same = ra["same_tissue"]
    P(f"The age slopes in HCATA agree little between studies of the same tissue. For the panel genes, "
      f"the three lung studies correlate at {', '.join(f2(v) for v in same['lung'])}, the two colon studies "
      f"at {f2(same['colon'][0])} and the two brain studies at {f2(same['brain'][0])}. The median over all "
      f"study pairs is {f2(ra['median_all_pairs'])}. A model that predicted one tissue's ageing exactly could "
      f"not correlate highly with another study of the same tissue. Weak agreement with any single "
      f"reference therefore cannot by itself diagnose model failure.")
    figure(figs[1], "Figure 3. Left: Spearman correlation of HCATA annual log2 fold changes between studies, "
           "over shared panel genes ('All Cells' contexts). Right: Spearman correlation between the GI "
           "predicted expression level at each study's young age and HCATA's fitted expression at that "
           "age (intercept + slope × age). The right panel is a baseline expression diagnostic, not an "
           "age test.")
    base = res["baseline_level_spearman"]
    P(f"The models' predicted levels track HCATA's fitted expression levels moderately: median Spearman "
      f"{f2(np.median([base[p][a] for p in order]))} (g0-expression) and "
      f"{f2(np.median([base[p][b] for p in order]))} (g0-expression-8192) across studies. The models "
      f"therefore capture which genes are highly or lowly expressed in these tissues, which makes the "
      f"absence of age agreement a statement about the age response rather than about expression levels "
      f"in general.")

    P("Interpretation and limits", "Heading 1")
    P("This study did not find a reproducible positive association between GI's age-conditioned "
      "predictions and HCATA's reported age effects. This holds for both models, and for every tissue "
      "with more than one study. The negative associations in the multiple sclerosis brain study and "
      "in bone marrow are not explained, and both cohorts have unusual composition.")
    P("What this study can support", "Heading 2")
    P("The workflow scores a fixed gene panel against 12 independent single-cell age references with two "
      "models, reproducibly and without failures. It also reports reference consistency and expression "
      "level diagnostics alongside the age comparison. The same code can extend the study to HCATA's "
      "full gene universe (29,166 genes). Its effects crawl was tested on 50 genes; the full crawl "
      "would take an estimated 30–40 hours against HCATA's public API.")
    P("What remains unresolved", "Heading 2")
    P("Reference strength. HCATA's age effects are mostly non-significant point estimates from small "
      "cohorts with disease samples and uncertain cell composition. The same-tissue studies disagree "
      "with one another. Intervals resample genes, not donors.")
    P("Prompt specificity. No same-age wording controls were run here. In the earlier muscle pilot, "
      "rewording an age-25 description shifted predictions more than the age contrast did. Any age "
      "response measured here is therefore of uncertain specificity.")
    P("Training independence. GI's catalogue lists ENCODE, GTEx/recount3 and cellxgene pseudobulk among "
      "its training data. Overlap with the HCATA source studies has not been checked.")
    P("Scope. The text changes while the reference DNA stays fixed, so this probes context conditioning, "
      "not age-related genomic change. Only the combined 'All Cells' context was scored; HCATA's "
      "cell-type effects were not used.")
    P("Recommended next experiment", "Heading 2")
    P("Score cell-type contexts with matching cell-type descriptions, and add same-age wording and "
      "irrelevant-context controls before inference. Restrict the primary comparison to genes with "
      "precise and cross-study consistent age effects, and verify that the evaluation cohorts are "
      "excluded from GI's training data.")

    P("Reproducibility and sources", "Heading 1")
    P("Everything is in this repository. scripts/hcata_full.py builds the tables (genes, effect-tables, "
      "sequences, requests), runs the predictions (predict, cached per request) and computes the "
      "per-study evaluation (evaluate). scripts/hcata_panel_report.py produces this report, its figures "
      "and stats JSON. data/hcata_panel150_v1 holds the panel, study curation (studies.tsv), primary "
      "targets, sequence metadata and hashes, requests, predictions for both models (with GI request "
      "IDs) and per-study pair tables. Raw GI responses are cached outside the repository.")
    P(f"Panel: 150 genes; 1,552 gene-study pairs; 3,104 requests per model; 0 failed requests. "
      f"Bootstrap and permutation seed {meta['seed']}, {meta['perms']} iterations each.", "Caption")
    P("References", "Heading 2")
    for r in meta["references"]:
        P(r)
    P(f"Sources accessed 3–4 October 2026. Data and model contributions are credited to their providers. "
      f"This computational analysis includes no new human sampling or laboratory experiments.", "Caption")
    d.save(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("data/hcata_panel150_v1"))
    ap.add_argument("--raw", type=Path, default=Path("data/hcata_muscle_v1/source"))
    ap.add_argument("--template", type=Path, default=Path("reports/gi_age_expression_150_gene_report.docx"))
    ap.add_argument("--out", type=Path, default=Path("reports"))
    ap.add_argument("--perms", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261004)
    a = ap.parse_args()
    fig_dir = a.out / "gi_hcata_12_study_figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    res = analyse(a.root, a.raw, a.perms, a.seed)
    figs = figures(res, fig_dir)
    samples = json.loads(gzip.decompress((a.raw / "hcata_samples.json.gz").read_bytes()))
    donor = {}
    for p in res["studies"]:
        ages = [float(s["age"]) for s in samples if str(s.get("study_id")) == p
                and str(s.get("age", "")).replace(".", "", 1).isdigit()]
        donor[p] = f"{int(min(ages))}–{int(max(ages))}" if ages else "?"
    meta = dict(date="4 October 2026", seed=a.seed, perms=a.perms, donor_ranges=donor,
                runs={"g0-expression-8192": "2026-10-03 19:17–19:35 UTC",
                      "g0-expression": "2026-10-03 19:41–19:44 UTC"},
                references=[
                    "[1] Bartz J et al. Human Cell Aging Transcriptome Atlas (HCATA): a single-cell atlas of "
                    "age-associated transcriptomic alterations across human tissues. Communications Biology "
                    "(2025). https://doi.org/10.1038/s42003-025-08845-8  Portal: http://hcata-xiaodonglab.org:3304",
                    "[2] Genomic Intelligence. REST API contract (OpenAPI 2026.10.02.6) and expression model "
                    "catalogue. https://api.genomicintelligence.ai/v1/openapi.json",
                    "[3] Ensembl release 115, GRCh38 primary assembly and gene annotation. https://ftp.ensembl.org/pub/release-115/",
                    "[4] Human Ageing Genomic Resources. Ageing expression signatures (panel selection). "
                    "https://www.genomics.senescence.info/genes/microarray.php",
                    "[5] Source studies (HCATA 'All Cells' contexts): PubMed IDs 36108172, 32832599, 30554520, "
                    "30518681, 35021063, 36516485, 34450029, 34428183, 32888429, 31316211, 31924475, 31178122. "
                    "https://pubmed.ncbi.nlm.nih.gov/",
                    "[6] Earlier report in this repository: GI expression predictions of age-related changes in "
                    "human muscle, preliminary 150-gene study (reports/gi_age_expression_150_gene_report.docx)."])
    out_docx = a.out / "gi_hcata_12_study_150_gene_report.docx"
    build_docx(res, figs, a.template, out_docx, meta)
    stats = {k: v for k, v in res.items() if not k.startswith("_")}
    (fig_dir / "stats.json").write_text(json.dumps(stats, indent=1, default=float))
    print(f"wrote {out_docx} and {fig_dir}")
    for m in MODELS:
        print(m, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in res["models"][m].items()})
    print("agreement", {k: round(v, 3) for k, v in res["model_agreement"].items()})
    print("same-tissue reference", res["reference_agreement"]["same_tissue"], "median",
          round(res["reference_agreement"]["median_all_pairs"], 3))
    print("baseline", {p: {m: round(v, 2) for m, v in d.items()} for p, d in res["baseline_level_spearman"].items()})


if __name__ == "__main__":
    main()
