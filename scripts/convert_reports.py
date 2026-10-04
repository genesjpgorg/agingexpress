#!/usr/bin/env python3
"""Convert the DOCX reports in reports/ to Markdown and PDF.

Markdown: pandoc -> GitHub-flavoured Markdown, images extracted to reports/<name>_media/.
PDF: pandoc -> standalone HTML (images embedded) styled like the DOCX (Letter, 0.75 in margins,
Calibri-metric Carlito font if --fonts is given) -> WeasyPrint.

Requires pypandoc_binary and weasyprint (pip). The DOCX files stay the source of truth.

Usage: python scripts/convert_reports.py [--fonts DIR_WITH_CARLITO_TTF] [reports/*.docx]
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import pypandoc
import weasyprint

CSS = """
@page { size: Letter; margin: 0.67in 0.75in 0.75in 0.75in;
        @bottom-center { content: counter(page); font-size: 8.5pt; color: #666; } }
body { font-family: {font}; font-size: 10.5pt; line-height: 1.32; color: #000; }
h1.title { font-size: 22pt; margin: 0 0 4pt 0; }
p.subtitle { font-size: 11pt; color: #444; margin: 0 0 12pt 0; }
h1 { font-size: 15pt; margin: 16pt 0 6pt 0; }
h2 { font-size: 12pt; margin: 12pt 0 4pt 0; }
p { margin: 0 0 6pt 0; }
img { max-width: 100%; display: block; margin: 6pt auto 2pt auto; }
table { border-collapse: collapse; margin: 6pt auto 10pt auto; font-size: 8.8pt; }
tr { break-inside: avoid; }
td { hyphens: manual; }
th, td { border-top: 0.5pt solid #bfbfbf; border-bottom: 0.5pt solid #bfbfbf; padding: 2.5pt 6pt;
         vertical-align: top; text-align: left; }
th, thead td { background: #ededea; font-weight: bold; }
.caption, p.caption { font-size: 9pt; color: #333; }
a { color: #1f5fa8; text-decoration: none; word-break: break-all; }
"""


def font_faces(fonts: Path | None) -> tuple[str, str]:
    if not fonts:
        return "", "'DejaVu Sans', sans-serif"
    faces = []
    for style, weight, italic in (("Regular", 400, "normal"), ("Bold", 700, "normal"),
                                  ("Italic", 400, "italic"), ("BoldItalic", 700, "italic")):
        f = fonts / f"Carlito-{style}.ttf"
        if f.exists():
            faces.append(f"@font-face {{ font-family: 'Carlito'; src: url('{f.resolve().as_uri()}'); "
                         f"font-weight: {weight}; font-style: {italic}; }}")
    return "\n".join(faces), "'Carlito', 'DejaVu Sans', sans-serif"


def mark_captions(html: str, docx: Path) -> str:
    """pandoc drops Word paragraph styles; re-tag Caption-style paragraphs by their text."""
    import docx as docxlib

    caps = [p.text.strip() for p in docxlib.Document(docx).paragraphs
            if p.style.name == "Caption" and p.text.strip()]
    sub = [p.text.strip() for p in docxlib.Document(docx).paragraphs if p.style.name == "Subtitle"]
    for text in caps:
        key = re.escape(text[:40].replace("&", "&amp;"))
        html = re.sub(rf"<p>(?=({key}))", '<p class="caption">', html, count=1)
    for text in sub:
        key = re.escape(text[:30].replace("&", "&amp;"))
        html = re.sub(rf"<p>(?=({key}))", '<p class="subtitle">', html, count=1)
    return html


def convert(docx: Path, fonts: Path | None) -> None:
    stem, out = docx.stem, docx.parent
    media = out / f"{stem}_media"
    if media.exists():
        shutil.rmtree(media)
    # Markdown (images extracted next to it, paths relative to reports/)
    md = pypandoc.convert_file(str(docx), "gfm", extra_args=["--wrap=none", f"--extract-media={media}"])
    md = md.replace(str(media) + "/", f"{stem}_media/").replace(str(media.resolve()) + "/", f"{stem}_media/")
    for p in sorted(media.rglob("*")) if media.exists() else []:
        if p.is_file() and p.parent != media:  # pandoc nests media/<sub>/; flatten
            p.rename(media / p.name)
            md = md.replace(f"{stem}_media/{p.parent.name}/{p.name}", f"{stem}_media/{p.name}")
    for d in sorted(media.glob("*"), reverse=True) if media.exists() else []:
        if d.is_dir():
            shutil.rmtree(d)
    import docx as docxlib  # Word Title/Subtitle become pandoc metadata, not body text

    paras = docxlib.Document(docx).paragraphs
    title = next((x.text.strip() for x in paras if x.style.name == "Title"), "")
    subtitle = next((x.text.strip() for x in paras if x.style.name == "Subtitle"), "")
    head = (f"# {title}\n\n" if title else "") + (f"*{subtitle}*\n\n" if subtitle else "")
    head += f"*Converted from [{docx.name}]({docx.name}); the DOCX is the source of record.*\n\n"
    (out / f"{stem}.md").write_text(head + md)
    # PDF
    faces, family = font_faces(fonts)
    html = pypandoc.convert_file(str(docx), "html5", extra_args=["--standalone", "--embed-resources",
                                                                 "--metadata", f"pagetitle={stem}"])
    html = re.sub(r"<style>.*?</style>", "", html, flags=re.S)  # drop pandoc's default css
    html = mark_captions(html, docx)
    html = re.sub(r"<colgroup>.*?</colgroup>", "", html, flags=re.S)  # size columns by content
    # first H1 is the document title (Word Title style)
    html = re.sub(r"<h1([^>]*)>", r'<h1 class="title"\1>', html, count=1)
    css = faces + CSS.replace("{font}", family)
    from weasyprint.text.fonts import FontConfiguration

    fc = FontConfiguration()  # needed for @font-face (Carlito) to take effect
    weasyprint.HTML(string=html, base_url=str(out.resolve())).write_pdf(
        out / f"{stem}.pdf", stylesheets=[weasyprint.CSS(string=css, font_config=fc)], font_config=fc)
    print(f"{docx} -> {stem}.md, {stem}.pdf")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("docx", nargs="*", type=Path)
    ap.add_argument("--fonts", type=Path, default=None)
    a = ap.parse_args()
    for d in a.docx or sorted(Path("reports").glob("*.docx")):
        convert(d, a.fonts)


if __name__ == "__main__":
    main()
