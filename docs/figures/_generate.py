#!/usr/bin/env python3
"""Regenerate the primer SVGs in this directory from the pilot numbers
embedded below. Stdlib only.  python docs/figures/_generate.py"""
from pathlib import Path

OUT = Path(__file__).resolve().parent
OUT.mkdir(parents=True, exist_ok=True)

INK, DIM = "#1c231d", "#5c6657"
SYN, NON, POOL = "#2f6b4f", "#a8432e", "#31607d"
SANS = ("font-family='ui-sans-serif,-apple-system,Segoe UI,Roboto,Helvetica,"
        "Arial,sans-serif'")
MONO = "font-family='ui-monospace,SFMono-Regular,Menlo,Consolas,monospace'"

# real pilot data (Stage 7, two-ratio foreground branch):
# 20 bins over [0,1) then one overflow bin for omega >= 1
AHAL = [37,75,82,57,50,21,9,8,5,6,3,2,3,1,1,0,2,1,0,0,6]
CGRA = [99,33,35,30,20,10,6,5,4,4,1,2,1,4,4,1,0,1,1,0,108]
ST = {"a_halleri":    dict(mean="5.6", median=0.139, pooled=0.163, n=369),
      "c_grandiflora":dict(mean="226", median=0.177, pooled=0.151, n=369)}


def hist_panel(bins, name, x0, y0, w, h):
    d = ST[name]
    mx = max(bins)
    nb = len(bins)
    bw = w / nb
    s = [f"<text x='{x0}' y='{y0-9}' {SANS} font-size='12.5' font-weight='600' "
         f"fill='{INK}'>{name}  &#183;  {d['n']} genes</text>"]
    s.append(f"<line x1='{x0}' y1='{y0+h}' x2='{x0+w}' y2='{y0+h}' "
             f"stroke='{DIM}' stroke-width='1'/>")
    for i, c in enumerate(bins):
        bh = (c / mx) * h if mx else 0
        fill = "#cbd6cb" if i < nb - 1 else "#e3b7a8"
        s.append(f"<rect x='{x0+i*bw:.1f}' y='{y0+h-bh:.1f}' width='{bw-1.2:.1f}' "
                 f"height='{bh:.1f}' fill='{fill}'/>")
    span = w - bw            # x-range that maps [0,1)
    for frac, lab in [(0, "0"), (0.25, ".25"), (0.5, ".5"), (0.75, ".75")]:
        xx = x0 + frac * span
        s.append(f"<text x='{xx:.0f}' y='{y0+h+15}' {MONO} font-size='9.5' "
                 f"fill='{DIM}' text-anchor='middle'>{lab}</text>")
    s.append(f"<text x='{x0+w-bw/2:.0f}' y='{y0+h+15}' {MONO} font-size='9.5' "
             f"fill='{NON}' text-anchor='middle'>&#969;&#8805;1</text>")
    # dashed markers, labelled compactly above the panel
    for val, col in [(d['median'], INK), (d['pooled'], POOL)]:
        xx = x0 + min(val, 1.0) * span
        s.append(f"<line x1='{xx:.1f}' y1='{y0}' x2='{xx:.1f}' y2='{y0+h}' "
                 f"stroke='{col}' stroke-width='1.4' stroke-dasharray='3 2'/>")
    s.append(f"<text x='{x0+w:.0f}' y='{y0-9}' {MONO} font-size='10' fill='{DIM}' "
             f"text-anchor='end'>median {d['median']:.2f}  &#183;  pooled "
             f"{d['pooled']:.3f}  &#183;  mean {d['mean']} (off scale)</text>")
    return "\n".join(s)


def fig_estimators():
    W, H = 720, 400
    s = [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {W} {H}'>",
         f"<rect width='{W}' height='{H}' fill='#ffffff'/>",
         f"<text x='28' y='28' {SANS} font-size='15' font-weight='700' fill='{INK}'>"
         f"Per-gene &#969; on one lineage's branch &#8212; what each summary sees</text>",
         f"<text x='28' y='47' {SANS} font-size='11.5' fill='{DIM}'>"
         f"Pilot, 369 BUSCO genes. Bar height = number of genes in that &#969; bin. "
         f"Dashed lines: <tspan fill='{INK}'>median</tspan>, "
         f"<tspan fill='{POOL}'>count-pooled &#969;</tspan>.</text>"]
    s.append(hist_panel(AHAL, "a_halleri", 44, 92, 632, 110))
    s.append(hist_panel(CGRA, "c_grandiflora", 44, 268, 632, 110))
    s.append(f"<text x='28' y='{H-14}' {SANS} font-size='11' fill='{DIM}'>"
             f"a_halleri: short branch, one clean bulk. c_grandiflora: deep branch "
             f"&#8212; ~100 genes pinned at &#969;&#8776;0, ~110 at the ceiling. "
             f"Count-pooled &#969; (after the dS ceiling): a_halleri 0.163, "
             f"c_grandiflora 0.151.</text>")
    s.append("</svg>")
    (OUT / "omega_estimators.svg").write_text("\n".join(s), encoding="utf-8")


def _ticks(x0, y, xs, col):
    return "".join(f"<line x1='{x0+dx}' y1='{y-7}' x2='{x0+dx}' y2='{y+7}' "
                   f"stroke='{col}' stroke-width='2'/>" for dx in xs)


def fig_one_gene():
    W, H = 720, 288
    s = [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {W} {H}'>",
         f"<rect width='{W}' height='{H}' fill='#ffffff'/>",
         f"<text x='28' y='27' {SANS} font-size='15' font-weight='700' fill='{INK}'>"
         f"One gene, one branch: &#969; is a ratio of two small counts</text>",
         f"<text x='28' y='46' {SANS} font-size='11.5' fill='{DIM}'>"
         f"Each tick = one inferred substitution on the branch.  "
         f"<tspan fill='{SYN}' font-weight='700'>&#9632;</tspan> synonymous   "
         f"<tspan fill='{NON}' font-weight='700'>&#9632;</tspan> nonsynonymous</text>"]
    rows = [("Gene A", [18, 55, 92, 118, 150], [40], "&#969; &#8776; 0.20",
             "usable"),
            ("Gene B", [70], [], "&#969; = 0 / tiny",
             "&#8594; pinned at 0"),
            ("Gene C", [], [50, 110], "&#969; = large / 0",
             "&#8594; pinned at ceiling")]
    y = 88
    for lab, syn, non, om, tag in rows:
        s.append(f"<text x='28' y='{y+4}' {SANS} font-size='12' font-weight='600' "
                 f"fill='{INK}'>{lab}</text>")
        s.append(f"<line x1='92' y1='{y}' x2='262' y2='{y}' stroke='{DIM}' "
                 f"stroke-width='1.5'/>")
        s.append(_ticks(92, y, syn, SYN))
        s.append(_ticks(92, y, non, NON))
        s.append(f"<text x='300' y='{y+4}' {MONO} font-size='11' fill='{INK}'>{om}</text>")
        s.append(f"<text x='470' y='{y+4}' {MONO} font-size='11' fill='{DIM}'>{tag}</text>")
        y += 58
    s.append(f"<text x='28' y='{H-20}' {SANS} font-size='11' fill='{DIM}'>"
             f"Expected counts per branch are single digits, so Poisson noise is "
             f"&#177;50&#8211;100%. A conserved gene often shows no nonsynonymous "
             f"change (&#969;&#8594;0); a short branch often shows no synonymous "
             f"change (&#969;&#8594;&#8734;).</text>")
    s.append("</svg>")
    (OUT / "omega_one_gene.svg").write_text("\n".join(s), encoding="utf-8")


def fig_regimes():
    W, H = 720, 292
    s = [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {W} {H}'>",
         f"<rect width='{W}' height='{H}' fill='#ffffff'/>",
         f"<text x='28' y='27' {SANS} font-size='15' font-weight='700' fill='{INK}'>"
         f"Three regimes a foreground branch can be in</text>",
         f"<text x='28' y='46' {SANS} font-size='11.5' fill='{DIM}'>"
         f"Real pilot means for the highlighted branch.</text>"]
    panels = [
        (40,  "healthy",   "a_halleri",              "t = 0.20", "dS = 0.25", 4,
         ["unsaturated, constrained by its", "sister a_thaliana &#8594; &#969; estimable"]),
        (270, "saturated",  "c_grandiflora, ok genes", "t = 0.73", "dS = 0.83", 7,
         ["deep branch, synonymous sites", "hit repeatedly &#8594; dS unreliable"]),
        (500, "collapsed",  "c_grandiflora, flagged",  "t &#8776; 0.001", "dS &#8776; 0", 1.2,
         ["3-taxon model under-determined", "&#8594; branch driven to zero"]),
    ]
    for x0, tag, name, tt, dss, lw, note in panels:
        cx = x0 + 95
        s.append(f"<text x='{cx}' y='78' {SANS} font-size='12' font-weight='700' "
                 f"fill='{INK}' text-anchor='middle'>{tag}</text>")
        s.append(f"<line x1='{cx}' y1='92' x2='{cx}' y2='114' stroke='{DIM}' stroke-width='1.5'/>")
        s.append(f"<line x1='{cx}' y1='114' x2='{cx-32}' y2='146' stroke='{DIM}' stroke-width='1.5'/>")
        s.append(f"<line x1='{cx}' y1='114' x2='{cx+42}' y2='152' stroke='{NON}' stroke-width='{lw}'/>")
        s.append(f"<text x='{cx+48}' y='164' {MONO} font-size='9' fill='{NON}'>{name}</text>")
        s.append(f"<text x='{cx-88}' y='196' {MONO} font-size='11' fill='{INK}'>{tt}</text>")
        s.append(f"<text x='{cx-88}' y='212' {MONO} font-size='11' fill='{INK}'>{dss}</text>")
        for k, line in enumerate(note):
            s.append(f"<text x='{cx-88}' y='{236+k*14}' {SANS} font-size='10' "
                     f"fill='{DIM}'>{line}</text>")
    s.append("</svg>")
    (OUT / "omega_regimes.svg").write_text("\n".join(s), encoding="utf-8")


def fig_interpretation():
    W, H = 720, 300
    s = [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {W} {H}'>",
         f"<rect width='{W}' height='{H}' fill='#ffffff'/>",
         f"<text x='28' y='27' {SANS} font-size='15' font-weight='700' fill='{INK}'>"
         f"Reading the Ne proxy: one pooled &#969; per lineage on a shared axis</text>",
         f"<text x='28' y='46' {SANS} font-size='11.5' fill='{DIM}'>"
         f"Illustrative. Lower genome-wide &#969; = selection removes slightly-"
         f"deleterious changes more efficiently = larger Ne.</text>"]

    def axis(x0, y, w, lab):
        out = [f"<line x1='{x0}' y1='{y}' x2='{x0+w}' y2='{y}' stroke='{DIM}' "
               f"stroke-width='1.5'/>"]
        for frac, t in [(0, "0.05"), (0.5, "0.20"), (1.0, "0.35")]:
            xx = x0 + frac * w
            out.append(f"<line x1='{xx}' y1='{y-4}' x2='{xx}' y2='{y+4}' "
                       f"stroke='{DIM}'/>")
            out.append(f"<text x='{xx:.0f}' y='{y+18}' {MONO} font-size='9.5' "
                       f"fill='{DIM}' text-anchor='middle'>{t}</text>")
        out.append(f"<text x='{x0-8}' y='{y-10}' {SANS} font-size='10' "
                   f"fill='{SYN}' text-anchor='end'>&#8592; larger Ne</text>")
        out.append(f"<text x='{x0+w+8}' y='{y-10}' {SANS} font-size='10' "
                   f"fill='{NON}'>smaller Ne &#8594;</text>")
        out.append(f"<text x='{x0}' y='{y+40}' {SANS} font-size='11' "
                   f"font-weight='600' fill='{INK}'>{lab}</text>")
        return "".join(out)

    def dot(x0, w, val, y, name, col):
        xx = x0 + (val - 0.05) / 0.30 * w
        return (f"<circle cx='{xx:.0f}' cy='{y}' r='5' fill='{col}'/>"
                f"<text x='{xx:.0f}' y='{y-10}' {MONO} font-size='9' fill='{INK}' "
                f"text-anchor='middle'>{name}</text>")

    x0, w = 90, 520
    s.append(axis(x0, 110, w, "pilot &#8212; 3 Brassicaceae, &#969; all &#8776; 0.16"))
    for v, n in [(0.163, "a_hal"), (0.170, "a_tha"), (0.151, "c_gra")]:
        s.append(dot(x0, w, v, 110, n, POOL))
    s.append(f"<text x='{x0+w+14}' y='114' {SANS} font-size='10' fill='{DIM}'>"
             f"CIs overlap &#8594;<tspan x='{x0+w+14}' dy='12'>similar Ne</tspan></text>")

    s.append(axis(x0, 210, w, "hypothetical &#8212; a real Ne contrast"))
    for v, n, c in [(0.11, "sp A", POOL), (0.17, "sp B", POOL), (0.30, "sp C", NON)]:
        s.append(dot(x0, w, v, 210, n, c))
    s.append(f"<text x='{x0+w+14}' y='214' {SANS} font-size='10' fill='{DIM}'>"
             f"sp C spread out &#8594;<tspan x='{x0+w+14}' dy='12'>smaller Ne in C</tspan></text>")

    s.append(f"<text x='28' y='{H-14}' {SANS} font-size='11' fill='{DIM}'>"
             f"The bootstrap CI on each dot decides whether an apparent gap is "
             f"real. In the pilot the CIs overlap; the result is “no "
             f"detectable Ne difference among the three,” not three "
             f"separate numbers to report.</text>")
    s.append("</svg>")
    (OUT / "omega_interpretation.svg").write_text("\n".join(s), encoding="utf-8")


def fig_pipeline():
    W, H = 900, 772
    ACC2 = "#2f5233"
    BG = "#f5f8f3"
    SOFT = "#dfe8db"
    s = [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {W} {H}'>",
         f"<rect width='{W}' height='{H}' fill='#ffffff'/>",
         f"<text x='40' y='34' {SANS} font-size='18' font-weight='800' "
         f"fill='{INK}'>BUSCOmega</text>",
         f"<text x='40' y='54' {SANS} font-size='12' fill='{DIM}'>"
         f"BUSCO single-copy orthologs &#8594; per-lineage dN/dS &#8594; "
         f"a per-species N&#8203;e proxy. Stdlib Python; MAFFT / PAL2NAL / "
         f"codeml do the heavy lifting.</text>"]

    cx, bw = 250, 340       # stage column
    ix = 40                 # input labels (left)
    ox = cx + bw + 20       # outputs (right)

    def box(y, h, num, title, tool, fill, num_fill):
        g = [f"<rect x='{cx}' y='{y}' width='{bw}' height='{h}' rx='9' "
             f"fill='{fill}' stroke='{SOFT}'/>"]
        if num:
            g.append(f"<circle cx='{cx+24}' cy='{y+h/2:.0f}' r='13' "
                     f"fill='{num_fill}'/>")
            g.append(f"<text x='{cx+24}' y='{y+h/2+4:.0f}' {MONO} "
                     f"font-size='12' font-weight='700' fill='#fff' "
                     f"text-anchor='middle'>{num}</text>")
        gx = cx + (46 if num else 16)
        g.append(f"<text x='{gx}' y='{y+ (h/2-3 if tool else h/2+4):.0f}' "
                 f"{SANS} font-size='12.5' font-weight='700' "
                 f"fill='{INK}'>{title}</text>")
        if tool:
            g.append(f"<text x='{gx}' y='{y+h/2+13:.0f}' {MONO} "
                     f"font-size='9.5' fill='{DIM}'>{tool}</text>")
        return "".join(g)

    def arrow(y1, y2):
        return (f"<line x1='{cx+bw/2:.0f}' y1='{y1}' x2='{cx+bw/2:.0f}' "
                f"y2='{y2}' stroke='{DIM}' stroke-width='1.5' "
                f"marker-end='url(#a)'/>")

    def side(x, y, text, anchor, col=DIM):
        return (f"<text x='{x}' y='{y}' {MONO} font-size='9.5' fill='{col}' "
                f"text-anchor='{anchor}'>{text}</text>")

    def feed(fx, fy, ty):
        return (f"<path d='M {fx} {fy} H {cx-6}' stroke='{DIM}' "
                f"stroke-width='1.2' fill='none' marker-end='url(#a)'/>")

    s.append("<defs><marker id='a' markerWidth='7' markerHeight='7' "
             "refX='6' refY='3.5' orient='auto'>"
             f"<path d='M0 0 L7 3.5 L0 7 z' fill='{DIM}'/></marker></defs>")

    rows = [
        # y, h, num, title, tool, output-label
        (74,  40, "",  "isoform-cleaned proteomes + CDS", "one protein per gene (input contract)", None),
        (132, 44, "0", "BUSCO  (prerequisite)", "prep_optional/run_busco.sh -- protein mode", "full_table.tsv  x N species"),
        (204, 44, "1", "common single-copy set", "01_common_scos.py -- strict intersection", "common_scos.tsv"),
        (276, 44, "2", "extract per-gene sequences", "02_extract_sco_seqs.py -- exact ID + suffix-strip", "prt/&lt;id&gt;.faa , cds/&lt;id&gt;.fna"),
        (348, 44, "3", "codon alignment", "03_codon_align.py -- MAFFT --auto &#8594; PAL2NAL", "codon_aln/&lt;id&gt;.pml"),
        (420, 44, "4", "codeml analysis set-up", "04_codeml_control.py -- M0 + 2-ratio templates + labelled trees", "ctl/ , trees/ , stage4_analyses.tsv"),
        (492, 44, "5", "run codeml", "05_run_codeml.py -- batched, --jobs, per-gene retry", "&lt;analysis&gt;/batch_*/batch.mlc"),
        (564, 44, "6", "parse output", "06_parse_codeml_output.py -- + QC flags", "&lt;analysis&gt;_records.tsv"),
        (636, 48, "7", "the Ne proxy", "07_ne_proxy.py -- count-pooled omega + gene bootstrap", "ne_proxy.tsv  +  plots/*.svg"),
    ]
    # tree input feeds stage 4
    for i, (y, h, num, title, tool, out) in enumerate(rows):
        fill = BG if num else "#ffffff"
        s.append(box(y, h, num, title, tool, fill, ACC2))
        if out:
            s.append(f"<text x='{ox}' y='{y+h/2+3:.0f}' {MONO} font-size='9.5' "
                     f"fill='{ACC2}'>&#8594; {out}</text>")
        if 0 < i < len(rows):
            prev = rows[i - 1]
            s.append(arrow(prev[0] + prev[1], y))

    # species tree feeds stage 4
    s.append(f"<rect x='{ix}' y='426' width='150' height='32' rx='7' "
             f"fill='#ffffff' stroke='{SOFT}'/>")
    s.append(side(ix + 10, 446, "species tree (Newick)", "start", INK))
    s.append(f"<path d='M {ix+150} 442 H {cx-6}' stroke='{DIM}' "
             f"stroke-width='1.2' fill='none' marker-end='url(#a)'/>")
    s.append(side(ix + 10, 472, "prep_optional/species_tree.py", "start"))

    s.append(f"<text x='40' y='{H-34}' {SANS} font-size='10.5' fill='{DIM}'>"
             f"Every stage writes its own <tspan {MONO}>stageN.log</tspan> "
             f"(timestamp, params, tool version, elapsed) and a manifest of "
             f"what passed and what dropped, with the reason.</text>")
    s.append(f"<text x='40' y='{H-18}' {SANS} font-size='10.5' fill='{DIM}'>"
             f"Stdlib-only Python; runs from 3 to hundreds of taxa, any "
             f"clade, any annotation source. Two-ratio branch model, one "
             f"run per focal species.</text>")
    s.append("</svg>")
    (OUT / "pipeline.svg").write_text("\n".join(s), encoding="utf-8")


fig_estimators()
fig_one_gene()
fig_regimes()
fig_interpretation()
fig_pipeline()
print("wrote", *(p.name for p in sorted(OUT.glob("*.svg"))))
