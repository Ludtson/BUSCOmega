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
    """Hero illustration: many genomes narrow to a conserved working set,
    pass through the tree, and fan back out to one number per species."""
    W, H = 1200, 430
    BG = "#FBFAF6"
    G1, G2, GS = "#2F5233", "#5B8A63", "#DDE8DC"
    RU, RUS = "#A8432E", "#ECD5CD"
    GD, GY = "#8A6D1D", "#A8B0A4"
    cxs = [128, 278, 428, 578, 728, 878, 1028]
    yc = 214                                        # glyph band centre
    R = 52                                          # glyph half-size

    s = [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {W} {H}'>",
         "<defs>",
         f"<linearGradient id='rib' x1='0' y1='0' x2='1' y2='0'>"
         f"<stop offset='0' stop-color='{GS}'/><stop offset='0.5' "
         f"stop-color='{GS}'/><stop offset='1' stop-color='{RUS}'/>"
         f"</linearGradient>",
         "<marker id='ah' markerWidth='8' markerHeight='8' refX='6' refY='4' "
         f"orient='auto'><path d='M0 0 L8 4 L0 8 z' fill='{G1}'/></marker>",
         "</defs>",
         f"<rect width='{W}' height='{H}' fill='{BG}'/>",
         f"<text x='44' y='42' {SANS} font-size='21' font-weight='800' "
         f"fill='{INK}' letter-spacing='-0.5'>BUSCOmega</text>",
         f"<text x='44' y='64' {SANS} font-size='12.5' fill='{DIM}'>"
         f"Many genomes narrow to the orthologs they share, pass through the "
         f"tree one lineage at a time, and fan back out to one N&#8203;e "
         f"proxy per species.</text>"]

    # the ribbon: wide (all genes) -> pinched (working set) -> flared (per species)
    top, bot = yc - 66, yc + 66
    p3, p5 = cxs[2], cxs[4]
    s.append(
        f"<path d='M {cxs[0]-70} {top} "
        f"C {p3-90} {top}, {p3-30} {yc-26}, {p3} {yc-26} "
        f"C {p5-30} {yc-26}, {p5-60} {top-4}, {cxs[6]+70} {top-14} "
        f"L {cxs[6]+70} {bot+14} "
        f"C {p5-60} {bot+4}, {p5-30} {yc+26}, {p5} {yc+26} "
        f"C {p3-30} {yc+26}, {p3-90} {bot}, {cxs[0]-70} {bot} Z' "
        f"fill='url(#rib)' opacity='0.8'/>")

    # ---- glyphs -----------------------------------------------------------
    def circle_venn(cx):
        o = []
        for dx in (-15, 15, 0):
            dy = 0 if dx else 16
            o.append(f"<circle cx='{cx+dx}' cy='{yc-6+dy}' r='24' "
                     f"fill='{G2}' fill-opacity='0.28' stroke='{G2}'/>")
        o.append(f"<circle cx='{cx}' cy='{yc+2}' r='7' fill='{G1}'/>")
        return "".join(o)

    def seq_stack(cx):
        o = []
        for i, col in enumerate((G1, RU, G2)):
            y = yc - 24 + i * 17
            o.append(f"<path d='M {cx-30} {y} h 8' stroke='{col}' "
                     f"stroke-width='3'/>")
            o.append(f"<path d='M {cx-18} {y} q 8 -6 16 0 t 16 0 t 12 0' "
                     f"stroke='{col}' stroke-width='2' fill='none'/>")
        return "".join(o)

    def codon_grid(cx):
        o = []
        pat = ["GGGrGGGG", "GGrGG-GG", "GGGGGGrG", "GrGG-GGG"]
        for r, row in enumerate(pat):
            for c, ch in enumerate(row):
                col = {"G": G2, "r": RU, "-": GY}[ch]
                op = "0.35" if ch == "-" else "0.9"
                o.append(f"<rect x='{cx-32+c*8.2:.1f}' y='{yc-22+r*11:.1f}' "
                         f"width='6.6' height='9' rx='1.5' fill='{col}' "
                         f"fill-opacity='{op}'/>")
        return "".join(o)

    def tree(cx, fg=None):
        # 4-tip tree, root left; fg = index 0..3 of the highlighted tip
        tips = [yc - 27, yc - 9, yc + 9, yc + 27]
        o = [f"<path d='M {cx-34} {yc} h 14' stroke='{G1}' stroke-width='2' "
             f"fill='none'/>"]
        o.append(f"<path d='M {cx-20} {yc} V {tips[0]} V {tips[3]}' "
                 f"stroke='{G1}' stroke-width='2' fill='none'/>")
        o.append(f"<path d='M {cx-20} {(tips[0]+tips[1])/2:.0f} h 10 "
                 f"V {tips[0]} M {cx-10} {(tips[0]+tips[1])/2:.0f} V {tips[1]}' "
                 f"stroke='{G1}' stroke-width='2' fill='none'/>")
        for i, ty in enumerate(tips):
            x0 = cx - 20 if i in (0, 3) else cx - 10
            hot = (i == fg)
            o.append(f"<path d='M {x0} {ty} H {cx+28}' stroke="
                     f"'{RU if hot else G1}' stroke-width='{4 if hot else 2}' "
                     f"fill='none'/>")
            o.append(f"<circle cx='{cx+28}' cy='{ty}' r='{3.5 if hot else 2.5}' "
                     f"fill='{RU if hot else G2}'/>")
        if fg is not None:
            o.append(f"<text x='{cx+34}' y='{tips[fg]+3}' {MONO} "
                     f"font-size='9' font-weight='700' fill='{RU}'>#1</text>")
        return "".join(o)

    def doc_tree(cx):
        o = [f"<rect x='{cx-34}' y='{yc-30}' width='30' height='40' rx='3' "
             f"fill='#fff' stroke='{G2}'/>"]
        for k in range(4):
            o.append(f"<path d='M {cx-28} {yc-22+k*8} h 18' stroke='{GY}' "
                     f"stroke-width='1.5'/>")
        o.append(tree(cx + 16))
        return "".join(o)

    def table(cx):
        o = []
        for r in range(3):
            for c in range(3):
                o.append(f"<rect x='{cx-30+c*20}' y='{yc-22+r*15}' width='17' "
                         f"height='12' rx='1.5' fill='#fff' stroke='{GY}'/>")
        o.append(f"<rect x='{cx-30}' y='{yc-22}' width='57' height='12' "
                 f"fill='{G2}' fill-opacity='0.25'/>")
        return "".join(o)

    def forest(cx):
        o = [f"<line x1='{cx-34}' y1='{yc+26}' x2='{cx+34}' y2='{yc+26}' "
             f"stroke='{GY}'/>"]
        for i, (mx, w) in enumerate([(-6, 12), (4, 10), (-14, 11)]):
            y = yc - 16 + i * 16
            o.append(f"<line x1='{cx+mx-w}' y1='{y}' x2='{cx+mx+w}' y2='{y}' "
                     f"stroke='{RU}' stroke-width='2'/>")
            o.append(f"<circle cx='{cx+mx}' cy='{y}' r='3.5' fill='{RU}'/>")
        return "".join(o)

    glyphs = [circle_venn, seq_stack, codon_grid, doc_tree,
              lambda cx: tree(cx, fg=1) +
              f"<text x='{cx-2}' y='{yc+34}' {MONO} font-size='10' "
              f"fill='{G1}' text-anchor='middle'>&#969;</text>",
              table, forest]
    caps = [
        ("1", "common set", "01_common_scos.py"),
        ("2", "sequences", "02_extract_sco_seqs.py"),
        ("3", "codon align", "03_codon_align.py"),
        ("4", "codeml set-up", "04_codeml_control.py"),
        ("5", "run codeml", "05_run_codeml.py"),
        ("6", "parse", "06_parse_codeml_output.py"),
        ("7", "Ne proxy", "07_ne_proxy.py"),
    ]
    for i, cx in enumerate(cxs):
        s.append(f"<circle cx='{cx}' cy='{yc}' r='{R}' fill='#fff' "
                 f"stroke='{GS}' stroke-width='2'/>")
        s.append(glyphs[i](cx))
        if i:
            s.append(f"<path d='M {cxs[i-1]+R+3} {yc} H {cx-R-6}' "
                     f"stroke='{G1}' stroke-width='1.6' fill='none' "
                     f"marker-end='url(#ah)' opacity='0.8'/>")
        num, name, script = caps[i]
        by = yc + R + 34
        s.append(f"<circle cx='{cx}' cy='{by-6}' r='11' fill='{G1}'/>")
        s.append(f"<text x='{cx}' y='{by-2}' {MONO} font-size='11' "
                 f"font-weight='700' fill='#fff' text-anchor='middle'>{num}</text>")
        s.append(f"<text x='{cx}' y='{by+20}' {SANS} font-size='11.5' "
                 f"font-weight='700' fill='{INK}' text-anchor='middle'>{name}</text>")
        s.append(f"<text x='{cx}' y='{by+35}' {MONO} font-size='8.5' "
                 f"fill='{DIM}' text-anchor='middle'>{script}</text>")

    # left inflow: genomes + BUSCO
    gx = cxs[0] - 88
    for k in range(3):
        for j in range(3):
            s.append(f"<rect x='{gx-9+j*9}' y='{yc-30+k*10}' width='7' "
                     f"height='5' rx='1' fill='{G2}' fill-opacity='0.85'/>")
    s.append(f"<text x='{gx+4}' y='{yc+31}' {SANS} font-size='9.5' fill='{DIM}' "
             f"text-anchor='middle'>genomes</text>")
    s.append(f"<text x='{gx+4}' y='{yc+43}' {SANS} font-size='9.5' fill='{DIM}' "
             f"text-anchor='middle'>+ BUSCO</text>")
    s.append(f"<path d='M {gx+22} {yc} H {cxs[0]-R-4}' stroke='{G1}' "
             f"stroke-width='1.6' marker-end='url(#ah)' opacity='0.8'/>")

    # tree feed into stage 4
    s.append(f"<path d='M {cxs[3]} {yc-R-30} V {yc-R-4}' stroke='{G1}' "
             f"stroke-width='1.5' fill='none' marker-end='url(#ah)' "
             f"opacity='0.7'/>")
    s.append(f"<text x='{cxs[3]}' y='{yc-R-38}' {SANS} font-size='9.5' "
             f"fill='{DIM}' text-anchor='middle'>species tree</text>")

    s.append(f"<text x='44' y='{H-46}' {SANS} font-size='11' fill='{DIM}'>"
             f"<tspan font-weight='700' fill='{INK}'>N&#8203;e proxy</tspan> = "
             f"count-pooled dN/dS on each species' own branch "
             f"(two-ratio branch model), with a gene-bootstrap CI.</text>")
    s.append(f"<text x='44' y='{H-28}' {SANS} font-size='11' fill='{DIM}'>"
             f"Standard-library Python; MAFFT / PAL2NAL / codeml do the "
             f"alignment and model fitting. Every stage: a timestamped log "
             f"and a pass/drop manifest.</text>")
    s.append(f"<text x='44' y='{H-10}' {MONO} font-size='9' fill='{GY}'>"
             f"run_buscomega.py orchestrates all seven &#183; runs 3 to "
             f"hundreds of taxa, any clade</text>")
    s.append("</svg>")
    (OUT / "pipeline.svg").write_text("\n".join(s), encoding="utf-8")


fig_estimators()
fig_one_gene()
fig_regimes()
fig_interpretation()
fig_pipeline()
print("wrote", *(p.name for p in sorted(OUT.glob("*.svg"))))
