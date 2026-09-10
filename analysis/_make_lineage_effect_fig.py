#!/usr/bin/env python3
"""Two candidate plots for 'is the per-species omega different from M0?',
built from examples/pilot_stage7_out/. Stdlib only.

  lineage_effect_lrt.svg   per-gene LRT vs the chi-square(1) null  (recommended)
  lineage_effect_dist.svg  per-gene 2-ratio omega distribution + M0 & pooled lines
"""
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
S7 = HERE.parent / "examples" / "pilot_stage7_out"
SP = ["a_halleri", "a_thaliana", "c_grandiflora"]
SANS = "font-family='ui-sans-serif,-apple-system,Segoe UI,Helvetica,Arial,sans-serif'"
MONO = "font-family='ui-monospace,Menlo,Consolas,monospace'"
INK, DIM, GRN, RUS, GOLD = "#1c231d", "#5c6657", "#2f5233", "#a8432e", "#8a6d1d"

# ---- load per-gene values -------------------------------------------------
def col(path, name):
    lines = Path(path).read_text().splitlines()
    h = lines[0].split("\t"); i = h.index(name); gi = h.index("gene_id")
    out = {}
    for ln in lines[1:]:
        c = ln.split("\t")
        try:
            out[c[gi]] = float(c[i])
        except ValueError:
            pass
    return out

m0_lnl = col(S7 / "m0_records.tsv", "lnL")
m0_om = col(S7 / "m0_records.tsv", "omega")
data = {}
for s in SP:
    p = S7 / f"two_ratio.{s}_records.tsv"
    tr_lnl = col(p, "lnL")
    # foreground = the row whose tip == s
    lines = Path(p).read_text().splitlines()
    h = lines[0].split("\t"); ti, oi, gi = h.index("tip"), h.index("omega"), h.index("gene_id")
    tr_om = {}
    for ln in lines[1:]:
        c = ln.split("\t")
        if c[ti] == s:
            try:
                tr_om[c[gi]] = float(c[oi])
            except ValueError:
                pass
    genes = [g for g in tr_lnl if g in m0_lnl]
    lrt = {g: 2 * (tr_lnl[g] - m0_lnl[g]) for g in genes}
    data[s] = dict(lrt=lrt, tr_om=tr_om)

ne = {}
for ln in (S7 / "ne_proxy.tsv").read_text().splitlines()[1:]:
    c = ln.split("\t"); ne[c[0]] = c
NEH = (S7 / "ne_proxy.tsv").read_text().splitlines()[0].split("\t")
def nef(sp, k):
    return float(ne[sp][NEH.index(k)])


def chi2_1_pdf(x):
    if x <= 0:
        return 0.0
    return math.exp(-x / 2) / math.sqrt(2 * math.pi * x)


# ---- plot 1: LRT vs chi2(1) --------------------------------------------
def fig_lrt():
    W, ph, pad = 760, 150, 60
    H = pad + ph * 3 + 40
    x0, x1 = 60, W - 24
    XMAX = 12
    NB = 30
    s = [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {W} {H}'>",
         f"<rect width='{W}' height='{H}' fill='#ffffff'/>",
         f"<text x='24' y='26' {SANS} font-size='14' font-weight='700' fill='{INK}'>"
         f"Lineage effect: per-gene LRT (2-ratio vs M0) against the &#967;&#178;(1) null</text>",
         f"<text x='24' y='44' {SANS} font-size='11' fill='{DIM}'>"
         f"If species has no distinct &#969;, the bars follow the curve and "
         f"~5% of genes exceed 3.84.</text>"]
    for k, sp in enumerate(SP):
        top = pad + k * ph
        base = top + ph - 34
        lrt = [max(0, v) for v in data[sp]["lrt"].values()]
        n = len(lrt)
        bw = (x1 - x0) / NB
        bins = [0] * NB
        for v in lrt:
            b = min(NB - 1, int(v / XMAX * NB))
            bins[b] += 1
        mx = max(bins) or 1
        # expected chi2(1) density scaled to counts
        s.append(f"<text x='{x0}' y='{top+2}' {MONO} font-size='10.5' "
                 f"font-weight='700' fill='{INK}'>{sp}</text>")
        for i, c in enumerate(bins):
            hh = c / mx * (ph - 46)
            s.append(f"<rect x='{x0+i*bw:.1f}' y='{base-hh:.1f}' "
                     f"width='{bw-1:.1f}' height='{hh:.1f}' fill='#cbd6cb'/>")
        # chi2 curve
        pts = []
        for j in range(1, 200):
            xx = j / 200 * XMAX
            dens = chi2_1_pdf(xx) * n * (XMAX / NB)
            yy = base - dens / mx * (ph - 46)
            pts.append(f"{x0 + xx/XMAX*(x1-x0):.1f},{max(top+8, yy):.1f}")
        s.append(f"<polyline points='{' '.join(pts)}' fill='none' "
                 f"stroke='{GRN}' stroke-width='1.8'/>")
        # 3.84 line
        xt = x0 + 3.841 / XMAX * (x1 - x0)
        s.append(f"<line x1='{xt:.1f}' y1='{top+8}' x2='{xt:.1f}' y2='{base}' "
                 f"stroke='{RUS}' stroke-width='1.4' stroke-dasharray='3 2'/>")
        s.append(f"<text x='{xt+4:.0f}' y='{top+18}' {MONO} font-size='8.5' "
                 f"fill='{RUS}'>p=0.05</text>")
        s.append(f"<line x1='{x0}' y1='{base}' x2='{x1}' y2='{base}' stroke='{DIM}'/>")
        for gx in (0, 3.841, 6, 9, 12):
            xx = x0 + gx / XMAX * (x1 - x0)
            s.append(f"<text x='{xx:.0f}' y='{base+14}' {MONO} font-size='8.5' "
                     f"fill='{DIM}' text-anchor='middle'>{gx:g}</text>")
        nsig = int(nef(sp, "n_lrt_p05")); ntot = int(nef(sp, "n_lrt"))
        ml = nef(sp, "mean_lrt")
        s.append(f"<text x='{x1}' y='{top+2}' {MONO} font-size='9' fill='{INK}' "
                 f"text-anchor='end'>{nsig}/{ntot} genes p&lt;0.05 "
                 f"({nsig/ntot:.1%}, exp 5%) &#183; mean LRT {ml:.2f} (exp 1.0)</text>")
    s.append(f"<text x='{(x0+x1)/2:.0f}' y='{H-8}' {SANS} font-size='10' "
             f"fill='{DIM}' text-anchor='middle'>2 &#215; (lnL 2-ratio &#8722; lnL M0)</text>")
    s.append("</svg>")
    (HERE / "lineage_effect_lrt.svg").write_text("\n".join(s), encoding="utf-8")


# ---- plot 2: per-gene omega dist + M0 & pooled lines -------------------
def fig_dist():
    W, ph, pad = 760, 150, 56
    H = pad + ph * 3 + 34
    x0, x1 = 60, W - 24
    NB = 40
    s = [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {W} {H}'>",
         f"<rect width='{W}' height='{H}' fill='#ffffff'/>",
         f"<text x='24' y='26' {SANS} font-size='14' font-weight='700' fill='{INK}'>"
         f"Per-gene 2-ratio &#969; per species, with M0 and the pooled estimate</text>",
         f"<text x='24' y='44' {SANS} font-size='11' fill='{DIM}'>"
         f"Bars: per-gene foreground &#969; (&#8805;1 lumped at the right). "
         f"<tspan fill='{GOLD}'>M0</tspan>, <tspan fill='{RUS}'>pooled 2-ratio "
         f"(95% CI band)</tspan>.</text>"]
    for k, sp in enumerate(SP):
        top = pad + k * ph
        base = top + ph - 30
        oms = [v for v in data[sp]["tr_om"].values()]
        n = len(oms)
        bw = (x1 - x0) / (NB + 1)
        bins = [0] * (NB + 1)
        for o in oms:
            bins[NB if o >= 1 else max(0, int(o * NB))] += 1
        mx = max(bins) or 1
        s.append(f"<text x='{x0}' y='{top+2}' {MONO} font-size='10.5' "
                 f"font-weight='700' fill='{INK}'>{sp}  (n={n})</text>")
        for i, c in enumerate(bins):
            hh = c / mx * (ph - 44)
            fill = "#cbd6cb" if i < NB else "#e3b7a8"
            s.append(f"<rect x='{x0+i*bw:.1f}' y='{base-hh:.1f}' "
                     f"width='{bw-1:.1f}' height='{hh:.1f}' fill='{fill}'/>")
        span = x1 - x0 - bw
        pooled, lo, hi = (nef(sp, "omega_pooled"), nef(sp, "ci_lo"),
                          nef(sp, "ci_hi"))
        m0v = nef(sp, "omega_M0")
        s.append(f"<rect x='{x0+lo*span:.1f}' y='{top+8}' "
                 f"width='{(hi-lo)*span:.1f}' height='{base-top-8:.1f}' "
                 f"fill='{RUS}' fill-opacity='0.10'/>")
        for v, col_, lab in [(m0v, GOLD, "M0"), (pooled, RUS, "pooled")]:
            xx = x0 + min(v, 1.0) * span
            s.append(f"<line x1='{xx:.1f}' y1='{top+8}' x2='{xx:.1f}' "
                     f"y2='{base}' stroke='{col_}' stroke-width='1.6' "
                     f"stroke-dasharray='4 2'/>")
            s.append(f"<text x='{xx+3:.1f}' y='{top+18 if lab=='M0' else top+30}' "
                     f"{MONO} font-size='8.5' fill='{col_}'>{lab} {v:.3f}</text>")
        s.append(f"<line x1='{x0}' y1='{base}' x2='{x1}' y2='{base}' stroke='{DIM}'/>")
        for frac, lab in [(0, "0"), (0.25, "0.25"), (0.5, "0.5"),
                          (0.75, "0.75")]:
            s.append(f"<text x='{x0+frac*span:.0f}' y='{base+13}' {MONO} "
                     f"font-size='8.5' fill='{DIM}' text-anchor='middle'>{lab}</text>")
        s.append(f"<text x='{x0+NB*bw+bw/2:.0f}' y='{base+13}' {MONO} "
                 f"font-size='8.5' fill='{RUS}' text-anchor='middle'>&#8805;1</text>")
    s.append("</svg>")
    (HERE / "lineage_effect_dist.svg").write_text("\n".join(s), encoding="utf-8")


fig_lrt()
fig_dist()
print("wrote lineage_effect_lrt.svg, lineage_effect_dist.svg")
