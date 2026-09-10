#!/usr/bin/env python3
"""Dot plot of pilot_method_comparison.tsv. Stdlib only."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
SANS = "font-family='ui-sans-serif,-apple-system,Segoe UI,Helvetica,Arial,sans-serif'"
MONO = "font-family='ui-monospace,Menlo,Consolas,monospace'"
INK, DIM = "#1c231d", "#5c6657"

# (label, colour, marker)  -- order matters for the legend
METHODS = [
    ("M0 per-species (dS-filtered)",            "#8A6D1D", "d"),
    ("2-ratio per-gene pooled (BUSCOmega)",     "#2F5233", "o"),
    ("free-ratio per-gene pooled",              "#5B8A63", "s"),
    ("free-ratio concat, misalignment-filtered (approx eLife)", "#A8432E", "^"),
]
SPECIES = ["a_halleri", "a_thaliana", "c_grandiflora"]

rows, ngenes = {}, {}
for ln in (HERE / "pilot_method_comparison.tsv").read_text().splitlines()[1:]:
    c = ln.split("\t")
    rows[c[0]] = {sp: float(v) for sp, v in zip(SPECIES, c[1:4])}
    ngenes[c[0]] = {sp: int(v) for sp, v in zip(SPECIES, c[4:7])}

NG = {
    "M0 per-species (dS-filtered)": ngenes["M0 per-species (dS-filtered)"],
    "2-ratio per-gene pooled (BUSCOmega)":
        ngenes["2-ratio per-gene pooled -- BUSCOmega"],
    "free-ratio per-gene pooled": ngenes["free-ratio per-gene pooled"],
    "free-ratio concat, misalignment-filtered (approx eLife)":
        ngenes["free-ratio concatenate -- misalignment-filtered (approx eLife)"],
}

VALS = {
    "M0 per-species (dS-filtered)": rows["M0 per-species (dS-filtered)"],
    "2-ratio per-gene pooled (BUSCOmega)":
        rows["2-ratio per-gene pooled -- BUSCOmega"],
    "free-ratio per-gene pooled": rows["free-ratio per-gene pooled"],
    "free-ratio concat, misalignment-filtered (approx eLife)":
        rows["free-ratio concatenate -- misalignment-filtered (approx eLife)"],
}

W, H = 720, 340
x0, x1 = 150, W - 40
xmin, xmax = 0.12, 0.22
row_y = {sp: 70 + i * 74 for i, sp in enumerate(SPECIES)}


def px(v):
    return x0 + (v - xmin) / (xmax - xmin) * (x1 - x0)


s = [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {W} {H}'>",
     f"<rect width='{W}' height='{H}' fill='#ffffff'/>",
     f"<text x='28' y='26' {SANS} font-size='14' font-weight='700' fill='{INK}'>"
     f"Per-species &#969; on the 3-taxon pilot: four estimation routes</text>",
     f"<text x='28' y='44' {SANS} font-size='11' fill='{DIM}'>"
     f"369 BUSCO genes. The branch model barely moves the estimate; "
     f"gene filtering does.</text>"]

for frac in (0, 0.25, 0.5, 0.75, 1.0):
    xv = xmin + frac * (xmax - xmin)
    s.append(f"<line x1='{px(xv):.0f}' y1='58' x2='{px(xv):.0f}' y2='{H-46}' "
             f"stroke='#eee'/>")
    s.append(f"<text x='{px(xv):.0f}' y='{H-30}' {MONO} font-size='9' "
             f"fill='{DIM}' text-anchor='middle'>{xv:.2f}</text>")
s.append(f"<text x='{(x0+x1)/2:.0f}' y='{H-12}' {SANS} font-size='10.5' "
         f"fill='{DIM}' text-anchor='middle'>count-pooled &#969; (= dN/dS)</text>")

for sp in SPECIES:
    y = row_y[sp]
    s.append(f"<text x='{x0-12}' y='{y+4}' {MONO} font-size='11' fill='{INK}' "
             f"text-anchor='end'>{sp}</text>")
    s.append(f"<line x1='{x0}' y1='{y}' x2='{x1}' y2='{y}' stroke='#f0f0f0'/>")
    for j, (label, col, mk) in enumerate(METHODS):
        v = VALS[label][sp]
        cx, cy = px(v), y - 18 + j * 12
        if mk == "o":
            s.append(f"<circle cx='{cx:.1f}' cy='{cy}' r='4.5' fill='{col}'/>")
        elif mk == "s":
            s.append(f"<rect x='{cx-4:.1f}' y='{cy-4}' width='8' height='8' "
                     f"fill='{col}'/>")
        elif mk == "d":
            s.append(f"<path d='M {cx:.1f} {cy-5} l 5 5 l -5 5 l -5 -5 z' "
                     f"fill='{col}'/>")
        else:
            s.append(f"<path d='M {cx:.1f} {cy-5} l 5 9 l -10 0 z' fill='{col}'/>")
        s.append(f"<text x='{cx+8:.1f}' y='{cy+3}' {MONO} font-size='8' "
                 f"fill='{col}'>{v:.3f}  <tspan fill='{DIM}'>n={NG[label][sp]}</tspan></text>")

# legend
ly = H - 46
for j, (label, col, mk) in enumerate(METHODS):
    lx = 28 + j * 168
    s.append(f"<rect x='{lx}' y='{ly-7}' width='9' height='9' fill='{col}'/>")
    s.append(f"<text x='{lx+13}' y='{ly+1}' {SANS} font-size='8.5' "
             f"fill='{DIM}'>{label}</text>")

s.append("</svg>")
(HERE / "method_comparison.svg").write_text("\n".join(s), encoding="utf-8")
print("wrote method_comparison.svg")
