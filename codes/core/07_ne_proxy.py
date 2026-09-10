#!/usr/bin/env python3
"""
07_ne_proxy.py

Stage 7 of BUSCOmega: turn the per-gene, per-branch dN/dS records from
Stage 6 into one Ne proxy per species.

The proxy for species X is the **count-pooled omega** on X's own terminal
branch, taken from X's two-ratio run. It is the concatenation-equivalent
estimator: pool synonymous and nonsynonymous *sites* and *substitutions*
across genes, then take the ratio once.

    dN_pool = sum_g( N_g * dN_{g,X} ) / sum_g( N_g )      # subs / sites
    dS_pool = sum_g( S_g * dS_{g,X} ) / sum_g( S_g )
    omega_X = dN_pool / dS_pool

This is *not* the mean or median of the per-gene ratios, which are
boundary-dominated and uninterpretable (see docs/primer.md section 8).
Mean and median are written alongside only to show the difference.

Gene filtering before pooling (docs/primer.md section 8.5):
  - keep genes flagged `omega_boundary` (they add ~0 to both sums)
  - drop genes flagged `ds_floor` (no synonymous signal on that branch ->
    unbalanced contribution)
  - drop a gene whose dS on that branch exceeds `--ds-ceiling` (default 1.5;
    synonymous saturation -> unreliable denominator)

Uncertainty: a gene bootstrap (`--bootstrap`, default 1000) -- resample
genes with replacement, recompute pooled omega, report the 2.5/97.5 percentiles.

Input (one of):
    --records-dir DIR   directory of Stage 6 tables: m0_records.tsv and
                        one two_ratio.<species>_records.tsv per focal species
    --stage5-dir DIR    a Stage 5 out-dir; Stage 6 is run on it first
                        (needs 06_parse_codeml_output.py alongside this file)

Outputs (into --out-dir, default 07_ne_proxy/):
    ne_proxy.tsv        species  omega_pooled  ci_lo  ci_hi  n_genes
                        n_excl_ds_floor  n_excl_ds_ceiling  mean_omega
                        median_omega  mean_t  omega_M0
    per_gene_omega.tsv  species  gene_id  N  S  dN  dS  omega  t  qc_flag
                        used  excl_reason
    plots/*.svg         forest plot, per-gene distributions, omega vs divergence
    plots/*.png         only with --png (needs rsvg-convert / inkscape / cairosvg)
    stage7.log

Stdlib only.

Author: Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)
Version: 0.1.0
"""
from __future__ import annotations

__version__ = "0.1.0"
__author__ = "Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)"

import argparse
import random
import shutil
import statistics as st
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

DS_CEILING_DEFAULT = 1.5
BOOTSTRAP_DEFAULT = 1000


def stage_log(logpath: Path, stage: str, event: str, **fields) -> None:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")
    kv = "  ".join(f"{k}={v}" for k, v in fields.items())
    with logpath.open("a", encoding="utf-8") as fh:
        fh.write(f"[{ts}] {stage} {event}  {kv}\n")


# --------------------------------------------------------------------------
# records
# --------------------------------------------------------------------------
class Row:
    __slots__ = ("gene_id", "tip", "t", "N", "S", "omega", "dN", "dS", "qc",
                 "lnL")

    def __init__(self, d: dict):
        self.gene_id = d["gene_id"]
        self.tip = d.get("tip") or ""
        self.t = float(d["t"])
        self.N = float(d["N"])
        self.S = float(d["S"])
        self.omega = float(d["omega"])
        self.dN = float(d["dN"])
        self.dS = float(d["dS"])
        self.qc = d.get("qc_flag", "ok")
        v = d.get("lnL", "")
        self.lnL = float(v) if v not in ("", "None") else None


def read_records(path: Path) -> list[Row]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    out = []
    for ln in lines[1:]:
        if ln.strip():
            out.append(Row(dict(zip(header, ln.split("\t")))))
    return out


# --------------------------------------------------------------------------
# the estimator
# --------------------------------------------------------------------------
def pooled_omega(rows: list[Row]) -> float:
    """Concatenation-equivalent omega: pool sites and substitutions, divide once.

        dN_pool = sum(N*dN)/sum(N) ; dS_pool = sum(S*dS)/sum(S) ; omega = dN_pool/dS_pool
    """
    sN = sum(r.N for r in rows)
    sS = sum(r.S for r in rows)
    if sN <= 0 or sS <= 0:
        return float("nan")
    dn = sum(r.N * r.dN for r in rows) / sN
    ds = sum(r.S * r.dS for r in rows) / sS
    return dn / ds if ds > 0 else float("nan")


def bootstrap_ci(rows: list[Row], b: int, seed: int) -> tuple[float, float]:
    """Resample genes with replacement, recompute pooled omega, 2.5/97.5 pct."""
    by_gene: dict[str, list[Row]] = {}
    for r in rows:
        by_gene.setdefault(r.gene_id, []).append(r)
    genes = list(by_gene)
    if len(genes) < 2 or b < 2:
        return (float("nan"), float("nan"))
    rng = random.Random(seed)
    vals = []
    for _ in range(b):
        nsub = ssub = nsite = ssite = 0.0
        for _ in range(len(genes)):
            for r in by_gene[genes[rng.randrange(len(genes))]]:
                nsub += r.N * r.dN
                ssub += r.S * r.dS
                nsite += r.N
                ssite += r.S
        if nsite > 0 and ssite > 0 and ssub > 0:
            vals.append((nsub / nsite) / (ssub / ssite))
    if len(vals) < 2:
        return (float("nan"), float("nan"))
    vals.sort()
    lo = vals[int(0.025 * (len(vals) - 1))]
    hi = vals[int(round(0.975 * (len(vals) - 1)))]
    return (lo, hi)


def filter_focal(rows: list[Row], species: str, ds_ceiling: float):
    """Return (used_rows, per_gene_annotations, n_ds_floor, n_ceiling)."""
    used, annot = [], []
    n_floor = n_ceil = 0
    for r in rows:
        if r.tip != species:
            continue
        reason = ""
        if r.qc == "ds_floor":
            reason = "ds_floor"
            n_floor += 1
        elif r.dS > ds_ceiling:
            reason = f"ds>{ds_ceiling}"
            n_ceil += 1
        if reason:
            annot.append((r, False, reason))
        else:
            used.append(r)
            annot.append((r, True, ""))
    return used, annot, n_floor, n_ceil


# --------------------------------------------------------------------------
# minimal SVG plotting
# --------------------------------------------------------------------------
INK, DIM, ACC = "#1c231d", "#5c6657", "#31607d"
SANS = "font-family='ui-sans-serif,-apple-system,Segoe UI,Helvetica,Arial,sans-serif'"
MONO = "font-family='ui-monospace,Menlo,Consolas,monospace'"


def _svg(w, h, body: str) -> str:
    return (f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {w} {h}'>"
            f"<rect width='{w}' height='{h}' fill='#ffffff'/>{body}</svg>")


def forest_plot(rows_by_sp: dict, m0: float, path: Path) -> None:
    sp = list(rows_by_sp)
    vals = [rows_by_sp[s] for s in sp]           # (omega, lo, hi)
    xs = [v for tup in vals for v in tup if v == v]        # drop NaN
    if m0 == m0:
        xs.append(m0)
    if not xs:
        path.write_text(_svg(400, 80,
            f"<text x='20' y='45' {SANS} font-size='13' fill='{DIM}'>"
            f"no finite pooled omega to plot</text>"), encoding="utf-8")
        return
    xmin, xmax = min(xs) * 0.9, max(xs) * 1.1
    if xmax - xmin < 1e-9:
        xmin, xmax = xmin * 0.5, xmax * 1.5 + 1e-6
    W, H = 620, 60 + 34 * len(sp)
    x0, x1 = 150, W - 30

    def px(v):
        return x0 + (v - xmin) / (xmax - xmin) * (x1 - x0)

    b = [f"<text x='24' y='28' {SANS} font-size='14' font-weight='700' "
         f"fill='{INK}'>Ne proxy per species (count-pooled &#969;, 95% CI)</text>"]
    if m0 == m0:
        b.append(f"<line x1='{px(m0):.1f}' y1='42' x2='{px(m0):.1f}' y2='{H-16}' "
                 f"stroke='{DIM}' stroke-dasharray='3 3'/>")
        b.append(f"<text x='{px(m0):.1f}' y='{H-4}' {MONO} font-size='9' "
                 f"fill='{DIM}' text-anchor='middle'>M0 {m0:.3f}</text>")
    y = 58
    for s in sp:
        om, lo, hi = rows_by_sp[s]
        b.append(f"<text x='{x0-10}' y='{y+4}' {MONO} font-size='11' fill='{INK}' "
                 f"text-anchor='end'>{s}</text>")
        if lo == lo and hi == hi:
            b.append(f"<line x1='{px(lo):.1f}' y1='{y}' x2='{px(hi):.1f}' y2='{y}' "
                     f"stroke='{ACC}' stroke-width='2'/>")
        b.append(f"<circle cx='{px(om):.1f}' cy='{y}' r='4' fill='{ACC}'/>")
        b.append(f"<text x='{px(om):.1f}' y='{y-9}' {MONO} font-size='9' "
                 f"fill='{INK}' text-anchor='middle'>{om:.3f}</text>")
        y += 34
    for frac in (0, 0.5, 1.0):
        xv = xmin + frac * (xmax - xmin)
        b.append(f"<text x='{px(xv):.0f}' y='{H-16}' {MONO} font-size='9' "
                 f"fill='{DIM}' text-anchor='middle'>{xv:.3f}</text>")
    path.write_text(_svg(W, H, "".join(b)), encoding="utf-8")


def dist_plot(annot_by_sp: dict, pooled_by_sp: dict, path: Path) -> None:
    sp = list(annot_by_sp)
    W, ph = 620, 96
    H = 50 + ph * len(sp)
    x0, x1 = 44, W - 24
    NB = 20
    b = [f"<text x='24' y='26' {SANS} font-size='14' font-weight='700' "
         f"fill='{INK}'>Per-gene &#969; by species (dashed = count-pooled)</text>"]
    y0 = 44
    for s in sp:
        oms = [a[0].omega for a in annot_by_sp[s] if a[0].omega == a[0].omega]
        bins = [0] * (NB + 1)
        for o in oms:
            bins[NB if o >= 1 else max(0, int(o * NB))] += 1
        mx = max(bins) or 1
        bw = (x1 - x0) / (NB + 1)
        b.append(f"<text x='{x0}' y='{y0-4}' {MONO} font-size='10' fill='{INK}'>"
                 f"{s}  (n={len(oms)})</text>")
        for i, c in enumerate(bins):
            hh = c / mx * (ph - 26)
            fill = "#cbd6cb" if i < NB else "#e3b7a8"
            b.append(f"<rect x='{x0+i*bw:.1f}' y='{y0+ph-26-hh:.1f}' "
                     f"width='{bw-1:.1f}' height='{hh:.1f}' fill='{fill}'/>")
        pv = pooled_by_sp[s]
        if pv == pv:
            xx = x0 + min(pv, 1.0) * (x1 - x0 - bw)
            b.append(f"<line x1='{xx:.1f}' y1='{y0}' x2='{xx:.1f}' "
                     f"y2='{y0+ph-26}' stroke='{ACC}' stroke-width='1.4' "
                     f"stroke-dasharray='3 2'/>")
        b.append(f"<line x1='{x0}' y1='{y0+ph-26}' x2='{x1}' y2='{y0+ph-26}' "
                 f"stroke='{DIM}'/>")
        for frac, lab in [(0, "0"), (0.5, ".5"), (1.0, "&#8805;1")]:
            b.append(f"<text x='{x0+frac*(x1-x0-bw):.0f}' y='{y0+ph-12}' {MONO} "
                     f"font-size='8.5' fill='{DIM}' text-anchor='middle'>{lab}</text>")
        y0 += ph
    path.write_text(_svg(W, H, "".join(b)), encoding="utf-8")


def omega_vs_t_plot(points: list, path: Path) -> None:
    # points: (species, mean_t, omega_pooled)
    W, H = 560, 360
    x0, x1, y0, y1 = 70, W - 24, 40, H - 44
    ts = [p[1] for p in points if p[1] == p[1]]
    os_ = [p[2] for p in points if p[2] == p[2]]
    if not ts or not os_:
        path.write_text(_svg(400, 80,
            f"<text x='20' y='45' {SANS} font-size='13' fill='{DIM}'>"
            f"not enough finite points to plot</text>"), encoding="utf-8")
        return
    tmin, tmax = min(ts) * 0.9, max(ts) * 1.1
    omin, omax = min(os_) * 0.9, max(os_) * 1.1
    if tmax - tmin < 1e-9:
        tmin, tmax = tmin * 0.5, tmax * 1.5 + 1e-6
    if omax - omin < 1e-9:
        omin, omax = omin * 0.5, omax * 1.5 + 1e-6

    def px(t):
        return x0 + (t - tmin) / (tmax - tmin) * (x1 - x0)

    def py(o):
        return y1 - (o - omin) / (omax - omin) * (y1 - y0)

    b = [f"<text x='24' y='24' {SANS} font-size='14' font-weight='700' "
         f"fill='{INK}'>Ne proxy vs branch divergence</text>",
         f"<line x1='{x0}' y1='{y1}' x2='{x1}' y2='{y1}' stroke='{DIM}'/>",
         f"<line x1='{x0}' y1='{y0}' x2='{x0}' y2='{y1}' stroke='{DIM}'/>",
         f"<text x='{(x0+x1)/2:.0f}' y='{H-10}' {SANS} font-size='11' "
         f"fill='{DIM}' text-anchor='middle'>mean terminal branch length t</text>",
         f"<text x='16' y='{(y0+y1)/2:.0f}' {SANS} font-size='11' fill='{DIM}' "
         f"transform='rotate(-90 16 {(y0+y1)/2:.0f})' text-anchor='middle'>"
         f"pooled &#969;</text>"]
    for s, t, o in points:
        if o != o:
            continue
        b.append(f"<circle cx='{px(t):.1f}' cy='{py(o):.1f}' r='4' fill='{ACC}'/>")
        b.append(f"<text x='{px(t)+7:.1f}' y='{py(o)+3:.1f}' {MONO} font-size='9' "
                 f"fill='{INK}'>{s}</text>")
    for frac in (0, 0.5, 1.0):
        b.append(f"<text x='{px(tmin+frac*(tmax-tmin)):.0f}' y='{y1+16}' {MONO} "
                 f"font-size='9' fill='{DIM}' text-anchor='middle'>"
                 f"{tmin+frac*(tmax-tmin):.2f}</text>")
        b.append(f"<text x='{x0-8}' y='{py(omin+frac*(omax-omin))+3:.0f}' {MONO} "
                 f"font-size='9' fill='{DIM}' text-anchor='end'>"
                 f"{omin+frac*(omax-omin):.3f}</text>")
    path.write_text(_svg(W, H, "".join(b)), encoding="utf-8")


def to_png(svg: Path) -> bool:
    png = svg.with_suffix(".png")
    for cmd in (["rsvg-convert", "-o", str(png), str(svg)],
                ["inkscape", str(svg), "--export-type=png",
                 f"--export-filename={png}"],
                ["cairosvg", str(svg), "-o", str(png)]):
        if shutil.which(cmd[0]):
            if subprocess.run(cmd, capture_output=True).returncode == 0:
                return True
    return False


# --------------------------------------------------------------------------
def run_stage6(stage5_dir: Path, out_dir: Path) -> Path:
    parser = Path(__file__).with_name("06_parse_codeml_output.py")
    if not parser.exists():
        sys.exit(f"--stage5-dir needs {parser.name} next to this script")
    parsed = out_dir / "06_parsed"
    subprocess.run([sys.executable, str(parser), "--stage5-dir",
                    str(stage5_dir), "-o", str(parsed)], check=True)
    return parsed


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--records-dir", type=Path,
                     help="dir of Stage 6 *_records.tsv tables")
    src.add_argument("--stage5-dir", type=Path,
                     help="a Stage 5 out-dir; Stage 6 is run on it first")
    ap.add_argument("-o", "--out-dir", type=Path, default=Path("07_ne_proxy"))
    ap.add_argument("--ds-ceiling", type=float, default=DS_CEILING_DEFAULT,
                    help=f"drop a gene's contribution to a lineage when its dS "
                         f"there exceeds this (default {DS_CEILING_DEFAULT})")
    ap.add_argument("--bootstrap", type=int, default=BOOTSTRAP_DEFAULT,
                    help=f"gene-bootstrap replicates (default {BOOTSTRAP_DEFAULT})")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--png", action="store_true",
                    help="also write PNG copies of the plots")
    args = ap.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    log = args.out_dir / "stage7.log"
    log.write_text("", encoding="utf-8")
    t0 = time.perf_counter()

    rec_dir = (run_stage6(args.stage5_dir, args.out_dir)
               if args.stage5_dir else args.records_dir)

    m0_path = rec_dir / "m0_records.tsv"
    tr_paths = sorted(rec_dir.glob("two_ratio.*_records.tsv"))
    if not tr_paths and not m0_path.exists():
        ap.error(f"no m0_records.tsv or two_ratio.*_records.tsv in {rec_dir}")
    stage_log(log, "stage7", "start", focal_species=len(tr_paths),
              ds_ceiling=args.ds_ceiling, bootstrap=args.bootstrap,
              m0=m0_path.exists())

    omega_m0 = float("nan")
    if m0_path.exists():
        # apply the same dS filter as the per-species numbers, so omega_M0
        # is comparable to them (an unfiltered M0 is dragged around by the
        # handful of mis-aligned genes the ceiling exists to remove)
        m0_all = read_records(m0_path)
        m0_kept = [r for r in m0_all
                   if r.qc != "ds_floor" and r.dS <= args.ds_ceiling]
        omega_m0 = pooled_omega(m0_kept)

    if not tr_paths:
        # M0-only: no per-species proxy is possible, but the genome-wide
        # pooled omega is still a usable number.
        m0_rows = m0_kept
        lo, hi = bootstrap_ci(m0_rows, args.bootstrap, args.seed)
        (args.out_dir / "ne_proxy.tsv").write_text(
            "scope\tomega_pooled\tci_lo\tci_hi\tn_genes\n"
            f"genome_wide_M0\t{omega_m0:.6f}\t{lo:.6f}\t{hi:.6f}\t"
            f"{len({r.gene_id for r in m0_rows})}\n", encoding="utf-8")
        stage_log(log, "stage7", "done", mode="m0_only",
                  omega_M0=f"{omega_m0:.4f}",
                  elapsed_s=f"{time.perf_counter() - t0:.2f}")
        print(f"M0-only: genome-wide pooled omega = {omega_m0:.4f} "
              f"[{lo:.4f}, {hi:.4f}]\n"
              "no two_ratio.*_records.tsv -> no per-species proxy "
              "(run Stage 4/5 with the two_ratio analysis for that).")
        return 0

    # one M0 log-likelihood per gene, for the LRT vs the two-ratio model
    m0_lnl: dict[str, float] = {}
    if m0_path.exists():
        for r in m0_all:
            if r.lnL is not None:
                m0_lnl.setdefault(r.gene_id, r.lnL)
    CHI2_1DF_P05 = 3.841

    proxy_rows = ["species\tomega_pooled\tci_lo\tci_hi\tn_genes\t"
                  "n_excl_ds_floor\tn_excl_ds_ceiling\tmean_omega\t"
                  "median_omega\tmean_t\tomega_M0\tn_lrt\tn_lrt_p05\t"
                  "frac_lrt_p05\tmean_lrt"]
    per_gene = ["species\tgene_id\tN\tS\tdN\tdS\tomega\tt\tqc_flag\tused\texcl_reason"]
    forest, dists, pooled_by_sp, tvst = {}, {}, {}, []

    for p in tr_paths:
        species = p.name[len("two_ratio."):-len("_records.tsv")]
        rows = read_records(p)
        used, annot, n_floor, n_ceil = filter_focal(rows, species,
                                                    args.ds_ceiling)
        if not used:
            print(f"NOTE: {species}: 0 usable genes after filtering",
                  file=sys.stderr)
        om = pooled_omega(used)
        lo, hi = bootstrap_ci(used, args.bootstrap, args.seed)
        oms = [r.omega for r in used]
        mean_o = st.fmean(oms) if oms else float("nan")
        med_o = st.median(oms) if oms else float("nan")
        mean_t = st.fmean([r.t for r in used]) if used else float("nan")

        # LRT: 2-ratio vs M0, per gene (2-ratio adds one omega parameter -> 1 df)
        tr_lnl: dict[str, float] = {}
        for r in rows:
            if r.lnL is not None:
                tr_lnl.setdefault(r.gene_id, r.lnL)
        lrts = [2.0 * (tr_lnl[g] - m0_lnl[g])
                for g in tr_lnl if g in m0_lnl]
        n_lrt = len(lrts)
        n_sig = sum(1 for x in lrts if x > CHI2_1DF_P05)
        frac_sig = n_sig / n_lrt if n_lrt else float("nan")
        mean_lrt = st.fmean(lrts) if lrts else float("nan")

        proxy_rows.append(
            f"{species}\t{om:.6f}\t{lo:.6f}\t{hi:.6f}\t{len(used)}\t"
            f"{n_floor}\t{n_ceil}\t{mean_o:.6f}\t{med_o:.6f}\t{mean_t:.6f}\t"
            f"{omega_m0:.6f}\t{n_lrt}\t{n_sig}\t{frac_sig:.4f}\t{mean_lrt:.4f}")
        for r, is_used, reason in annot:
            per_gene.append(
                f"{species}\t{r.gene_id}\t{r.N:.1f}\t{r.S:.1f}\t{r.dN:.5f}\t"
                f"{r.dS:.5f}\t{r.omega:.5f}\t{r.t:.5f}\t{r.qc}\t"
                f"{int(is_used)}\t{reason}")
        forest[species] = (om, lo, hi)
        dists[species] = annot
        pooled_by_sp[species] = om
        tvst.append((species, mean_t, om))
        print(f"{species}: pooled omega = {om:.4f}  "
              f"[{lo:.4f}, {hi:.4f}]  (n={len(used)}, "
              f"-{n_floor} ds_floor, -{n_ceil} ds>{args.ds_ceiling})"
              + (f"  |  LRT vs M0: {n_sig}/{n_lrt} genes p<0.05 "
                 f"({frac_sig:.1%}, exp. 5%)" if n_lrt else ""))

    (args.out_dir / "ne_proxy.tsv").write_text(
        "\n".join(proxy_rows) + "\n", encoding="utf-8")
    (args.out_dir / "per_gene_omega.tsv").write_text(
        "\n".join(per_gene) + "\n", encoding="utf-8")

    plots = args.out_dir / "plots"
    plots.mkdir(exist_ok=True)
    forest_plot(forest, omega_m0, plots / "ne_proxy_forest.svg")
    dist_plot(dists, pooled_by_sp, plots / "omega_per_gene_dist.svg")
    omega_vs_t_plot(tvst, plots / "omega_vs_divergence.svg")
    n_png = 0
    if args.png:
        for svg in sorted(plots.glob("*.svg")):
            n_png += to_png(svg)
        if n_png == 0:
            print("NOTE: --png set but no rasteriser found "
                  "(rsvg-convert / inkscape / cairosvg). SVGs written only.",
                  file=sys.stderr)

    elapsed = time.perf_counter() - t0
    stage_log(log, "stage7", "done", species=len(tr_paths),
              omega_M0=f"{omega_m0:.4f}", png=n_png, elapsed_s=f"{elapsed:.2f}")
    print(f"\nM0 genome-wide omega = {omega_m0:.4f}")
    print(f"wrote ne_proxy.tsv, per_gene_omega.tsv, plots/ -> {args.out_dir}/ "
          f"[{elapsed:.2f}s]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
