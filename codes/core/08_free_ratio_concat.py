#!/usr/bin/env python3
"""
08_free_ratio_concat.py

Optional 8th analysis for BUSCOmega: free-ratio (codeml model=1) on ONE
concatenated alignment of every gene, rather than pooling per-gene results.
This is the eLife / Galtier-lab-style layout -- one big alignment, one
codeml run, each species' terminal branch read directly off that run.

Not one of the core 7 stages, and not run by default. It exists as a
second, structurally-independent way to get a free-ratio number per
species, for comparing against Stage 7's `omega_free_ratio` (free-ratio,
per-gene pooled) and `omega_pooled` (2-ratio, per-gene pooled). Pilot
validation (analysis/compare_omega_methods.py) found the per-gene-pooled
and concatenate routes agree once the same genes are used -- disagreement
here usually means gene filtering differs between the two, not that the
model itself is unstable.

Because there is only one "gene" (the whole concatenation), there is
nothing left to pool: the omega codeml reports for a species' branch on
that single run IS the number, with no further arithmetic. That is also
why there is no bootstrap CI here -- Stage 7's gene-bootstrap resamples
genes, and a concatenate has exactly one.

Input:
    --aln-dir     Stage 3 codon_aln/ directory (<busco_id>.pml)
    --tree        species tree, Newick; topology only (same rule as
                  every other stage) -- no foreground label needed, model=1
                  gives every branch its own omega already

A gene is included in the concatenation only if its alignment has exactly
the tree's tip set (no more, no fewer) and every taxon's sequence is the
same length -- both are logged and reported, not silently dropped.

Outputs (into --out-dir, default 08_free_ratio_concat/):
    concat.pml               the supermatrix, PAML sequential format
    concat.ctl, concat.mlc   the single codeml run
    free_ratio_concat.tsv    species  omega  N  S  dN  dS  t
    skipped_genes.tsv        busco_id  reason
    stage8.log

Feed free_ratio_concat.tsv to 07_ne_proxy.py's --concat-tsv to add an
omega_free_ratio_concat column to ne_proxy.tsv alongside the per-gene-pooled
numbers.

Stdlib only. codeml (PAML) must be on PATH.

Author: Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)
Version: 0.1.0
"""
from __future__ import annotations

__version__ = "0.1.0"
__author__ = "Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)"

import argparse
import importlib.util
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _load(name: str, filename: str):
    """Modules here are numeric-prefixed (04_..., 06_...) so they can't be
    imported with a normal `import` statement -- same pattern used by
    analysis/compare_omega_methods.py."""
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


CTL = _load("codeml_control", "04_codeml_control.py")
PARSE = _load("parse_codeml", "06_parse_codeml_output.py")


def stage_log(logpath: Path, stage: str, event: str, **fields) -> None:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")
    kv = "  ".join(f"{k}={v}" for k, v in fields.items())
    with logpath.open("a", encoding="utf-8") as fh:
        fh.write(f"[{ts}] {stage} {event}  {kv}\n")


def read_pml(path: Path) -> tuple[list[str], dict[str, str]]:
    """PAML sequential -> (taxa in file order, {taxon: sequence})."""
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    seqchars = set("ACGTNRYSWKMBDHVacgtn-?")
    taxa, seqs, cur = [], {}, None
    for ln in lines[1:]:
        s = ln.strip()
        if set(s) <= seqchars and cur is not None:
            seqs[cur] = seqs.get(cur, "") + s
        else:
            cur = s.split()[0]
            taxa.append(cur)
    return taxa, seqs


def concatenate(aln_dir: Path, dest: Path,
                tips: set[str]) -> tuple[list[str], int, int, list[tuple[str, str]]]:
    """Concatenate every *.pml in aln_dir whose taxon set exactly matches
    `tips` and whose per-taxon sequence lengths agree. Returns
    (taxa_order, n_used, n_sites, skipped) -- skipped is [(busco_id, reason)],
    never silently dropped."""
    genes = sorted(aln_dir.glob("*.pml"))
    taxa_order = sorted(tips)
    cat: dict[str, list[str]] = {t: [] for t in taxa_order}
    used = 0
    skipped: list[tuple[str, str]] = []
    for g in genes:
        taxa, seqs = read_pml(g)
        if set(taxa) != tips:
            skipped.append((g.stem, "tip set does not match the tree"))
            continue
        lengths = {len(v) for v in seqs.values()}
        if len(lengths) != 1:
            skipped.append((g.stem, "unequal sequence lengths within the gene"))
            continue
        for t in taxa_order:
            cat[t].append(seqs[t])
        used += 1
    total = sum(len(x) for x in cat[taxa_order[0]]) if taxa_order else 0
    with dest.open("w", encoding="utf-8") as fh:
        fh.write(f" {len(taxa_order)} {total}\n")
        for t in taxa_order:
            fh.write(f"{t}\n{''.join(cat[t])}\n")
    return taxa_order, used, total, skipped


def run_codeml(workdir: Path, seqfile: str, treefile: str) -> Path:
    workdir.mkdir(parents=True, exist_ok=True)
    ctl = CTL.write_control(seqfile=seqfile, treefile=treefile,
                            outfile="concat.mlc", ndata=1, model=1,
                            nssites="0", cleandata=1)
    (workdir / "concat.ctl").write_text(ctl, encoding="utf-8")
    proc = subprocess.run(["codeml", "concat.ctl"], cwd=workdir,
                          capture_output=True, text=True)
    (workdir / "codeml.log").write_text(
        f"$ codeml concat.ctl\n\n--- stdout ---\n{proc.stdout}\n"
        f"--- stderr ---\n{proc.stderr}\n", encoding="utf-8")
    return workdir / "concat.mlc"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln-dir", type=Path, required=True,
                    help="Stage 3 codon_aln/ directory")
    ap.add_argument("--tree", type=Path, required=True)
    ap.add_argument("-o", "--out-dir", type=Path, default=Path("08_free_ratio_concat"))
    args = ap.parse_args(argv)

    if not args.tree.exists():
        ap.error(f"tree not found: {args.tree}")
    topo = CTL.strip_to_topology(args.tree.read_text(encoding="utf-8"))
    tips = CTL.tip_labels(topo)
    if len(set(tips)) != len(tips):
        ap.error("tree has duplicate tip labels")

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    log = out / "stage8.log"
    log.write_text("", encoding="utf-8")
    t0 = time.perf_counter()
    stage_log(log, "stage8", "start", tips=len(tips), aln_dir=str(args.aln_dir))

    taxa_order, n_used, n_sites, skipped = concatenate(
        args.aln_dir, out / "concat.pml", set(tips))
    (out / "skipped_genes.tsv").write_text(
        "busco_id\treason\n" + "\n".join(f"{g}\t{r}" for g, r in skipped)
        + ("\n" if skipped else ""), encoding="utf-8")
    if n_used < 2:
        ap.error(f"only {n_used} genes had the full {len(tips)}-tip set "
                 f"with matching lengths -- nothing to concatenate "
                 f"(see {out}/skipped_genes.tsv)")
    print(f"concatenated {n_used} genes ({len(skipped)} skipped, "
          f"see skipped_genes.tsv) -> {n_sites} codon sites, "
          f"{len(taxa_order)} taxa")

    header = f"{len(tips)} 1\n"
    tree_path = out / "tree.nwk"
    tree_path.write_text(header + topo + "\n", encoding="utf-8")

    mlc = run_codeml(out, "concat.pml", "tree.nwk")
    records = PARSE.parse_mlc_file(mlc)
    by_tip = {r.tip: r for r in records if r.tip}

    missing = set(tips) - set(by_tip)
    if missing:
        print(f"NOTE: codeml did not report a branch for: "
              f"{', '.join(sorted(missing))} -- check {out}/codeml.log",
              file=sys.stderr)

    rows = ["species\tomega\tN\tS\tdN\tdS\tt"]
    for sp in tips:
        r = by_tip.get(sp)
        if r is None:
            continue
        rows.append(f"{sp}\t{r.omega:.6f}\t{r.N:.2f}\t{r.S:.2f}\t"
                    f"{r.dN:.6f}\t{r.dS:.6f}\t{r.t:.6f}")
    (out / "free_ratio_concat.tsv").write_text(
        "\n".join(rows) + "\n", encoding="utf-8")

    elapsed = time.perf_counter() - t0
    stage_log(log, "stage8", "done", n_used=n_used, n_skipped=len(skipped),
              n_sites=n_sites, species_reported=len(by_tip),
              elapsed_s=f"{elapsed:.2f}")
    print(f"wrote {out}/free_ratio_concat.tsv "
          f"({len(by_tip)}/{len(tips)} species) [{elapsed:.1f}s]")
    if len(by_tip) < len(tips):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
