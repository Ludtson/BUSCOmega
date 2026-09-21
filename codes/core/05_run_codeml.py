#!/usr/bin/env python3
"""
05_run_codeml.py

Stage 5 of BUSCOmega: run codeml for every analysis in the Stage 4 plan.

For each row of `stage4_analyses.tsv` (M0, and one 2-ratio per focal
species) and each batch of `--batch-size` genes:

  1. concatenate the batch's codon alignments into one PAML file
     (`ndata = batch size` -- codeml resets tree, parameters and likelihood
     between datasets, so a batch gives byte-identical results to running
     the genes one at a time; the Stage 6 parser splits on `Dataset N`)
  2. copy the analysis's tree in (it already carries the `<ntax> 1` header
     Stage 4 wrote, so codeml reuses it across the batch)
  3. fill the four placeholders in the Stage 4 control template
  4. run codeml in the batch's own directory, so its fixed-name scratch
     files (`2NG.*`, `rst`, `rst1`, `rub`, `lnf`) cannot collide with a
     concurrent worker

**Success is decided by the output, not the exit code.** codeml in this
PAML build returns non-zero even on clean runs; a batch counts as done only
if its `.mlc` contains one result block per gene. A batch that comes up
short is retried one gene at a time to isolate the failure; genes that
still fail are logged and skipped, the rest are kept.

Input:
    <stage4_dir>       Stage 4 out-dir (stage4_analyses.tsv, ctl/, trees/)
    --aln-dir          Stage 3 codon_aln/ directory (<busco_id>.pml)

Outputs (into --out-dir, default 05_codeml_out/):
    <analysis>/batch_NNNN/   batch.pml, batch.ctl, tree.nwk, batch.mlc,
                             codeml.log  (+ codeml scratch files)
    <analysis>/batch_NNNN/retry/<busco_id>.mlc   salvaged single-gene runs
    stage5_manifest.tsv      analysis batch n_genes n_ok n_datasets
                             status wall_s mlc
    stage5_failed.tsv        analysis busco_id reason
    stage5.log

Stdlib only. codeml (PAML) must be on PATH.

Author: Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)
Version: 0.1.0
"""
from __future__ import annotations

__version__ = "0.1.0"
__author__ = "Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)"

import argparse
import concurrent.futures as cf
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

CODEML = "codeml"

PH_SEQ = "__SEQFILE__"
PH_TREE = "__TREEFILE__"
PH_OUT = "__OUTFILE__"
PH_NDATA = "__NDATA__"

# one of these appears once per dataset in an mlc file
_DATASET = re.compile(r"^Dataset \d+$", re.M)
_CODONML = re.compile(r"^CODONML \(in paml", re.M)
_PAML_VER = re.compile(r"CODONML \(in paml version ([^,)]+)")

FAIL_RATE_WARN = 0.02


def stage_log(logpath: Path, stage: str, event: str, **fields) -> None:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")
    kv = "  ".join(f"{k}={v}" for k, v in fields.items())
    with logpath.open("a", encoding="utf-8") as fh:
        fh.write(f"[{ts}] {stage} {event}  {kv}\n")


# --------------------------------------------------------------------------
def n_datasets(mlc_text: str) -> int:
    """How many completed result blocks an mlc file holds."""
    return max(len(_DATASET.findall(mlc_text)), len(_CODONML.findall(mlc_text)))


def fill_template(template: str, *, seqfile: str, treefile: str,
                  outfile: str, ndata: int) -> str:
    return (template.replace(PH_SEQ, seqfile).replace(PH_TREE, treefile)
            .replace(PH_OUT, outfile).replace(PH_NDATA, str(ndata)))


def concat_pml(genes: list[str], aln_dir: Path, dest: Path) -> None:
    # newline="\n": dest is real alignment data fed straight to PAML's codeml
    # (a Linux C binary) -- Windows-Python text mode would write CRLF here,
    # and codeml's own block parser isn't guaranteed to treat a stray \r as
    # a safe no-op.
    with dest.open("w", encoding="utf-8", newline="\n") as out:
        for g in genes:
            out.write((aln_dir / f"{g}.pml").read_text(encoding="utf-8").rstrip())
            out.write("\n\n")


# --------------------------------------------------------------------------
# the one external call -- isolated so tests can replace it
# --------------------------------------------------------------------------
def run_codeml(ctl_path: Path, workdir: Path) -> tuple[int, str]:
    """Run codeml in `workdir`. Returns (returncode, stdout+stderr).

    The return code is unreliable in PAML 4.10.x; the caller judges success
    from the .mlc contents.
    """
    proc = subprocess.run([CODEML, ctl_path.name], cwd=workdir,
                          capture_output=True, text=True)
    return proc.returncode, proc.stdout + proc.stderr


# --------------------------------------------------------------------------
@dataclass
class BatchResult:
    analysis: str
    batch: str
    n_genes: int
    n_ok: int
    n_datasets: int
    status: str            # ok | partial | failed
    wall_s: float
    mlc: str
    failed: list[tuple[str, str]]   # (busco_id, reason)


def _one_codeml(workdir: Path, template: str, tree_src: Path,
                genes: list[str], aln_dir: Path, out_name: str) -> int:
    """Set up and run one codeml call in `workdir`. Returns dataset count."""
    workdir.mkdir(parents=True, exist_ok=True)
    concat_pml(genes, aln_dir, workdir / "batch.pml")
    # record the gene order so Stage 6 can map `Dataset N` -> busco_id
    (workdir / "genes.txt").write_text(
        "\n".join(genes) + "\n", encoding="utf-8")
    shutil.copyfile(tree_src, workdir / "tree.nwk")
    ctl = fill_template(template, seqfile="batch.pml", treefile="tree.nwk",
                        outfile=out_name, ndata=len(genes))
    (workdir / "batch.ctl").write_text(ctl, encoding="utf-8")
    rc, log = run_codeml(workdir / "batch.ctl", workdir)
    (workdir / "codeml.log").write_text(log, encoding="utf-8")
    mlc = workdir / out_name
    return n_datasets(mlc.read_text(encoding="utf-8")) if mlc.exists() else 0


def run_batch(analysis: str, batch_idx: int, genes: list[str],
              template: str, tree_src: Path, aln_dir: Path,
              out_root: Path, resume: bool = False) -> BatchResult:
    t0 = time.perf_counter()
    bname = f"batch_{batch_idx:04d}"
    wd = out_root / analysis / bname
    failed: list[tuple[str, str]] = []

    if resume:
        mlc = wd / "batch.mlc"
        if mlc.exists() and n_datasets(
                mlc.read_text(encoding="utf-8", errors="replace")) == len(genes):
            return BatchResult(analysis, bname, len(genes), len(genes),
                               len(genes), "ok", 0.0, str(mlc), [])

    try:
        got = _one_codeml(wd, template, tree_src, genes, aln_dir, "batch.mlc")
    except Exception as exc:                       # noqa: BLE001 - report, continue
        got = 0
        (wd).mkdir(parents=True, exist_ok=True)
        (wd / "error.txt").write_text(repr(exc), encoding="utf-8")

    if got == len(genes):
        return BatchResult(analysis, bname, len(genes), len(genes), got, "ok",
                           time.perf_counter() - t0, str(wd / "batch.mlc"), [])

    # short -- retry each gene alone to find the culprit(s)
    n_ok = 0
    retry = wd / "retry"
    for g in genes:
        try:
            g_got = _one_codeml(retry / g, template, tree_src, [g], aln_dir,
                                f"{g}.mlc")
        except Exception as exc:                   # noqa: BLE001
            g_got = 0
            (retry / g).mkdir(parents=True, exist_ok=True)
            (retry / g / "error.txt").write_text(repr(exc), encoding="utf-8")
        if g_got == 1:
            n_ok += 1
        else:
            failed.append((g, "codeml produced no result block"))

    status = "ok" if n_ok == len(genes) else ("partial" if n_ok else "failed")
    # when the retry recovered everything, the per-gene mlc files are the
    # usable output; Stage 6 reads whatever mlc files exist under the analysis
    return BatchResult(analysis, bname, len(genes), n_ok,
                       got if got else n_ok, status,
                       time.perf_counter() - t0,
                       str(wd / "batch.mlc"), failed)


# --------------------------------------------------------------------------
def read_plan(stage4_dir: Path) -> list[dict]:
    path = stage4_dir / "stage4_analyses.tsv"
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    rows = []
    for ln in lines[1:]:
        if ln.strip():
            rows.append(dict(zip(header, ln.split("\t"))))
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage4_dir", type=Path, help="Stage 4 out-dir")
    ap.add_argument("--aln-dir", type=Path, required=True,
                    help="Stage 3 codon_aln/ directory")
    ap.add_argument("-o", "--out-dir", type=Path, default=Path("05_codeml_out"))
    ap.add_argument("--jobs", type=int, default=1,
                    help="codeml batches run concurrently (default 1)")
    ap.add_argument("--batch-size", type=int, default=1,
                    help="genes per codeml call, via PAML ndata (default 1)")
    ap.add_argument("--analyses", default="",
                    help="comma list of analysis names or types to run "
                         "(default: all in the plan)")
    ap.add_argument("--genes-file", type=Path,
                    help="restrict to the busco_ids listed here (one per line)")
    ap.add_argument("--resume", action="store_true",
                    help="skip a batch whose .mlc already has one result "
                         "block per gene (for restarting a long run)")
    ap.add_argument("--skip-tool-check", action="store_true",
                    help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    if not args.skip_tool_check and shutil.which(CODEML) is None:
        sys.exit(f"Stage 5 needs {CODEML} on PATH (PAML - "
                 "http://abacus.gene.ucl.ac.uk/software/paml.html; "
                 "conda: conda install -c bioconda paml)")

    plan = read_plan(args.stage4_dir)
    if args.analyses:
        want = {a.strip() for a in args.analyses.split(",") if a.strip()}
        plan = [r for r in plan if r["analysis"] in want or r["type"] in want]
    if not plan:
        ap.error("no analyses selected")

    genes = sorted(p.stem for p in args.aln_dir.glob("*.pml"))
    if args.genes_file:
        keep = {ln.strip() for ln in
                args.genes_file.read_text(encoding="utf-8").splitlines()
                if ln.strip()}
        genes = [g for g in genes if g in keep]
    if not genes:
        ap.error(f"no .pml genes found in {args.aln_dir}")

    bs = max(1, args.batch_size)
    batches = [genes[i:i + bs] for i in range(0, len(genes), bs)]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    log = args.out_dir / "stage5.log"
    log.write_text("", encoding="utf-8")
    t0 = time.perf_counter()
    stage_log(log, "stage5", "start", analyses=len(plan), genes=len(genes),
              batch_size=bs, batches_per_analysis=len(batches), jobs=args.jobs)
    print(f"{len(plan)} analyses x {len(batches)} batches "
          f"(batch-size {bs}, {len(genes)} genes)  jobs={args.jobs}")

    jobs = []
    for row in plan:
        template = (args.stage4_dir / row["ctl_template"]).read_text(
            encoding="utf-8")
        tree_src = args.stage4_dir / row["treefile"]
        for bi, gsub in enumerate(batches):
            jobs.append((row["analysis"], bi, gsub, template, tree_src))

    results: list[BatchResult] = []
    if args.jobs <= 1:
        for a, bi, gsub, tmpl, tsrc in jobs:
            results.append(run_batch(a, bi, gsub, tmpl, tsrc, args.aln_dir,
                                     args.out_dir, args.resume))
    else:
        with cf.ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = [ex.submit(run_batch, a, bi, gsub, tmpl, tsrc,
                              args.aln_dir, args.out_dir, args.resume)
                    for a, bi, gsub, tmpl, tsrc in jobs]
            for fut in cf.as_completed(futs):
                results.append(fut.result())

    results.sort(key=lambda r: (r.analysis, r.batch))
    man = ["analysis\tbatch\tn_genes\tn_ok\tn_datasets\tstatus\twall_s\tmlc"]
    failed_rows = ["analysis\tbusco_id\treason"]
    n_ok = n_fail = 0
    for r in results:
        man.append(f"{r.analysis}\t{r.batch}\t{r.n_genes}\t{r.n_ok}\t"
                   f"{r.n_datasets}\t{r.status}\t{r.wall_s:.1f}\t{r.mlc}")
        n_ok += r.n_ok
        for g, reason in r.failed:
            failed_rows.append(f"{r.analysis}\t{g}\t{reason}")
            n_fail += 1
    (args.out_dir / "stage5_manifest.tsv").write_text(
        "\n".join(man) + "\n", encoding="utf-8")
    (args.out_dir / "stage5_failed.tsv").write_text(
        "\n".join(failed_rows) + "\n", encoding="utf-8")

    # PAML version, read from the actual output (codeml has no --version)
    paml_ver = "unknown"
    for r in results:
        p = Path(r.mlc)
        if p.exists():
            m = _PAML_VER.search(p.read_text(encoding="utf-8", errors="replace"))
            if m:
                paml_ver = m.group(1).strip()
                break

    elapsed = time.perf_counter() - t0
    stage_log(log, "stage5", "done", gene_analyses_ok=n_ok, failed=n_fail,
              paml_version=paml_ver, codeml=shutil.which(CODEML) or CODEML,
              elapsed_s=f"{elapsed:.1f}")
    total = len(plan) * len(genes)
    print(f"\n{n_ok} / {total} gene-analyses done, {n_fail} failed, "
          f"in {elapsed:.1f}s -> {args.out_dir}/")
    by_status: dict[str, int] = {}
    for r in results:
        by_status[r.status] = by_status.get(r.status, 0) + 1
    print("  batches: " + ", ".join(f"{k}={v}" for k, v in sorted(by_status.items())))

    if total and n_fail / total > FAIL_RATE_WARN:
        print(f"\nNOTE: {n_fail / total:.1%} of gene-analyses failed. Check "
              "stage5_failed.tsv and a retry/<gene>/codeml.log - usually a "
              "gene with too few informative sites after cleandata, or a "
              "convergence problem.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
