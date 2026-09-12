#!/usr/bin/env python3
"""
run_buscomega.py -- the BUSCOmega orchestrator.

Runs stages 1 -> 7 into one numbered run directory, then writes run.log,
run_summary.tsv, and the regenerable gene-tracking audit.

    run_buscomega.py \
        --busco-dir  <dir with one BUSCO output folder per species> \
        --prt-dir    <dir with one <species>.faa per species> \
        --cds-dir    <dir with one <species>.fna per species> \
        --tree       <species_tree.nwk> \
        -o           <run_dir> \
        [--focal SP ...]   (default: every tip in the tree) \
        [--exclude SP ...] (drops it from Stage 1 AND prunes its tree tip) \
        [--jobs N] [--threads N] [--batch-size N] \
        [--ds-ceiling 1.5] [--bootstrap 1000] [--png] \
        [--cds-suffix-strip .p ...] \
        [--from-stage N]  [--to-stage N]  [--resume]  [--dry-run]

Each stage is the standalone script in codes/core/, called as a subprocess
(so a stage failure stops the run cleanly and can be re-run with
--from-stage). Nothing here re-implements a stage.

Stdlib only.

Author: Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)
Version: 0.1.0
"""
from __future__ import annotations

__version__ = "0.1.0"
__author__ = "Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)"

import argparse
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

CORE = Path(__file__).resolve().parent / "codes" / "core"
PREP = Path(__file__).resolve().parent / "codes" / "prep_optional"
STAGE_DIRS = {
    1: "01_common_scos", 2: "02_sequences", 3: "03_alignments",
    4: "04_codeml", 5: "05_codeml_out", 6: "06_parsed", 7: "07_ne_proxy",
}


def _log(run_dir: Path, msg: str) -> None:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with (run_dir / "run.log").open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def _run(run_dir: Path, stage: int, cmd: list[str], dry: bool) -> None:
    _log(run_dir, f"stage {stage}: {' '.join(cmd)}")
    if dry:
        return
    t0 = time.perf_counter()
    proc = subprocess.run([sys.executable, *cmd], capture_output=True, text=True)
    (run_dir / f"stage{stage}.stdout").write_text(
        proc.stdout + proc.stderr, encoding="utf-8")
    dt = time.perf_counter() - t0
    if proc.returncode != 0:
        _log(run_dir, f"stage {stage}: FAILED (rc={proc.returncode}, {dt:.0f}s) "
                      f"-- see {run_dir}/stage{stage}.stdout")
        sys.stdout.write(proc.stdout[-2000:])
        sys.stderr.write(proc.stderr[-2000:])
        raise SystemExit(proc.returncode)
    _log(run_dir, f"stage {stage}: done ({dt:.0f}s)")


# --------------------------------------------------------------------------
# audit: left-join the per-stage manifests into one gene-tracking table
# --------------------------------------------------------------------------
def _read_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
    if not path.exists():
        return [], []
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return [], []
    return lines[0].split("\t"), [ln.split("\t") for ln in lines[1:] if ln.strip()]


def build_audit(run_dir: Path) -> None:
    d1 = run_dir / STAGE_DIRS[1]
    d2 = run_dir / STAGE_DIRS[2]
    d3 = run_dir / STAGE_DIRS[3]
    d5 = run_dir / STAGE_DIRS[5]
    d7 = run_dir / STAGE_DIRS[7]

    h, rows = _read_tsv(d1 / "common_scos.tsv")
    genes = [r[0] for r in rows]
    if not genes:
        return
    track: dict[str, dict] = {g: {"stage1": "common"} for g in genes}

    _, drop2 = _read_tsv(d2 / "stage2_dropped.tsv")
    for r in drop2:
        track.setdefault(r[0], {})["stage2"] = f"dropped:{r[1]}" if len(r) > 1 else "dropped"
    h2, man2 = _read_tsv(d2 / "stage2_manifest.tsv")
    for g in {r[0] for r in man2}:
        track.setdefault(g, {}).setdefault("stage2", "ok")

    _, drop3 = _read_tsv(d3 / "stage3_dropped.tsv")
    for r in drop3:
        track.setdefault(r[0], {})["stage3"] = (
            f"dropped:{r[1]}:{r[2]}" if len(r) > 2 else "dropped")
    h3, man3 = _read_tsv(d3 / "stage3_manifest.tsv")
    for r in man3:
        track.setdefault(r[0], {}).setdefault("stage3", "aligned")

    _, fail5 = _read_tsv(d5 / "stage5_failed.tsv")
    for r in fail5:                       # analysis, busco_id, reason
        if len(r) >= 2:
            track.setdefault(r[1], {}).setdefault("stage5_failed", set())
            track[r[1]]["stage5_failed"].add(r[0])

    hp, pg = _read_tsv(d7 / "per_gene_omega.tsv")
    if hp:
        gi, ui, si = hp.index("gene_id"), hp.index("used"), hp.index("species")
        used_by_gene: dict[str, set] = {}
        for r in pg:
            if r[ui] == "1":
                used_by_gene.setdefault(r[gi], set()).add(r[si])
        for g, sps in used_by_gene.items():
            track.setdefault(g, {})["stage7_used_in"] = ",".join(sorted(sps))

    cols = ["busco_id", "stage1", "stage2", "stage3", "stage5_failed_for",
            "stage7_used_in"]
    out = ["\t".join(cols)]
    for g in sorted(track):
        t = track[g]
        sf = t.get("stage5_failed")
        out.append("\t".join([
            g, t.get("stage1", "-"), t.get("stage2", "-"), t.get("stage3", "-"),
            ",".join(sorted(sf)) if sf else "-", t.get("stage7_used_in", "-")]))
    (run_dir / "gene_tracking.tsv").write_text("\n".join(out) + "\n",
                                               encoding="utf-8")
    _log(run_dir, f"audit: gene_tracking.tsv ({len(track)} genes)")


def build_summary(run_dir: Path) -> None:
    rows = ["stage\tlog_line"]
    for n in range(1, 8):
        for cand in (run_dir / STAGE_DIRS[n] / f"stage{n}.log",):
            if cand.exists():
                for ln in cand.read_text(encoding="utf-8").splitlines():
                    rows.append(f"{n}\t{ln}")
    (run_dir / "run_summary.tsv").write_text("\n".join(rows) + "\n",
                                             encoding="utf-8")


# --------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--busco-dir", type=Path, required=True)
    ap.add_argument("--prt-dir", type=Path, required=True)
    ap.add_argument("--cds-dir", type=Path, required=True)
    ap.add_argument("--tree", type=Path, required=True)
    ap.add_argument("-o", "--out-dir", type=Path, required=True,
                    help="the run directory (created)")
    ap.add_argument("--focal", action="append", default=[], metavar="SPECIES")
    ap.add_argument("--exclude", action="append", default=[], metavar="SPECIES",
                    help="drop this species from the whole run (repeatable) "
                         "-- e.g. a polyploid whose homeologs make "
                         "single-copy orthology a poor fit. Excludes it from "
                         "Stage 1 AND prunes the matching tip from --tree, "
                         "so the two never drift out of sync.")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--batch-size", type=int, default=1)
    ap.add_argument("--ds-ceiling", type=float, default=1.5)
    ap.add_argument("--bootstrap", type=int, default=1000)
    ap.add_argument("--png", action="store_true")
    ap.add_argument("--cds-suffix-strip", action="append", default=[],
                    metavar="SUFFIX")
    ap.add_argument("--from-stage", type=int, default=1, choices=range(1, 8))
    ap.add_argument("--to-stage", type=int, default=7, choices=range(1, 8))
    ap.add_argument("--resume", action="store_true",
                    help="pass --resume to Stage 5 (skip finished batches)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    run = args.out_dir
    run.mkdir(parents=True, exist_ok=True)
    if args.from_stage == 1 and not args.dry_run:
        (run / "run.log").write_text("", encoding="utf-8")
    _log(run, f"BUSCOmega {__version__}  stages {args.from_stage}-{args.to_stage}"
              f"  jobs={args.jobs} threads={args.threads} "
              f"batch-size={args.batch_size}")

    need_bins = {"mafft", "pal2nal.pl", "codeml"}
    have = {b for b in need_bins if shutil.which(b)}
    _log(run, f"tools on PATH: {', '.join(sorted(have)) or 'none'}"
              + (f"  MISSING: {', '.join(sorted(need_bins - have))}"
                 if need_bins - have else ""))

    d = {n: run / STAGE_DIRS[n] for n in STAGE_DIRS}
    S = args.from_stage, args.to_stage

    def want(n):
        return S[0] <= n <= S[1]

    tree = args.tree
    if args.exclude:
        tree = run / "pruned_tree.nwk"
        if not args.dry_run and (want(1) or want(4)):
            prune_cmd = [sys.executable, str(PREP / "species_tree.py"), "prune",
                        str(args.tree), "-o", str(tree)]
            for sp in args.exclude:
                prune_cmd += ["--drop", sp]
            _log(run, f"pruning --exclude species from tree: "
                      f"{', '.join(args.exclude)}")
            subprocess.run(prune_cmd, check=True)
        _log(run, f"--exclude {args.exclude}: dropped from Stage 1 input "
                  f"and pruned from tree -> {tree}")

    if want(1):
        cmd = [str(CORE / "01_common_scos.py"), str(args.busco_dir),
               "-o", str(d[1])]
        for sp in args.exclude:
            cmd += ["--exclude", sp]
        _run(run, 1, cmd, args.dry_run)
    if want(2):
        cmd = [str(CORE / "02_extract_sco_seqs.py"),
               str(d[1] / "common_scos.tsv"),
               "--prt-dir", str(args.prt_dir), "--cds-dir", str(args.cds_dir),
               "-o", str(d[2])]
        for suf in args.cds_suffix_strip:
            cmd += ["--cds-suffix-strip", suf]
        _run(run, 2, cmd, args.dry_run)
    if want(3):
        _run(run, 3, [str(CORE / "03_codon_align.py"), str(d[2]),
                      "-o", str(d[3]), "--jobs", str(args.jobs),
                      "--threads", str(args.threads)], args.dry_run)
    if want(4):
        cmd = [str(CORE / "04_codeml_control.py"), "--tree", str(tree),
               "--stage3-dir", str(d[3]), "-o", str(d[4])]
        for sp in args.focal:
            cmd += ["--focal", sp]
        _run(run, 4, cmd, args.dry_run)
    if want(5):
        cmd = [str(CORE / "05_run_codeml.py"), str(d[4]),
               "--aln-dir", str(d[3] / "codon_aln"), "-o", str(d[5]),
               "--jobs", str(args.jobs), "--batch-size", str(args.batch_size)]
        if args.resume:
            cmd.append("--resume")
        _run(run, 5, cmd, args.dry_run)
    if want(6):
        _run(run, 6, [str(CORE / "06_parse_codeml_output.py"),
                      "--stage5-dir", str(d[5]), "-o", str(d[6])], args.dry_run)
    if want(7):
        cmd = [str(CORE / "07_ne_proxy.py"), "--records-dir", str(d[6]),
               "-o", str(d[7]), "--ds-ceiling", str(args.ds_ceiling),
               "--bootstrap", str(args.bootstrap)]
        if args.png:
            cmd.append("--png")
        _run(run, 7, cmd, args.dry_run)

    if not args.dry_run:
        build_summary(run)
        build_audit(run)
        _log(run, f"run complete -> {run}/  "
                  f"(ne_proxy.tsv, run_summary.tsv, gene_tracking.tsv, run.log)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
