#!/usr/bin/env python3
"""
run_busco.py

Prerequisite step for BUSCOmega (not one of the 7 core stages): run BUSCO
on a directory of per-species FASTA files against a lineage dataset, in
parallel, with a manifest and a resume mode. Stage 1 starts from BUSCO's
own `full_table.tsv`, so this script's job ends there.

Two run_busco.sh variants existed in the original toolkit (a plain
sequential loop, and one with fixed-size `&`/`wait` background-job
batches). Rewritten in Python instead of adapting either, for the same
reasons the 7 core stages are Python and not shell:
  - a real work queue (ProcessPoolExecutor schedules the next species the
    moment a worker frees up; fixed-size wait-batches idle on the slowest
    job in each batch)
  - a manifest (species, status, elapsed, % complete), not just stdout
  - one species' BUSCO failure drops + logs and the run continues,
    instead of the whole script exiting
  - testable (subprocess is mocked in the test, same pattern as Stage 3/5)

Input:
    fasta_dir     one <species>.{faa,fa,fna,fasta} per species
    --lineage     a pre-downloaded BUSCO lineage dataset directory (BUSCO
                  runs --offline against it; download it yourself once --
                  see docs/install.md)
    --mode        genome | protein | transcriptome  ("protein", matching
                  isoform-cleaned proteomes, is what BUSCOmega expects)

Outputs (into --out-dir):
    <species>/                 BUSCO's own output tree -- full_table.tsv
                               under run_<lineage>/ is what Stage 1 reads
    logs/<species>.log         that species' BUSCO stdout+stderr
    busco_run_manifest.tsv     species, status, elapsed_s, complete_pct
    busco_summary.csv          via prep_optional/busco_summary.py, if present
    busco_run.log              timestamped start/done line, matching the
                               stageN.log convention

Stdlib only, so it runs fine under any env's Python. BUSCO itself must be
on PATH -- install it into its OWN conda env, not `buscomega` (two
different real solve failures doing otherwise; see docs/install.md):
`conda create -n busco -c bioconda -c conda-forge python=3.11 busco>=5.5`,
then `conda activate busco` for this one step.

Author: Adekolá Owoyemi (Casola Lab, ECCB, Texas A&M University)
Version: 0.1.0
"""
from __future__ import annotations

__version__ = "0.1.0"
__author__ = "Adekola Owoyemi (Casola Lab, ECCB, Texas A&M University)"

import argparse
import concurrent.futures as cf
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

FASTA_EXTS = (".faa", ".fa", ".fna", ".fasta")
MODES = ("genome", "protein", "transcriptome")


def stage_log(logpath: Path, stage: str, event: str, **fields) -> None:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")
    kv = "  ".join(f"{k}={v}" for k, v in fields.items())
    with logpath.open("a", encoding="utf-8") as fh:
        fh.write(f"[{ts}] {stage} {event}  {kv}\n")


def check_tools() -> None:
    if shutil.which("busco") is None:
        sys.exit("run_busco.py needs `busco` on PATH. Install it into its "
                 "own env (not buscomega -- see docs/install.md for why):\n"
                 "  conda create -n busco -c bioconda -c conda-forge "
                 "python=3.11 busco>=5.5\n"
                 "then `conda activate busco` before running this script, "
                 "or otherwise put busco on PATH.")


def find_fasta(fasta_dir: Path) -> list[tuple[str, Path]]:
    """[(species, path), ...], sorted, species = filename up to the first dot
    (matches the original toolkit's `basename | cut -d '.' -f1`)."""
    out: dict[str, Path] = {}
    for p in sorted(fasta_dir.iterdir()):
        if p.is_file() and p.suffix.lower() in FASTA_EXTS:
            out.setdefault(p.name.split(".")[0], p)
    return sorted(out.items())


def already_done(out_dir: Path, species: str) -> bool:
    """True if a full_table.tsv already exists for this species (--resume)."""
    return bool(list((out_dir / species).glob("*/full_table.tsv")))


def complete_pct(out_dir: Path, species: str) -> float | None:
    hits = list((out_dir / species).glob("*/full_table.tsv"))
    if not hits:
        return None
    total = complete = 0
    for line in hits[0].read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        total += 1
        if line.split("\t")[1] == "Complete":
            complete += 1
    return 100 * complete / total if total else None


def run_one(species: str, fasta: Path, out_dir: Path, lineage: Path,
           mode: str, threads: int, log_dir: Path) -> tuple[str, str, float, float | None]:
    t0 = time.perf_counter()
    cmd = ["busco", "--offline", "--in", str(fasta), "--cpu", str(threads),
           "--lineage_dataset", str(lineage), "--mode", mode,
           "--out", species, "--out_path", str(out_dir), "--force"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    (log_dir / f"{species}.log").write_text(
        f"$ {' '.join(cmd)}\n\n--- stdout ---\n{proc.stdout}\n"
        f"--- stderr ---\n{proc.stderr}\n", encoding="utf-8")
    elapsed = time.perf_counter() - t0
    ok = proc.returncode == 0 and already_done(out_dir, species)
    return (species, "ok" if ok else "failed", elapsed,
           complete_pct(out_dir, species) if ok else None)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta_dir", type=Path,
                    help="one <species>.{faa,fa,fna,fasta} per species")
    ap.add_argument("-o", "--out-dir", type=Path, default=Path("busco_out"))
    ap.add_argument("--lineage", type=Path, required=True,
                    help="pre-downloaded BUSCO lineage dataset directory "
                         "(BUSCO runs --offline against it)")
    ap.add_argument("--mode", choices=MODES, required=True)
    ap.add_argument("--jobs", type=int, default=1,
                    help="species run concurrently (default 1)")
    ap.add_argument("--threads", type=int, default=4,
                    help="--cpu passed to each BUSCO call (default 4). "
                         "Total cores used = jobs x threads.")
    ap.add_argument("--resume", action="store_true",
                    help="skip a species that already has a full_table.tsv")
    args = ap.parse_args(argv)

    check_tools()
    if not args.lineage.is_dir():
        ap.error(f"--lineage not found: {args.lineage}")
    lineage = args.lineage.resolve()

    species_list = find_fasta(args.fasta_dir)
    if not species_list:
        ap.error(f"no {'/'.join(FASTA_EXTS)} files in {args.fasta_dir}")

    out_dir = args.out_dir.resolve()
    log_dir = out_dir / "logs"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(exist_ok=True)
    log = out_dir / "busco_run.log"
    t0 = time.perf_counter()
    stage_log(log, "busco", "start", species=len(species_list),
              mode=args.mode, jobs=args.jobs, threads=args.threads,
              lineage=lineage.name)

    todo = []
    skipped = []
    for species, fasta in species_list:
        if args.resume and already_done(out_dir, species):
            skipped.append(species)
        else:
            todo.append((species, fasta))
    if skipped:
        print(f"--resume: skipping {len(skipped)} species with an existing "
              f"full_table.tsv: {', '.join(skipped)}")

    print(f"{len(todo)} species to run  |  jobs={args.jobs}  "
          f"threads={args.threads}  mode={args.mode}")

    results: list[tuple[str, str, float, float | None]] = []
    if args.jobs <= 1:
        for species, fasta in todo:
            r = run_one(species, fasta, out_dir, lineage, args.mode,
                       args.threads, log_dir)
            results.append(r)
            print(f"  {r[0]}: {r[1]}  ({r[2]:.0f}s)"
                  + (f"  {r[3]:.1f}% complete" if r[3] is not None else ""))
    else:
        with cf.ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(run_one, s, f, out_dir, lineage, args.mode,
                              args.threads, log_dir): s for s, f in todo}
            for fut in cf.as_completed(futs):
                r = fut.result()
                results.append(r)
                print(f"  {r[0]}: {r[1]}  ({r[2]:.0f}s)"
                      + (f"  {r[3]:.1f}% complete" if r[3] is not None else ""))

    manifest = ["species\tstatus\telapsed_s\tcomplete_pct"]
    for species, status, elapsed, pct in sorted(results):
        pct_str = f"{pct:.2f}" if pct is not None else ""
        manifest.append(f"{species}\t{status}\t{elapsed:.1f}\t{pct_str}")
    for species in skipped:
        pct = complete_pct(out_dir, species)
        pct_str = f"{pct:.2f}" if pct is not None else ""
        manifest.append(f"{species}\tskipped\t0.0\t{pct_str}")
    (out_dir / "busco_run_manifest.tsv").write_text(
        "\n".join(manifest) + "\n", encoding="utf-8")

    n_ok = sum(1 for _, s, _, _ in results if s == "ok")
    n_failed = sum(1 for _, s, _, _ in results if s == "failed")
    if n_failed:
        print(f"\nNOTE: {n_failed} species failed BUSCO -- see "
              f"{log_dir}/<species>.log. Continuing is fine; Stage 1 just "
              f"needs >=2 successful species.", file=sys.stderr)

    summary_script = Path(__file__).with_name("busco_summary.py")
    if summary_script.exists() and (n_ok or skipped):
        subprocess.run([sys.executable, str(summary_script), str(out_dir),
                        "--outfile", str(out_dir / "busco_summary.csv")])

    elapsed = time.perf_counter() - t0
    stage_log(log, "busco", "done", ok=n_ok, failed=n_failed,
              skipped=len(skipped), elapsed_s=f"{elapsed:.1f}")
    print(f"\n{n_ok} ok, {n_failed} failed, {len(skipped)} skipped "
          f"[{elapsed:.1f}s] -> {out_dir}/")
    usable = n_ok + len(skipped)
    if usable < 2:
        print(f"\nNOTE: only {usable} usable species -- Stage 1 needs >=2.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
