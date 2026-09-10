#!/usr/bin/env python3
"""
03_codon_align.py

Stage 3 of BUSCOmega: turn the Stage 2 per-gene protein + CDS files into
per-gene codon alignments ready for codeml.

Per common gene:
  - strip a trailing stop (`*`) from each protein sequence (MAFFT chokes on
    it and it is not part of the aligned protein)
  - align the proteins with MAFFT (`--auto --inputorder`)
  - project that protein alignment onto the CDS with PAL2NAL
    (`-output paml`), so every codon column corresponds to one aligned
    amino-acid column

PAL2NAL is run **without** `-nogap` and **without** `-nomismatch`:
  - gap columns are kept here and removed once, consistently, inside codeml
    (`cleandata = 1` at Stage 4/5) -- see docs/methods.md "How codeml
    handles gaps".
  - a protein/CDS translation mismatch is a real error (wrong CDS paired,
    frameshift, wrong genetic code). PAL2NAL exits non-zero and writes no
    usable output; this script drops that gene and logs it, rather than
    letting `-nomismatch` silently delete the offending codons.

A gene is dropped (never silently truncated) if MAFFT or PAL2NAL fails,
produces empty output, or the codon alignment has fewer taxa than went in.
Dropping keeps the invariant that every surviving gene still has every
species.

Input:
    <stage2_dir>            the Stage 2 --out-dir, containing:
        prt/<busco_id>.faa
        cds/<busco_id>.fna

Outputs (into --out-dir):
    prot_aln/<busco_id>.aln     MAFFT protein alignment (FASTA)
    codon_aln/<busco_id>.pml    PAL2NAL codon alignment (PAML sequential)
    stage3_manifest.tsv         busco_id  n_taxa  prot_aln_aa  codon_aln_bp
                                gap_frac  mafft_rc  pal2nal_rc
    stage3_dropped.tsv          busco_id  step  reason

Stdlib only. MAFFT and PAL2NAL (Perl) must be on PATH.

Author: Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)
Version: 0.1.0
"""
from __future__ import annotations

__version__ = "0.1.0"
__author__ = "Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)"

import argparse
import concurrent.futures as cf
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path


def stage_log(logpath: Path, stage: str, event: str, **fields) -> None:
    """One timestamped line appended to <out-dir>/stageN.log."""
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")
    kv = "  ".join(f"{k}={v}" for k, v in fields.items())
    with logpath.open("a", encoding="utf-8") as fh:
        fh.write(f"[{ts}] {stage} {event}  {kv}\n")

MAFFT = "mafft"
PAL2NAL = "pal2nal.pl"

# A drop rate above this suggests a systematic input problem (wrong CDS set,
# un-cleaned isoforms, wrong genetic code) rather than a few odd genes.
DROP_RATE_WARN = 0.05


# --------------------------------------------------------------------------
# small stdlib FASTA / PAML helpers
# --------------------------------------------------------------------------
def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Parse FASTA into [(header_without_'>', sequence), ...], order preserved."""
    recs: list[tuple[str, str]] = []
    hdr: str | None = None
    parts: list[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr is not None:
                    recs.append((hdr, "".join(parts)))
                hdr = line[1:].strip()
                parts = []
            elif line.strip():
                parts.append(line.strip())
    if hdr is not None:
        recs.append((hdr, "".join(parts)))
    return recs


def strip_stop(seq: str) -> str:
    """Drop a single trailing stop character ('*' or '.')."""
    return seq[:-1] if seq[-1:] in ("*", ".") else seq


def write_fasta(recs: list[tuple[str, str]], path: Path) -> None:
    path.write_text(
        "".join(f">{h}\n{s}\n" for h, s in recs), encoding="utf-8")


def paml_stats(pml_text: str) -> tuple[int, int, float]:
    """(n_taxa, n_sites, gap_fraction) from a PAML sequential alignment.

    Tolerant: the first non-blank line is the `<ntax> <nsites>` header; any
    later line made only of nucleotide/gap/ambiguity characters is treated
    as sequence, everything else (names, blank lines) is ignored.
    """
    seq_chars = set("ACGTNRYSWKMBDHVacgtnryswkmbdhv-?")
    n_taxa = n_sites = 0
    gaps = total = 0
    header_seen = False
    for raw in pml_text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if not header_seen:
            toks = line.split()
            if len(toks) >= 2 and toks[0].isdigit() and toks[1].isdigit():
                n_taxa, n_sites = int(toks[0]), int(toks[1])
            header_seen = True
            continue
        if line and set(line) <= seq_chars:
            gaps += line.count("-")
            total += len(line)
    gap_frac = gaps / total if total else 0.0
    return n_taxa, n_sites, gap_frac


# --------------------------------------------------------------------------
# external tools -- isolated here so tests can substitute them
# --------------------------------------------------------------------------
def run_mafft(faa_in: Path, aln_out: Path, threads: int) -> tuple[int, str]:
    """MAFFT protein alignment. Returns (returncode, stderr_text)."""
    cmd = [MAFFT, "--auto", "--inputorder", "--thread", str(threads), str(faa_in)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode == 0 and proc.stdout.strip():
        aln_out.write_text(proc.stdout, encoding="utf-8")
    return proc.returncode, proc.stderr


def run_pal2nal(prot_aln: Path, cds: Path, pml_out: Path) -> tuple[int, str]:
    """PAL2NAL codon projection (no -nogap, no -nomismatch).

    Returns (returncode, stderr_text).
    """
    cmd = [PAL2NAL, str(prot_aln), str(cds), "-output", "paml"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode == 0 and proc.stdout.strip():
        pml_out.write_text(proc.stdout, encoding="utf-8")
    return proc.returncode, proc.stderr


# --------------------------------------------------------------------------
# per-gene worker
# --------------------------------------------------------------------------
class GeneResult:
    __slots__ = ("busco_id", "ok", "row", "drop")

    def __init__(self, busco_id, ok, row=None, drop=None):
        self.busco_id = busco_id
        self.ok = ok
        self.row = row              # manifest line (str) when ok
        self.drop = drop            # (step, reason) when not ok


def align_one(busco_id: str, prt: Path, cds: Path,
              prot_aln_dir: Path, codon_aln_dir: Path,
              threads: int) -> GeneResult:
    prot_recs = read_fasta(prt)
    n_in = len(prot_recs)
    if n_in < 3:
        return GeneResult(busco_id, False,
                          drop=("input", f"only {n_in} protein records"))

    work = Path(tempfile.mkdtemp(prefix=f"s3_{busco_id}_"))
    try:
        tmp_faa = work / "prot.faa"
        write_fasta([(h, strip_stop(s)) for h, s in prot_recs], tmp_faa)

        aln_out = prot_aln_dir / f"{busco_id}.aln"
        rc_m, err_m = run_mafft(tmp_faa, aln_out, threads)
        if rc_m != 0 or not aln_out.exists():
            return GeneResult(busco_id, False,
                              drop=("mafft", _oneline(err_m) or f"rc={rc_m}"))

        pml_out = codon_aln_dir / f"{busco_id}.pml"
        rc_p, err_p = run_pal2nal(aln_out, cds, pml_out)
        if rc_p != 0 or not pml_out.exists():
            # PAL2NAL prints the mismatched sequence name(s) to stderr
            aln_out.unlink(missing_ok=True)   # no orphan protein aln for a drop
            return GeneResult(busco_id, False,
                              drop=("pal2nal", _oneline(err_p) or f"rc={rc_p}"))

        n_taxa, n_sites, gap_frac = paml_stats(pml_out.read_text(encoding="utf-8"))
        if n_taxa != n_in:
            aln_out.unlink(missing_ok=True)
            pml_out.unlink(missing_ok=True)
            return GeneResult(
                busco_id, False,
                drop=("pal2nal", f"codon aln has {n_taxa}/{n_in} taxa"))

        aln_recs = read_fasta(aln_out)
        if not aln_recs:
            aln_out.unlink(missing_ok=True)
            pml_out.unlink(missing_ok=True)
            return GeneResult(busco_id, False,
                              drop=("mafft", "empty / unparseable alignment"))
        prot_aa = len(aln_recs[0][1])
        row = (f"{busco_id}\t{n_taxa}\t{prot_aa}\t{n_sites}\t"
               f"{gap_frac:.4f}\t{rc_m}\t{rc_p}")
        return GeneResult(busco_id, True, row=row)
    except Exception as exc:                       # noqa: BLE001 - drop, don't crash the stage
        aln_out = prot_aln_dir / f"{busco_id}.aln"
        (codon_aln_dir / f"{busco_id}.pml").unlink(missing_ok=True)
        aln_out.unlink(missing_ok=True)
        return GeneResult(busco_id, False, drop=("error", repr(exc)[:200]))
    finally:
        shutil.rmtree(work, ignore_errors=True)


def _oneline(text: str) -> str:
    for ln in (text or "").splitlines():
        ln = ln.strip()
        if ln:
            return ln[:200]
    return ""


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------
def check_tools() -> None:
    missing = []
    if shutil.which(MAFFT) is None:
        missing.append(f"  {MAFFT}  -- https://mafft.cbrc.jp/alignment/software/"
                       "  (conda: conda install -c bioconda mafft)")
    if shutil.which(PAL2NAL) is None:
        missing.append(f"  {PAL2NAL}  -- http://www.bork.embl.de/pal2nal/"
                       "  (conda: conda install -c bioconda pal2nal)")
    if missing:
        sys.exit("Stage 3 needs these on PATH:\n" + "\n".join(missing))


def tool_versions() -> dict[str, str]:
    """Best-effort tool versions for the run log (audit trail)."""
    v: dict[str, str] = {}
    try:
        p = subprocess.run([MAFFT, "--version"], capture_output=True, text=True)
        line = (p.stderr or p.stdout).strip().splitlines()[0]
        v["mafft"] = "_".join(line.split())          # no spaces -> clean log line
    except Exception:                                    # noqa: BLE001
        v["mafft"] = "unknown"
    # PAL2NAL has no version flag; record the resolved path instead
    v["pal2nal"] = shutil.which(PAL2NAL) or PAL2NAL
    return v


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage2_dir", type=Path,
                    help="Stage 2 output dir (contains prt/ and cds/)")
    ap.add_argument("-o", "--out-dir", type=Path, default=Path("."))
    ap.add_argument("--jobs", type=int, default=1,
                    help="genes aligned concurrently (default 1)")
    ap.add_argument("--threads", type=int, default=4,
                    help="threads per MAFFT call (default 4)")
    ap.add_argument("--skip-tool-check", action="store_true",
                    help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    if not args.skip_tool_check:
        check_tools()

    prt_dir = args.stage2_dir / "prt"
    cds_dir = args.stage2_dir / "cds"
    if not prt_dir.is_dir() or not cds_dir.is_dir():
        ap.error(f"{args.stage2_dir} must contain prt/ and cds/ "
                 "(the Stage 2 --out-dir)")

    genes = sorted(p.stem for p in prt_dir.glob("*.faa"))
    if not genes:
        ap.error(f"no *.faa files in {prt_dir}")

    prot_aln_dir = args.out_dir / "prot_aln"
    codon_aln_dir = args.out_dir / "codon_aln"
    prot_aln_dir.mkdir(parents=True, exist_ok=True)
    codon_aln_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    log = args.out_dir / "stage3.log"
    log.write_text("", encoding="utf-8")
    vers = tool_versions() if not args.skip_tool_check else {}
    stage_log(log, "stage3", "start", genes=len(genes),
              jobs=args.jobs, threads=args.threads, **vers)

    print(f"{len(genes)} genes  |  jobs={args.jobs}  threads={args.threads}")

    tasks = []
    for g in genes:
        cds = cds_dir / f"{g}.fna"
        if not cds.exists():
            tasks.append(GeneResult(g, False, drop=("input", "no CDS file")))
            continue
        tasks.append((g, prt_dir / f"{g}.faa", cds))

    results: list[GeneResult] = [t for t in tasks if isinstance(t, GeneResult)]
    todo = [t for t in tasks if not isinstance(t, GeneResult)]

    if args.jobs <= 1:
        for g, prt, cds in todo:
            results.append(align_one(g, prt, cds, prot_aln_dir,
                                     codon_aln_dir, args.threads))
    else:
        with cf.ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(align_one, g, prt, cds, prot_aln_dir,
                              codon_aln_dir, args.threads): g
                    for g, prt, cds in todo}
            for fut in cf.as_completed(futs):
                results.append(fut.result())

    results.sort(key=lambda r: r.busco_id)
    ok = [r for r in results if r.ok]
    dropped = [r for r in results if not r.ok]

    (args.out_dir / "stage3_manifest.tsv").write_text(
        "busco_id\tn_taxa\tprot_aln_aa\tcodon_aln_bp\tgap_frac\t"
        "mafft_rc\tpal2nal_rc\n"
        + "".join(r.row + "\n" for r in ok), encoding="utf-8")
    (args.out_dir / "stage3_dropped.tsv").write_text(
        "busco_id\tstep\treason\n"
        + "".join(f"{r.busco_id}\t{r.drop[0]}\t{r.drop[1]}\n" for r in dropped),
        encoding="utf-8")

    elapsed = time.perf_counter() - t0
    stage_log(log, "stage3", "done", aligned=len(ok), dropped=len(dropped),
              elapsed_s=f"{elapsed:.1f}")
    print(f"\naligned {len(ok)} / {len(genes)} genes ({len(dropped)} dropped) "
          f"in {elapsed:.1f}s -> {prot_aln_dir}/ and {codon_aln_dir}/")
    if dropped:
        by_step: dict[str, int] = {}
        for r in dropped:
            by_step[r.drop[0]] = by_step.get(r.drop[0], 0) + 1
        print("  drops by step: "
              + ", ".join(f"{k}={v}" for k, v in sorted(by_step.items())))

    drop_rate = len(dropped) / len(genes) if genes else 0.0
    if drop_rate > DROP_RATE_WARN:
        print(f"\nNOTE: {drop_rate:.0%} of genes dropped at Stage 3. A few is "
              "normal; this many usually means the wrong CDS set, a wrong "
              "genetic code, or proteomes that were not reduced to one isoform "
              "per gene. Check stage3_dropped.tsv.", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
