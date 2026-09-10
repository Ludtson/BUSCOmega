#!/usr/bin/env python3
"""
02_extract_sco_seqs.py

Stage 2 of BUSCOmega: turn the Stage 1 common single-copy set into per-gene,
multi-species sequence files ready for alignment.

For each common BUSCO gene, Stage 1 has already named the exact protein ID
per species (BUSCO's "Sequence" column). This script:
  - looks up that protein sequence (exact ID match) in the species proteome
  - looks up the matching CDS (exact ID, then a small set of suffix-strip
    rules -- e.g. Phytozome protein IDs end ".p" while the CDS ID does not;
    never a fuzzy guess)
  - checks CDS length ~= 3 x protein length (a wrong-CDS catch; flagged, not
    dropped -- PAL2NAL at Stage 3 is the hard gate)
  - drops a whole gene if ANY species is missing its protein or CDS, to keep
    the invariant that every retained gene has every species

Input:
    common_scos.tsv    (Stage 1)  -- header: busco_id  <sp_A>  <sp_B> ...
    --prt-dir          one <species>.{faa,fa,fasta} per species column
    --cds-dir          one <species>.{fna,fa,fasta} per species column

Outputs (into --out-dir):
    prt/<busco_id>.faa       one record per species, header ">{sp} {busco_id}"
    cds/<busco_id>.fna       same, nucleotides
    stage2_manifest.tsv      long: busco_id  species  protein_id  cds_id
                             prt_len  cds_len  cds_len_ratio  match_rule
                             length_flag
    stage2_dropped.tsv       busco_id  reason

Stdlib only.

Author: Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)
Version: 0.1.0
"""
from __future__ import annotations

__version__ = "0.1.0"
__author__ = "Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)"

import argparse
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def stage_log(logpath: Path, stage: str, event: str, **fields) -> None:
    """One timestamped line appended to <out-dir>/stageN.log."""
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")
    kv = "  ".join(f"{k}={v}" for k, v in fields.items())
    with logpath.open("a", encoding="utf-8") as fh:
        fh.write(f"[{ts}] {stage} {event}  {kv}\n")


PRT_EXTS = (".faa", ".fa", ".fasta")
CDS_EXTS = (".fna", ".fa", ".fasta")
BUILTIN_STRIPS = (".p", ".P", "-p", "-P", ".cds", ".CDS")

# CDS length should be 3 x protein length, +3 if it carries the stop codon.
# Allow that stop plus a little slop for annotation quirks. Bigger deviation
# means the wrong CDS was probably paired.
def length_flag(prt_len: int, cds_len: int) -> int:
    expected = 3 * prt_len
    slop = max(3, int(0.02 * expected)) + 3   # +3 tolerates an optional stop
    return 0 if abs(cds_len - expected) <= slop else 1

# Heuristics for "your proteomes were probably not isoform-cleaned"
DROP_RATE_WARN = 0.05
FLAG_RATE_WARN = 0.10


def read_fasta(path: Path) -> dict[str, str]:
    """Parse a FASTA into {first_header_token: sequence}. Stdlib."""
    seqs: dict[str, str] = {}
    cur_id: str | None = None
    parts: list[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(parts)
                cur_id = line[1:].split()[0]
                parts = []
            elif line.strip():
                parts.append(line.strip())
    if cur_id is not None:
        seqs[cur_id] = "".join(parts)
    return seqs


def find_species_file(directory: Path, species: str, exts) -> Path | None:
    for ext in exts:
        p = directory / f"{species}{ext}"
        if p.exists():
            return p
    return None


def cds_lookup(cds_index: dict[str, str], protein_id: str, extra_strips):
    """Return (cds_id, rule) or (None, None). Exact, then suffix-strip; no fuzzy."""
    if protein_id in cds_index:
        return protein_id, "exact"
    for suf in (*extra_strips, *BUILTIN_STRIPS):
        if suf and protein_id.endswith(suf):
            cand = protein_id[: -len(suf)]
            if cand in cds_index:
                return cand, f"strip:{suf}"
    return None, None


def read_common_scos(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    assert header[0] == "busco_id", f"unexpected header in {path}: {header}"
    species = header[1:]
    rows = []
    for ln in lines[1:]:
        if not ln.strip():
            continue
        cols = ln.split("\t")
        rows.append((cols[0], dict(zip(species, cols[1:]))))
    return species, rows


def prot_len(seq: str) -> int:
    return len(seq.rstrip("*"))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("common_scos", type=Path, help="Stage 1 common_scos.tsv")
    ap.add_argument("--prt-dir", type=Path, required=True,
                    help="dir with one <species>.faa per species column")
    ap.add_argument("--cds-dir", type=Path, required=True,
                    help="dir with one <species>.fna per species column")
    ap.add_argument("-o", "--out-dir", type=Path, default=Path("."))
    ap.add_argument("--cds-suffix-strip", action="append", default=[],
                    metavar="SUFFIX",
                    help="extra suffix to try stripping from a protein ID to "
                         "match its CDS (repeatable; tried before the built-ins)")
    args = ap.parse_args(argv)

    t0 = time.perf_counter()
    species, rows = read_common_scos(args.common_scos)
    print(f"{len(rows)} common single-copy genes across {len(species)} species")

    # locate + load the per-species FASTAs
    prt_seqs: dict[str, dict[str, str]] = {}
    cds_seqs: dict[str, dict[str, str]] = {}
    for sp in species:
        pf = find_species_file(args.prt_dir, sp, PRT_EXTS)
        cf = find_species_file(args.cds_dir, sp, CDS_EXTS)
        if pf is None:
            ap.error(f"no protein FASTA for {sp!r} in {args.prt_dir} "
                     f"(tried {', '.join(sp + e for e in PRT_EXTS)})")
        if cf is None:
            ap.error(f"no CDS FASTA for {sp!r} in {args.cds_dir}")
        prt_seqs[sp] = read_fasta(pf)
        cds_seqs[sp] = read_fasta(cf)
        print(f"  {sp}: {len(prt_seqs[sp])} proteins, {len(cds_seqs[sp])} CDS")

    prt_out = args.out_dir / "prt"
    cds_out = args.out_dir / "cds"
    prt_out.mkdir(parents=True, exist_ok=True)
    cds_out.mkdir(parents=True, exist_ok=True)
    log = args.out_dir / "stage2.log"
    log.write_text("", encoding="utf-8")
    stage_log(log, "stage2", "start", genes_in=len(rows), species=len(species))

    manifest: list[str] = []
    dropped: list[tuple[str, str]] = []
    rule_by_species: dict[str, dict[str, int]] = {sp: {} for sp in species}
    n_flagged = 0

    for busco_id, sp_to_pid in rows:
        gene_prt: dict[str, str] = {}
        gene_cds: dict[str, str] = {}
        gene_rows: list[str] = []
        drop_reason: str | None = None

        for sp in species:
            pid = sp_to_pid[sp]
            pseq = prt_seqs[sp].get(pid)
            if pseq is None:
                drop_reason = f"{sp}:no_protein({pid})"
                break
            cds_id, rule = cds_lookup(cds_seqs[sp], pid, args.cds_suffix_strip)
            if cds_id is None:
                drop_reason = f"{sp}:no_cds({pid})"
                break
            cseq = cds_seqs[sp][cds_id]
            rule_by_species[sp][rule] = rule_by_species[sp].get(rule, 0) + 1

            plen, clen = prot_len(pseq), len(cseq)
            ratio = clen / (3 * plen) if plen else 0.0
            flag = length_flag(plen, clen)
            n_flagged += flag

            gene_prt[sp] = pseq
            gene_cds[sp] = cseq
            gene_rows.append(
                f"{busco_id}\t{sp}\t{pid}\t{cds_id}\t{plen}\t{clen}\t"
                f"{ratio:.4f}\t{rule}\t{flag}")

        if drop_reason:
            dropped.append((busco_id, drop_reason))
            continue

        (prt_out / f"{busco_id}.faa").write_text(
            "".join(f">{sp} {busco_id}\n{gene_prt[sp]}\n" for sp in species),
            encoding="utf-8")
        (cds_out / f"{busco_id}.fna").write_text(
            "".join(f">{sp} {busco_id}\n{gene_cds[sp]}\n" for sp in species),
            encoding="utf-8")
        manifest.extend(gene_rows)

    (args.out_dir / "stage2_manifest.tsv").write_text(
        "busco_id\tspecies\tprotein_id\tcds_id\tprt_len\tcds_len\t"
        "cds_len_ratio\tmatch_rule\tlength_flag\n" + "\n".join(manifest) + "\n",
        encoding="utf-8")
    (args.out_dir / "stage2_dropped.tsv").write_text(
        "busco_id\treason\n" + "".join(f"{b}\t{r}\n" for b, r in dropped),
        encoding="utf-8")

    n_written = len(rows) - len(dropped)
    elapsed = time.perf_counter() - t0
    stage_log(log, "stage2", "done", written=n_written, dropped=len(dropped),
              length_flags=n_flagged, elapsed_s=f"{elapsed:.2f}")
    print(f"\nwrote {n_written} genes ({len(dropped)} dropped) to "
          f"{prt_out}/ and {cds_out}/   [{elapsed:.2f}s]")
    print("CDS ID match rule per species:")
    for sp in species:
        rb = rule_by_species[sp]
        print(f"  {sp}: " + ", ".join(f"{k}={v}" for k, v in sorted(rb.items())))
    if n_flagged:
        print(f"length-flagged gene-species pairs (wrong-CDS suspects): {n_flagged}")

    drop_rate = len(dropped) / len(rows) if rows else 0
    flag_rate = n_flagged / (n_written * len(species)) if n_written else 0
    if drop_rate > DROP_RATE_WARN or flag_rate > FLAG_RATE_WARN:
        print("\nNOTE: high drop/mismatch rate. If your proteomes were not "
              "reduced to one representative isoform per gene, run "
              "protein-preprocessing-isoform-pipeline first and re-run BUSCO "
              "on the cleaned proteomes.", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
