#!/usr/bin/env python3
"""
01_common_scos.py

Stage 1 of BUSCOmega: from per-species BUSCO output, find the single-copy
orthologs (SCOs) common to every input species -- i.e. the BUSCO genes
reported as Complete (which, in BUSCO's own terms, means complete AND
single-copy) in all species.

Supersedes codes/legacy/merge_busco_table.py. Differences:
  - keeps Status == "Complete" only (single-copy); Duplicated/Fragmented/
    Missing all mean "this species does not contribute this gene"
  - computes the strict intersection here, rather than leaving it to a
    downstream script
  - emits a per-species status summary so an underperforming species is
    visible before it silently shrinks the common set
  - deterministic output filenames (reproducible, CI-diffable)

Input layout:
    <busco_parent>/
        <species_A>/ ... /full_table.tsv
        <species_B>/ ... /full_table.tsv
        ...
BUSCO is expected to have been run in protein mode on cleaned
(one-isoform-per-gene) proteomes, so the "Sequence" column is a protein ID
that will exact-match the proteome in Stage 2.

Outputs (into --out-dir):
    common_scos.tsv            busco_id  <species_A>  <species_B> ...
                               one row per SCO common to all species;
                               cell = that species' matched protein ID
    busco_status_summary.tsv   species  complete duplicated fragmented
                               missing  total  in_common_set
                               sole_block_dup  sole_block_frag
                               sole_block_missing  recoverable

                               sole_block_* / recoverable: genes that are
                               Complete in EVERY OTHER species and are kept
                               out of the common set only by this species'
                               Duplicated / Fragmented / Missing call.
                               "recoverable" = how many genes you would gain
                               back by fixing this one species' assembly or
                               annotation. A species with a low complete
                               count but a low recoverable count is not the
                               problem; a high recoverable count is the
                               signal to swap that genome.

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

STATUSES = ("Complete", "Duplicated", "Fragmented", "Missing")


def find_full_table(species_dir: Path) -> Path:
    """Locate the single full_table.tsv under a species' BUSCO output dir."""
    hits = sorted(species_dir.glob("*/full_table.tsv")) or \
        sorted(species_dir.glob("**/full_table.tsv"))
    if not hits:
        raise FileNotFoundError(f"no full_table.tsv found under {species_dir}")
    if len(hits) > 1:
        raise RuntimeError(
            f"{species_dir}: found {len(hits)} full_table.tsv files "
            f"({', '.join(str(h.relative_to(species_dir)) for h in hits)}). "
            "Clean the BUSCO output directory so there is exactly one."
        )
    return hits[0]


def parse_full_table(path: Path) -> tuple[dict[str, str], dict[str, int], dict[str, str]]:
    """Parse one full_table.tsv.

    Returns:
        complete:      {busco_id: sequence_id}  for Status == Complete only
        counts:        {status: n}   -- row counts across all four statuses
                       (Duplicated genes contribute >1 row)
        status_by_id:  {busco_id: status}  -- one entry per gene
    """
    complete: dict[str, str] = {}
    counts = {s: 0 for s in STATUSES}
    status_by_id: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            busco_id, status = cols[0], cols[1]
            if status not in counts:
                # unexpected status label -> count nothing, warn once per file
                print(f"warning: {path}: unrecognized status {status!r} "
                      f"for {busco_id}", file=sys.stderr)
                continue
            counts[status] += 1
            status_by_id[busco_id] = status
            if status == "Complete":
                seq_id = cols[2] if len(cols) > 2 else ""
                if not seq_id:
                    print(f"warning: {path}: {busco_id} is Complete but has "
                          f"no sequence ID", file=sys.stderr)
                    continue
                complete[busco_id] = seq_id
    return complete, counts, status_by_id


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("busco_parent", type=Path,
                    help="directory with one BUSCO output subdirectory per species")
    ap.add_argument("-o", "--out-dir", type=Path, default=Path("."),
                    help="where to write the two output tables (default: cwd)")
    ap.add_argument("--exclude", action="append", default=[], metavar="SPECIES",
                    help="drop this species from consideration entirely "
                         "(repeatable) -- e.g. a polyploid whose homeologs "
                         "make single-copy orthology a poor fit. Prune the "
                         "matching tip from your tree too (species_tree.py "
                         "prune), or use run_buscomega.py's own --exclude, "
                         "which does both.")
    args = ap.parse_args(argv)

    species_dirs = sorted(p for p in args.busco_parent.iterdir() if p.is_dir()
                          and p.name not in set(args.exclude))
    if len(species_dirs) < 2:
        ap.error(f"need >=2 species subdirectories in {args.busco_parent} "
                 f"after --exclude, found {len(species_dirs)}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    log = args.out_dir / "stage1.log"
    log.write_text("", encoding="utf-8")
    species_names = [d.name for d in species_dirs]
    stage_log(log, "stage1", "start", species=len(species_names))

    per_species_complete: dict[str, dict[str, str]] = {}
    per_species_counts: dict[str, dict[str, int]] = {}
    per_species_status: dict[str, dict[str, str]] = {}
    for d in species_dirs:
        ft = find_full_table(d)
        complete, counts, status_by_id = parse_full_table(ft)
        per_species_complete[d.name] = complete
        per_species_counts[d.name] = counts
        per_species_status[d.name] = status_by_id
        print(f"{d.name}: {counts['Complete']} Complete "
              f"(of {sum(counts.values())} BUSCOs), from {ft}")

    # A high Duplicated fraction usually means a multi-isoform proteome:
    # BUSCO matched two isoforms of the same gene and called it Duplicated.
    for s in species_names:
        st = per_species_status[s]
        n_dup = sum(1 for v in st.values() if v == "Duplicated")
        n_comp = sum(1 for v in st.values() if v == "Complete")
        if n_comp and n_dup / n_comp > 0.05:
            print(f"NOTE: {s} has {n_dup} duplicated BUSCO genes "
                  f"({100 * n_dup / n_comp:.0f}% of its Complete count) -- "
                  f"often a sign the proteome was not reduced to one isoform "
                  f"per gene. Consider isoform-cleaning "
                  f"(protein-preprocessing-isoform-pipeline) and re-running "
                  f"BUSCO on the cleaned proteome.", file=sys.stderr)

    n = len(species_names)
    all_ids = set().union(*(set(s) for s in per_species_status.values()))

    # strict intersection: Complete in every species
    common_sorted = sorted(
        bid for bid in all_ids
        if all(per_species_status[s].get(bid) == "Complete" for s in species_names)
    )
    print(f"\ncommon single-copy set: {len(common_sorted)} genes "
          f"in all {n} species")

    # sole-blocker accounting: a gene Complete in every species BUT ONE is
    # kept out of the common set solely by that one species. That gene is
    # "recoverable" -- fixing that species' assembly/annotation would add it.
    sole_block = {s: {"Duplicated": 0, "Fragmented": 0, "Missing": 0}
                  for s in species_names}
    for bid in all_ids:
        complete_in = [s for s in species_names
                       if per_species_status[s].get(bid) == "Complete"]
        if len(complete_in) == n - 1:
            blocker = next(s for s in species_names if s not in complete_in)
            st = per_species_status[blocker].get(bid, "Missing")
            if st in sole_block[blocker]:
                sole_block[blocker][st] += 1

    common_path = args.out_dir / "common_scos.tsv"
    with common_path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("busco_id\t" + "\t".join(species_names) + "\n")
        for bid in common_sorted:
            row = [bid] + [per_species_complete[s][bid] for s in species_names]
            fh.write("\t".join(row) + "\n")

    summary_path = args.out_dir / "busco_status_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("species\tcomplete\tduplicated\tfragmented\tmissing\ttotal\t"
                 "in_common_set\tsole_block_dup\tsole_block_frag\t"
                 "sole_block_missing\trecoverable\n")
        for s in species_names:
            c = per_species_counts[s]
            total = sum(c.values())
            sb = sole_block[s]
            recoverable = sb["Duplicated"] + sb["Fragmented"] + sb["Missing"]
            fh.write(f"{s}\t{c['Complete']}\t{c['Duplicated']}\t"
                     f"{c['Fragmented']}\t{c['Missing']}\t{total}\t"
                     f"{len(common_sorted)}\t{sb['Duplicated']}\t"
                     f"{sb['Fragmented']}\t{sb['Missing']}\t{recoverable}\n")

    elapsed = time.perf_counter() - t0
    stage_log(log, "stage1", "done", common_scos=len(common_sorted),
              elapsed_s=f"{elapsed:.2f}")
    print(f"wrote {common_path}\nwrote {summary_path}")
    print(f"{len(common_sorted)} common SCOs - done in {elapsed:.2f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
