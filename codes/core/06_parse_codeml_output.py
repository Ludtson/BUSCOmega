#!/usr/bin/env python3
"""
06_parse_codeml_output.py

Stage 6 of BUSCOmega. Parse codeml `mlc` output into a tidy per-gene,
per-branch table.

Covers the models this pipeline runs -- one-ratio (M0, model=0) and
2-ratio (branch, model=2) -- plus free-ratio (model=1). The 2-ratio and
free-ratio blocks share the "dN & dS for each branch" table, so one parser
path reads both; M0 is read from its "omega (dN/dS) = ..." line.

Why a hand-written parser instead of Bio.Phylo.PAML: the branch tables are
among the less consistently-handled cases across PAML versions, and this
parser is validated directly against real output -- both the historical
free-ratio pilot (`examples/3species_pilot/pilot_free_model.mlc`) and the
real Stage 5 run (M0 + 2-ratio, `examples/pilot_stage5_out/*_records.tsv`).

Handles both:
  - one gene per mlc file (one dir per gene), and
  - multiple "Dataset N" blocks concatenated in one mlc file (a batched
    Stage 5 run, codeml `ndata > 1`).

Usage:
    python 06_parse_codeml_output.py <mlc_file_or_dir> [--gene-id ID] -o out.tsv

Author: Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)
Version: 0.1.0
"""
from __future__ import annotations

__version__ = "0.1.0"
__author__ = "Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)"

import argparse
import csv
import re
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterator

DATASET_SPLIT = re.compile(r"\n\nDataset \d+\n")

BRANCH_TABLE_HEADER = re.compile(
    r"dN & dS for each branch\s*\n\s*\n"
    r"\s*branch\s+t\s+N\s+S\s+dN/dS\s+dN\s+dS\s+N\*dN\s+S\*dS\s*\n\s*\n"
)
BRANCH_ROW = re.compile(
    r"^\s*(?P<branch>\S+)\s+(?P<t>[\d.]+)\s+(?P<N>[\d.]+)\s+(?P<S>[\d.]+)\s+"
    r"(?P<omega>[\d.]+)\s+(?P<dN>[\d.]+)\s+(?P<dS>[\d.]+)\s+"
    r"(?P<NdN>[\d.]+)\s+(?P<SdS>[\d.]+)\s*$",
    re.MULTILINE,
)
SEQ_NUM_NAME = re.compile(r"^#(\d+):\s*(\S+)", re.MULTILINE)
NS_LINE = re.compile(r"ns\s*=\s*(\d+)")

LNL_LINE = re.compile(r"lnL\(ntime:\s*\d+\s*np:\s*(\d+)\):\s*(-?[\d.]+)")
KAPPA_LINE = re.compile(r"kappa \(ts/tv\)\s*=\s*([\d.]+)")
TREE_LEN_LINE = re.compile(r"tree length\s*=\s*([\d.]+)")
ONE_RATIO_LINE = re.compile(r"omega \(dN/dS\)\s*=\s*([\d.]+)")

OMEGA_LO = 1e-4   # PAML's near-zero floor for a degenerate branch
OMEGA_HI = 999.0  # PAML's ceiling for a degenerate branch (dS ~ 0)
DS_FLOOR = 1e-3   # below this, dN/dS is numerically unstable


@dataclass
class BranchRecord:
    gene_id: str
    dataset_index: int
    tip: str | None      # species/tip name if resolved, else None (internal branch)
    branch: str           # raw PAML branch label, e.g. "5..1"
    t: float
    N: float
    S: float
    omega: float
    dN: float
    dS: float
    lnL: float | None
    kappa: float | None
    is_terminal: bool
    qc_flag: str          # "ok" | "omega_boundary" | "ds_floor" | "unmatched_tip"


def _extract_lnl(block: str) -> float | None:
    m = LNL_LINE.search(block)
    return float(m.group(2)) if m else None


def _extract_kappa(block: str) -> float | None:
    m = KAPPA_LINE.search(block)
    return float(m.group(1)) if m else None


def _tip_number_map(block: str) -> dict[str, str]:
    """Map PAML's internal tip number ('1', '2', ...) -> sequence/species name.

    codeml numbers tips by their order of appearance in the alignment file
    (echoed per-block as '#1: name', '#2: name', ...), NOT by position in the
    tree string. Internal nodes are numbered above ntax. This mapping is
    deterministic for a fixed tree + fixed alignment order, so branch labels
    like '5..1' can be resolved to species exactly, instead of guessing from
    omega values (which can tie between branches).
    """
    return {num: name for num, name in SEQ_NUM_NAME.findall(block)}


def _ntax(block: str) -> int | None:
    m = NS_LINE.search(block)
    return int(m.group(1)) if m else None


def _qc_flag(omega: float, dS: float) -> str:
    if omega <= OMEGA_LO or omega >= OMEGA_HI:
        return "omega_boundary"
    if dS < DS_FLOOR:
        return "ds_floor"
    return "ok"


def parse_free_ratio_block(block: str, gene_id: str, dataset_index: int) -> list[BranchRecord]:
    """Parse one 'dN & dS for each branch' block (model=1, free-ratio)."""
    records: list[BranchRecord] = []
    header_m = BRANCH_TABLE_HEADER.search(block)
    if not header_m:
        return records
    tail = block[header_m.end():]
    # Row block ends at the next blank-line-preceded section header
    end_m = re.search(r"\n\s*\n\s*tree length for dN", tail)
    table_text = tail[: end_m.start()] if end_m else tail

    lnl = _extract_lnl(block)
    kappa = _extract_kappa(block)
    tip_names = _tip_number_map(block)   # e.g. {"1": "a_halleri", "2": "a_thaliana", "3": "c_grandiflora"}
    ntax = _ntax(block) or len(tip_names)

    for m in BRANCH_ROW.finditer(table_text):
        omega = float(m.group("omega"))
        dS = float(m.group("dS"))
        branch = m.group("branch")
        # branch label is "parent..child"; child <= ntax means it terminates
        # at a tip (deterministic PAML numbering), not internal-node math.
        child = branch.split("..")[-1]
        tip = tip_names.get(child) if child.isdigit() and int(child) <= ntax else None
        is_terminal = tip is not None
        qc = _qc_flag(omega, dS)
        records.append(
            BranchRecord(
                gene_id=gene_id,
                dataset_index=dataset_index,
                tip=tip,
                branch=m.group("branch"),
                t=float(m.group("t")),
                N=float(m.group("N")),
                S=float(m.group("S")),
                omega=omega,
                dN=float(m.group("dN")),
                dS=dS,
                lnL=lnl,
                kappa=kappa,
                is_terminal=is_terminal,
                qc_flag=qc,
            )
        )
    return records


def parse_one_ratio_block(block: str, gene_id: str, dataset_index: int) -> list[BranchRecord]:
    """Parse an M0 (model=0, single tree-wide omega) block into one record."""
    m = ONE_RATIO_LINE.search(block)
    if not m:
        return []
    omega = float(m.group(1))
    lnl = _extract_lnl(block)
    kappa = _extract_kappa(block)
    return [
        BranchRecord(
            gene_id=gene_id,
            dataset_index=dataset_index,
            tip=None,
            branch="tree-wide",
            t=float("nan"),
            N=float("nan"),
            S=float("nan"),
            omega=omega,
            dN=float("nan"),
            dS=float("nan"),
            lnL=lnl,
            kappa=kappa,
            is_terminal=False,
            qc_flag="ok" if OMEGA_LO < omega < OMEGA_HI else "omega_boundary",
        )
    ]


def iter_gene_blocks(mlc_text: str) -> Iterator[str]:
    """Split a possibly-batched mlc file into one block per gene/dataset."""
    parts = DATASET_SPLIT.split(mlc_text)
    for part in parts:
        if "CODONML" in part:
            yield part


def parse_mlc_file(path: Path, gene_id: str | None = None) -> list[BranchRecord]:
    text = path.read_text(errors="replace")
    gid_base = gene_id or path.stem
    records: list[BranchRecord] = []
    for i, block in enumerate(iter_gene_blocks(text), start=1):
        gid = gid_base if gene_id else f"{gid_base}_ds{i}"
        if "free dN/dS Ratios for branches" in block or "dN & dS for each branch" in block:
            records.extend(parse_free_ratio_block(block, gid, i))
        elif ONE_RATIO_LINE.search(block):
            records.extend(parse_one_ratio_block(block, gid, i))
    return records


def write_tsv(records: list[BranchRecord], out_path: Path) -> None:
    fields = list(asdict(records[0]).keys()) if records else [f.name for f in BranchRecord.__dataclass_fields__.values()]
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in records:
            w.writerow(asdict(r))


def parse_stage5_analysis(analysis_dir: Path) -> list[BranchRecord]:
    """Parse every batch mlc under one Stage 5 analysis dir, mapping each
    `Dataset N` back to its busco_id via the batch's genes.txt. Also picks
    up any salvaged retry/<busco_id>.mlc single-gene runs."""
    records: list[BranchRecord] = []
    for batch_dir in sorted(p for p in analysis_dir.iterdir()
                            if p.is_dir() and p.name.startswith("batch_")):
        mlc = batch_dir / "batch.mlc"
        genes_file = batch_dir / "genes.txt"
        if mlc.exists() and genes_file.exists():
            genes = [g for g in genes_file.read_text(
                encoding="utf-8").splitlines() if g.strip()]
            for rec in parse_mlc_file(mlc):
                idx = rec.dataset_index - 1
                if 0 <= idx < len(genes):
                    rec.gene_id = genes[idx]
                records.append(rec)
        retry = batch_dir / "retry"
        if retry.is_dir():
            for gmlc in sorted(retry.glob("*.mlc")):
                records.extend(parse_mlc_file(gmlc, gene_id=gmlc.stem))
    return records


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", type=Path, nargs="?",
                    help="mlc file, or a directory of <gene_id>/mlc files")
    ap.add_argument("-o", "--out", type=Path, required=True,
                    help="output TSV path (file), or output dir with --stage5-dir")
    ap.add_argument("--mlc-name", default="mlc",
                    help="codeml output filename inside each gene dir (default: mlc)")
    ap.add_argument("--stage5-dir", type=Path,
                    help="a Stage 5 out-dir (05_codeml_out/); parse every "
                         "<analysis>/batch_*/batch.mlc into "
                         "<out>/<analysis>_records.tsv")
    args = ap.parse_args(argv)

    if args.stage5_dir:
        args.out.mkdir(parents=True, exist_ok=True)
        analyses = sorted(p for p in args.stage5_dir.iterdir()
                          if p.is_dir() and not p.name.startswith("."))
        n_written = 0
        for adir in analyses:
            recs = parse_stage5_analysis(adir)
            if not recs:
                print(f"warning: no records for {adir.name}", file=sys.stderr)
                continue
            out = args.out / f"{adir.name}_records.tsv"
            write_tsv(recs, out)
            n_ok = sum(1 for r in recs if r.qc_flag == "ok")
            print(f"{adir.name}: {len(recs)} records ({n_ok} ok) -> {out}")
            n_written += 1
        if n_written == 0:
            print("error: parsed no analyses from "
                  f"{args.stage5_dir}", file=sys.stderr)
            return 1
        return 0

    if args.input is None:
        ap.error("need an mlc file/dir, or --stage5-dir")

    all_records: list[BranchRecord] = []
    if args.input.is_dir():
        for gene_dir in sorted(p for p in args.input.iterdir() if p.is_dir()):
            mlc_path = gene_dir / args.mlc_name
            if mlc_path.exists():
                all_records.extend(parse_mlc_file(mlc_path, gene_id=gene_dir.name))
            else:
                print(f"warning: no {args.mlc_name} in {gene_dir}", file=sys.stderr)
    else:
        all_records.extend(parse_mlc_file(args.input))

    write_tsv(all_records, args.out)
    n_ok = sum(1 for r in all_records if r.qc_flag == "ok")
    print(f"wrote {len(all_records)} branch records ({n_ok} passing QC) -> {args.out}")
    return 0 if all_records else 1


if __name__ == "__main__":
    raise SystemExit(main())
