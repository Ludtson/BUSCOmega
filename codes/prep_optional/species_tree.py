#!/usr/bin/env python3
"""
species_tree.py  --  get a codeml-ready species tree for BUSCOmega.

Stage 4 (`04_codeml_control.py`) needs a Newick species tree. Only the
topology is used -- codeml re-estimates every branch length -- but the tree
still has to have the right tip names and be readable. This optional helper
covers the two ways you get there.

  prepare   Clean an existing tree (from the literature or a previous
            analysis): strip branch lengths, internal-node support, inline
            [&...] comments and quoting; optionally rename tips to match
            your proteome names; check the tip set against your data.
            Stdlib only.

  infer     You have no tree. Concatenate the Stage 3 protein alignments
            into a supermatrix and (optionally) run IQ-TREE to get one.
            The concat is stdlib; the tree search needs IQ-TREE 2.

Neither subcommand adds the `<ntax> 1` header or the ` #1` foreground
labels -- those are run-specific and Stage 4 writes them. This helper's job
ends at "a plain topology Newick with the correct tip names".

Author: Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)
Version: 0.1.0
"""
from __future__ import annotations

__version__ = "0.1.0"
__author__ = "Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)"

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

ALN_EXTS = (".aln", ".fa", ".faa", ".fasta")

_COMMENT = re.compile(r"\[[^\]]*\]")          # [&label=...], [&&NHX:...]
_BRLEN = re.compile(r":[0-9.eE+\-]+")
# Anything right after a ')' up to the next structural character is an
# INTERNAL node label -- numeric support (`)95`, `)0.98`) or a named
# internal node (`)N11`, as RAxML/IQ-TREE "node_labels" output uses). A
# tip can never immediately follow ')' (a ')' always closes a clade), so
# stripping this unconditionally is safe and never touches a real tip name.
_SUPPORT = re.compile(r"\)[^,():;]+")


# --------------------------------------------------------------------------
# tree cleaning  (stdlib)
# --------------------------------------------------------------------------
def clean_newick(text: str) -> str:
    """Reduce any Newick string to a bare topology: names + parens + commas."""
    s = text.strip()
    # a NEXUS-ish wrapper: take the first parenthesised tree statement
    m = re.search(r"\(.*\)\s*;", s, flags=re.DOTALL)
    if m:
        s = m.group(0)
    s = _COMMENT.sub("", s)
    # quoted tip names may hold spaces -- Newick convention is space<->underscore
    s = re.sub(r"'([^']*)'", lambda m: m.group(1).replace(" ", "_"), s)
    s = re.sub(r'"([^"]*)"', lambda m: m.group(1).replace(" ", "_"), s)
    s = s.replace("\n", "").replace("\r", "").replace("\t", "")
    s = _BRLEN.sub("", s)
    # remove internal support repeatedly (nested clades)
    prev = None
    while prev != s:
        prev, s = s, _SUPPORT.sub(")", s)
    s = re.sub(r"\s+", "", s)
    if not s.endswith(";"):
        s += ";"
    return s


def tip_labels(topology: str) -> list[str]:
    body = topology.strip().rstrip(";")
    return [t for t in re.split(r"[(),]", body) if t]


def apply_rename(topology: str, mapping: dict[str, str]) -> str:
    def repl(m: re.Match) -> str:
        name = m.group(0)
        return mapping.get(name, name)
    # a tip token is bounded by ( , ) ;
    return re.sub(r"(?<=[(,])[^(),;]+(?=[,);])", repl, topology)


def read_rename(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for ln in path.read_text(encoding="utf-8").splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split("\t") if "\t" in ln else ln.split()
        if len(parts) != 2:
            raise SystemExit(f"rename file: expected 'old<TAB>new', got: {ln!r}")
        out[parts[0]] = parts[1]
    return out


def read_species(path: Path) -> set[str]:
    """Species names from a common_scos.tsv header, or a one-per-line list."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return set()
    head = lines[0].split("\t")
    if head and head[0] == "busco_id":
        return set(head[1:])
    return {ln.strip() for ln in lines if ln.strip() and not ln.startswith("#")}


def cmd_prepare(args) -> int:
    topo = clean_newick(args.in_tree.read_text(encoding="utf-8"))
    if args.rename:
        topo = apply_rename(topo, read_rename(args.rename))
    tips = tip_labels(topo)
    dupes = {t for t in tips if tips.count(t) > 1}
    if dupes:
        raise SystemExit(f"duplicate tip labels after cleaning: {sorted(dupes)}")
    if len(tips) < 3:
        raise SystemExit(f"only {len(tips)} tips; need >= 3")

    status = "ok"
    if args.match_to:
        want = read_species(args.match_to)
        if want:
            only_tree = sorted(set(tips) - want)
            only_data = sorted(want - set(tips))
            if only_tree or only_data:
                status = "MISMATCH"
                print("tip set does not match "
                      f"{args.match_to.name}:", file=sys.stderr)
                print(f"  only in tree: {only_tree or '-'}", file=sys.stderr)
                print(f"  only in data: {only_data or '-'}", file=sys.stderr)
                print("  (fix with --rename, or correct the tree)",
                      file=sys.stderr)

    args.out.write_text(topo + "\n", encoding="utf-8")
    print(f"{len(tips)} tips -> {args.out}   [{status}]")
    print("  " + ", ".join(tips))
    return 0 if status == "ok" or not args.match_to else 1


# --------------------------------------------------------------------------
# supermatrix + IQ-TREE  (infer)
# --------------------------------------------------------------------------
def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur = None
    parts: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(parts)
            cur = line[1:].split()[0]
            parts = []
        elif line.strip():
            parts.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(parts)
    return seqs


def build_supermatrix(aln_dir: Path, out_prefix: Path) -> tuple[Path, Path, int, int]:
    files = sorted(p for p in aln_dir.iterdir()
                   if p.suffix.lower() in ALN_EXTS)
    if not files:
        raise SystemExit(f"no alignment files ({'/'.join(ALN_EXTS)}) in {aln_dir}")

    taxa: list[str] | None = None
    blocks: list[tuple[str, int]] = []          # (gene_name, width)
    cat: dict[str, list[str]] = {}
    for f in files:
        seqs = read_fasta(f)
        if not seqs:
            continue
        names = sorted(seqs)
        width = len(next(iter(seqs.values())))
        if any(len(s) != width for s in seqs.values()):
            raise SystemExit(f"{f.name}: not aligned (ragged lengths)")
        if taxa is None:
            taxa = names
            cat = {t: [] for t in taxa}
        if names != taxa:
            missing = set(taxa) - set(names)
            extra = set(names) - set(taxa)
            raise SystemExit(
                f"{f.name}: taxon set differs from the first alignment "
                f"(missing {sorted(missing) or '-'}, extra {sorted(extra) or '-'}). "
                "BUSCOmega Stage 3 output should be uniform; check your inputs.")
        for t in taxa:
            cat[t].append(seqs[t])
        blocks.append((f.stem, width))

    total = sum(w for _, w in blocks)
    phy = out_prefix.with_suffix(".phy")
    with phy.open("w", encoding="utf-8") as fh:
        fh.write(f" {len(taxa)} {total}\n")
        for t in taxa:
            fh.write(f"{t}  {''.join(cat[t])}\n")

    parts = out_prefix.parent / "partitions.txt"
    with parts.open("w", encoding="utf-8") as fh:
        start = 1
        for name, w in blocks:
            fh.write(f"AA, {name} = {start}-{start + w - 1}\n")
            start += w

    return phy, parts, len(taxa), total


def cmd_infer(args) -> int:
    out_prefix = args.out_dir / "supermatrix"
    args.out_dir.mkdir(parents=True, exist_ok=True)
    phy, parts, ntax, nchar = build_supermatrix(args.aln_dir, out_prefix)
    print(f"supermatrix: {ntax} taxa x {nchar} aa  -> {phy}")
    print(f"partitions ({sum(1 for _ in parts.open())} genes) -> {parts}")

    iqcmd = ["iqtree2", "-s", str(phy), "-p", str(parts),
             "-m", "MFP", "--fast", "-B", "1000", "-T", str(args.threads),
             "--prefix", str(args.out_dir / "iqtree")]
    if not args.run_iqtree:
        print("\nnot running the tree search (pass --run-iqtree to run it).")
        print("command:\n  " + " ".join(iqcmd))
        return 0
    if shutil.which("iqtree2") is None:
        raise SystemExit("iqtree2 not on PATH -- install IQ-TREE 2 "
                         "(https://iqtree.github.io), or run the printed "
                         "command elsewhere.")
    print("\nrunning: " + " ".join(iqcmd))
    rc = subprocess.run(iqcmd).returncode
    if rc != 0:
        raise SystemExit(f"iqtree2 exited {rc}")
    treefile = args.out_dir / "iqtree.treefile"
    topo = clean_newick(treefile.read_text(encoding="utf-8"))
    (args.out_dir / "species_tree.nwk").write_text(topo + "\n", encoding="utf-8")
    print(f"\ntopology -> {args.out_dir / 'species_tree.nwk'}")
    print("  " + ", ".join(tip_labels(topo)))
    print("  feed this to 04_codeml_control.py --tree")
    return 0


# --------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("prepare", help="clean an existing Newick tree")
    p.add_argument("in_tree", type=Path)
    p.add_argument("-o", "--out", type=Path, default=Path("species_tree.nwk"))
    p.add_argument("--rename", type=Path,
                   help="tab-separated old<TAB>new tip-name map")
    p.add_argument("--match-to", type=Path,
                   help="common_scos.tsv (or a species list) to check tips against")
    p.set_defaults(func=cmd_prepare)

    q = sub.add_parser("infer", help="build a supermatrix and (optionally) run IQ-TREE")
    q.add_argument("aln_dir", type=Path,
                   help="Stage 3 prot_aln/ directory")
    q.add_argument("-o", "--out-dir", type=Path, default=Path("species_tree_infer"))
    q.add_argument("--run-iqtree", action="store_true")
    q.add_argument("--threads", type=int, default=4)
    q.set_defaults(func=cmd_infer)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
