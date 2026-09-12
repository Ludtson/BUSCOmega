#!/usr/bin/env python3
"""
species_tree.py  --  get a codeml-ready species tree for BUSCOmega.

Stage 4 (`04_codeml_control.py`) needs a Newick species tree. Only the
topology is used -- codeml re-estimates every branch length -- but the tree
still has to have the right tip names and be readable. This optional helper
covers the two ways you get there.

  prepare   Clean an existing tree (from the literature, OrthoFinder, or a
            previous analysis): strip branch lengths, internal-node labels
            (numeric support or named, e.g. RAxML/IQ-TREE/OrthoFinder's
            N0, N1, ...), inline [&...] comments and quoting; optionally
            rename tips to match your proteome names; check the tip set
            against your data. Stdlib only.

  match     Your tree's tip names don't match your species names (common
            with OrthoFinder's species tree -- tips are your input FASTA
            filenames, extension stripped, so a suffix like
            "_longest_isoforms" survives into the tree). Suggests a
            --rename map by exact / case-insensitive / suffix-stripped /
            unique-prefix matching against a directory of FASTA files, a
            common_scos.tsv, or a species list. Never auto-applies a
            guess -- anything not a confident match comes back flagged
            AMBIGUOUS or NO MATCH for you to resolve by hand.

  infer     You have no tree. Concatenate the Stage 3 protein alignments
            into a supermatrix and (optionally) run IQ-TREE to get one.
            The concat is stdlib; the tree search needs IQ-TREE 2.

  prune     Drop one or more named tips (e.g. a polyploid species you're
            excluding from Stage 1 with --exclude) and collapse any
            internal node left with a single child. run_buscomega.py's own
            --exclude does this automatically -- use this subcommand
            directly only if you're not going through the orchestrator.

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


def _parse_topology(body: str, i: int = 0):
    """Recursive-descent parse of a bare topology (names + parens + commas,
    no branch lengths/labels -- run clean_newick() first). A node is either
    a tip name (str) or a list of child nodes. Returns (node, next_index)."""
    if body[i] == "(":
        children = []
        i += 1
        while True:
            child, i = _parse_topology(body, i)
            children.append(child)
            if body[i] == ",":
                i += 1
                continue
            if body[i] == ")":
                i += 1
                break
        return children, i
    j = i
    while body[j] not in "(),;":
        j += 1
    return body[i:j], j


def _serialize_topology(node) -> str:
    if isinstance(node, str):
        return node
    return "(" + ",".join(_serialize_topology(c) for c in node) + ")"


def prune_tips(topology: str, drop: set[str]) -> str:
    """Drop one or more named tips from a bare topology, collapsing any
    internal node left with a single child so the result is still a valid
    binary/polytomous topology (codeml/tree tools reject a unary node)."""
    root, _ = _parse_topology(topology.strip().rstrip(";"))

    def _prune(node):
        if isinstance(node, str):
            return None if node in drop else node
        kept = [k for k in (_prune(c) for c in node) if k is not None]
        if not kept:
            return None
        if len(kept) == 1:
            return kept[0]              # unary node -> collapse
        return kept

    result = _prune(root)
    if result is None or isinstance(result, str):
        raise ValueError("prune_tips: fewer than 2 tips remain after pruning")
    return _serialize_topology(result) + ";"


def read_rename(path: Path) -> dict[str, str]:
    """old<TAB>new, one per line. A third+ column (e.g. `match`'s rule
    annotation) is tolerated and ignored, so its output is directly usable
    here without editing -- only the AMBIGUOUS / NO MATCH lines, which
    have no resolvable second column, need a human's attention first."""
    out: dict[str, str] = {}
    for ln in path.read_text(encoding="utf-8").splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split("\t") if "\t" in ln else ln.split()
        if len(parts) < 2:
            raise SystemExit(f"rename file: expected 'old<TAB>new', got: {ln!r}")
        out[parts[0]] = parts[1]
    return out


def read_species(path: Path) -> set[str]:
    """Species names from: a directory of FASTA files (name = filename up to
    the first dot, same rule as Stage 1/2), a common_scos.tsv header, or a
    one-per-line list."""
    if path.is_dir():
        names = set()
        for p in path.iterdir():
            if p.is_file() and p.suffix.lower() in (
                    ".faa", ".fa", ".fna", ".fasta"):
                names.add(p.name.split(".")[0])
        return names
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return set()
    head = lines[0].split("\t")
    if head and head[0] == "busco_id":
        return set(head[1:])
    return {ln.strip() for ln in lines if ln.strip() and not ln.startswith("#")}


# --------------------------------------------------------------------------
# matching tree tips to your actual species names -- never a silent guess
# --------------------------------------------------------------------------
# Suffixes seen in the wild on tree tips that came from a FASTA filename
# used verbatim as the tip label (OrthoFinder's species tree does exactly
# this: the tip is your input filename with its extension stripped, so
# whatever you named the file survives into the tree -- "_longest_isoforms"
# is common because that is what protein-preprocessing-isoform-pipeline
# appends). Not OrthoFinder-specific: the same stripping helps any tree
# whose tips carry a suffix your data's names don't.
COMMON_TIP_SUFFIXES = (
    "_longest_isoforms", "_longest_isoform", ".longest_isoforms",
    "_protein", "_proteins", "_pep", "_cds", "_cleaned",
)


def looks_like_orthofinder(tree_path: Path, raw_text: str) -> bool:
    """Informational only -- never gates behaviour. OrthoFinder's species
    tree (SpeciesTree_rooted[_node_labels].txt) numbers internal nodes
    N0, N1, N2, ... sequentially starting from the root; that pattern plus
    the filename is a strong enough signal to point the user at `match`."""
    if tree_path.name.startswith("SpeciesTree"):
        return True
    labels = re.findall(r"\)(N\d+)(?=[,():;])", raw_text)
    return len(labels) >= 3 and "N0" in labels


def guess_matches(tips: list[str], targets: set[str]
                  ) -> dict[str, tuple[str | None, str, list[str]]]:
    """Per tip: (best_guess_or_None, rule, other_candidates).

    Rules tried in order, first hit wins: exact / case-insensitive /
    suffix-stripped (either side) / unique-prefix. Never guesses when more
    than one target is equally plausible -- that comes back as "ambiguous"
    with the candidate list, for a human to resolve.
    """
    out: dict[str, tuple[str | None, str, list[str]]] = {}
    lower_targets = {t.lower(): t for t in targets}
    for tip in tips:
        if tip in targets:
            out[tip] = (tip, "exact", [])
            continue
        if tip.lower() in lower_targets:
            out[tip] = (lower_targets[tip.lower()], "case-insensitive", [])
            continue

        stripped_hit = None
        for suf in COMMON_TIP_SUFFIXES:
            if tip.endswith(suf) and tip[: -len(suf)] in targets:
                stripped_hit = tip[: -len(suf)]
                break
        if stripped_hit is None:
            for t in targets:
                for suf in COMMON_TIP_SUFFIXES:
                    if t.endswith(suf) and t[: -len(suf)] == tip:
                        stripped_hit = t
                        break
                if stripped_hit:
                    break
        if stripped_hit:
            out[tip] = (stripped_hit, "suffix-stripped", [])
            continue

        prefix_hits = [t for t in targets
                       if t.startswith(tip) or tip.startswith(t)]
        if len(prefix_hits) == 1:
            out[tip] = (prefix_hits[0], "unique-prefix", [])
        elif len(prefix_hits) > 1:
            out[tip] = (None, "ambiguous", sorted(prefix_hits))
        else:
            out[tip] = (None, "no-match", [])
    return out


def cmd_match(args) -> int:
    topo = clean_newick(args.in_tree.read_text(encoding="utf-8"))
    tips = tip_labels(topo)
    targets = read_species(args.to)
    if not targets:
        raise SystemExit(f"no species names found in {args.to}")

    if looks_like_orthofinder(args.in_tree,
                              args.in_tree.read_text(encoding="utf-8")):
        print("Looks like an OrthoFinder species tree (sequential N-numbered "
              "internal nodes, and/or the SpeciesTree filename). OrthoFinder "
              "tip names are your input FASTA filenames with the extension "
              "stripped -- if they don't match your species names, a suffix "
              "from an upstream step (e.g. _longest_isoforms) is the usual "
              "cause.\n")

    guesses = guess_matches(tips, targets)
    confident = {"exact", "case-insensitive", "suffix-stripped"}
    rows = []
    n_auto = n_review = 0
    for tip in tips:
        best, rule, cands = guesses[tip]
        if rule in confident:
            n_auto += 1
            rows.append(f"{tip}\t{best}\t{rule}")
        elif rule == "unique-prefix":
            n_review += 1
            rows.append(f"# VERIFY  {tip}\t{best}\tunique-prefix -- only "
                        f"plausible match, not exact or a known suffix. "
                        f"Uncomment (remove '# VERIFY  ') once you've "
                        f"checked it.")
        elif rule == "ambiguous":
            n_review += 1
            rows.append(f"# AMBIGUOUS  {tip}\t?\tcandidates: "
                        f"{', '.join(cands)}")
        else:
            n_review += 1
            rows.append(f"# NO MATCH   {tip}\t?")

    print(f"{len(tips)} tips, {len(targets)} target species: "
          f"{n_auto} confident (exact/case/suffix), {n_review} need review")
    for r in rows:
        print(f"  {r}")

    if args.out:
        header = ("# species_tree.py match -- old<TAB>new, one per tip.\n"
                  "# Lines starting with # are ignored by --rename.\n"
                  "# VERIFY / AMBIGUOUS / NO MATCH lines are commented out on\n"
                  "# purpose -- none of them is a confident match. Check each,\n"
                  "# then uncomment (or hand-write) a plain 'old<TAB>new' line.\n")
        args.out.write_text(header + "\n".join(rows) + "\n", encoding="utf-8")
        print(f"\n-> {args.out}"
              + (f"  ({n_review} line(s) need your review before --rename)"
                 if n_review else ""))

    return 0 if n_review == 0 else 1


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
                print("  fix with --rename, or correct the tree -- "
                      "`species_tree.py match` will suggest a rename map",
                      file=sys.stderr)

    args.out.write_text(topo + "\n", encoding="utf-8")
    print(f"{len(tips)} tips -> {args.out}   [{status}]")
    print("  " + ", ".join(tips))
    return 0 if status == "ok" or not args.match_to else 1


def cmd_prune(args) -> int:
    topo = clean_newick(args.in_tree.read_text(encoding="utf-8"))
    tips = set(tip_labels(topo))
    drop = set(args.drop)
    missing = drop - tips
    if missing:
        print(f"error: --drop name(s) not in the tree: {sorted(missing)}",
              file=sys.stderr)
        print(f"tree tips: {sorted(tips)}", file=sys.stderr)
        return 1
    try:
        pruned = prune_tips(topo, drop)
    except ValueError as e:
        raise SystemExit(str(e))
    out = args.out or args.in_tree.with_name(args.in_tree.stem + ".pruned.nwk")
    out.write_text(pruned + "\n", encoding="utf-8")
    kept = tip_labels(pruned)
    print(f"dropped {len(drop)} tip(s): {sorted(drop)}")
    print(f"{len(kept)} tips -> {out}")
    print("  " + ", ".join(kept))
    return 0


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

    m = sub.add_parser("match", help="suggest a --rename map against your "
                       "actual species names (never auto-applies a guess)")
    m.add_argument("in_tree", type=Path)
    m.add_argument("--to", type=Path, required=True,
                  help="a directory of FASTA files, a common_scos.tsv, or a "
                       "one-per-line species list -- your ground truth names")
    m.add_argument("-o", "--out", type=Path,
                  help="write the suggested map here (feed it to "
                       "`prepare --rename`); confident matches only need a "
                       "glance, AMBIGUOUS / NO MATCH lines need fixing by hand")
    m.set_defaults(func=cmd_match)

    r = sub.add_parser("prune", help="drop one or more tips from a Newick tree")
    r.add_argument("in_tree", type=Path)
    r.add_argument("--drop", action="append", required=True, metavar="TIP",
                   help="tip name to remove (repeatable)")
    r.add_argument("-o", "--out", type=Path,
                   help="default: <in_tree stem>.pruned.nwk next to in_tree")
    r.set_defaults(func=cmd_prune)

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
