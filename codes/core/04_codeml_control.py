#!/usr/bin/env python3
"""
04_codeml_control.py

Stage 4 of BUSCOmega: set up the codeml analyses.

Up to three analyses per gene (--analyses picks which):
  - M0          one omega tree-wide  (model=0, NSsites=0)
                -> the genome-wide dN/dS = the Ne proxy, and the null for a
                   likelihood-ratio test
  - 2-ratio     the focal species' terminal branch as foreground `#1`,
                everything else background  (model=2, NSsites=0)
                -> that lineage's own omega. One run per focal species.
                Default headline method: one focal species at a time keeps
                the estimate local to that lineage.
  - free_ratio  every branch gets its own omega  (model=1, NSsites=0).
                One run total (not one per focal species) -- every species'
                terminal-branch omega comes out of that same run, pooled
                across genes by Stage 7 exactly like 2-ratio is. OPT-IN,
                not the default: model=1 has more free parameters than
                2-ratio, so any one gene's per-branch estimate is noisier.
                Pilot validation (analysis/compare_omega_methods.py) found
                free-ratio and 2-ratio agree to ~2% once pooled the same
                way -- use it as a robustness check against 2-ratio, not
                as a replacement for it.

The seqfile changes every gene, but the tree and every model parameter are
constant across all genes of an analysis. So Stage 4 does NOT write a
control file per gene (369 x 4 for the pilot, ~100k for the full run --
that many near-identical tiny files is a filesystem problem, not a
feature). It writes:

    04_codeml/
      trees/
        m0.nwk                     topology only, no labels
        two_ratio/<species>.nwk    that tip tagged ` #1`
      ctl/
        m0.ctl                     one template, model params fixed,
        two_ratio.ctl              seqfile/treefile/outfile/ndata as
                                   __PLACEHOLDER__ tokens
      stage4_analyses.tsv          the plan: one row per analysis ->
                                   (type, model, NSsites, tree, template)
      stage4_params.tsv            the model-parameter choices, for the audit
      stage4.log

Stage 5 reads `stage4_analyses.tsv`, fills the four placeholders per gene
batch (via `write_control()` imported from this module -- not by editing
text), and runs codeml.

Input:
    --tree        species tree, Newick. TOPOLOGY is all that is used --
                  codeml re-estimates every branch length. Unrooted.
    --stage3-dir  (optional) Stage 3 out-dir; if given, the tree's tip set
                  is checked against stage3_manifest.tsv -- a mismatch means
                  the wrong tree and every downstream omega would be wrong.

Stdlib only. No tree library: tip extraction and `#1` tagging are string
operations on a topology, and the topology never changes.

Author: Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)
Version: 0.1.0
"""
from __future__ import annotations

__version__ = "0.1.0"
__author__ = "Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)"

import argparse
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

# Placeholder tokens Stage 5 substitutes per gene batch.
PH_SEQ = "__SEQFILE__"
PH_TREE = "__TREEFILE__"
PH_OUT = "__OUTFILE__"
PH_NDATA = "__NDATA__"

# analysis type -> (model, NSsites)
ANALYSIS_MODELS = {
    "m0": (0, "0"),
    "two_ratio": (2, "0"),
    "free_ratio": (1, "0"),
}


# --------------------------------------------------------------------------
# shared stage log (same shape in every stage script)
# --------------------------------------------------------------------------
def stage_log(logpath: Path, stage: str, event: str, **fields) -> None:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")
    kv = "  ".join(f"{k}={v}" for k, v in fields.items())
    with logpath.open("a", encoding="utf-8") as fh:
        fh.write(f"[{ts}] {stage} {event}  {kv}\n")


# --------------------------------------------------------------------------
# minimal Newick handling -- topology only
# --------------------------------------------------------------------------
_BRLEN = re.compile(r":[0-9.eE+\-]+")
# Anything right after a ')' up to the next structural character is an
# INTERNAL node label -- numeric support (`)95`) or a named internal node
# (`)N11`, as RAxML/IQ-TREE "node_labels" output uses). A tip can never
# immediately follow ')' (a ')' always closes a clade), so stripping this
# unconditionally is safe and never touches a real tip name.
_SUPPORT = re.compile(r"\)[^,():;]+")


def strip_to_topology(newick: str) -> str:
    """Drop branch lengths and internal-node support, keep topology + tips."""
    s = newick.strip()
    # a leading "3 1" style header (seen in some PAML-adjacent files)
    if "\n" in s:
        parts = [ln for ln in s.splitlines() if ln.strip()]
        s = parts[-1]
    s = _BRLEN.sub("", s)
    s = _SUPPORT.sub(")", s)
    s = s.replace(" ", "")
    if not s.endswith(";"):
        s += ";"
    return s


def tip_labels(topology: str) -> list[str]:
    """Tip names from a topology-only Newick, in string order."""
    body = topology.strip().rstrip(";")
    tips: list[str] = []
    for tok in re.split(r"[(),]", body):
        tok = tok.strip()
        if tok:
            tips.append(tok)
    return tips


def label_foreground(topology: str, species: str) -> str:
    """Return `topology` with ` #1` after the tip named `species` (exactly one)."""
    # match the bare tip token, bounded by ( , ) or ;  -- not a substring hit
    pat = re.compile(r"(?<=[(,])" + re.escape(species) + r"(?=[,);])")
    new, n = pat.subn(species + " #1", topology)
    if n != 1:
        raise ValueError(
            f"expected exactly one occurrence of tip {species!r}, found {n}")
    return new


# --------------------------------------------------------------------------
# control-file writer  (imported by Stage 5)
# --------------------------------------------------------------------------
def write_control(*, seqfile: str, treefile: str, outfile: str, ndata,
                  model: int, nssites: str, cleandata: int = 1,
                  codon_freq: int = 2, icode: int = 0,
                  fix_kappa: int = 0, kappa: float = 2.0,
                  fix_omega: int = 0, omega: float = 0.4) -> str:
    """Render a codeml control file. Any of seqfile/treefile/outfile/ndata
    may be a placeholder token; everything else is a concrete value."""
    return f"""      seqfile = {seqfile}
     treefile = {treefile}
      outfile = {outfile}

        noisy = 3      * screen output verbosity
      verbose = 0      * 0: concise output
      runmode = 0      * 0: user tree (topology fixed, codeml fits branch lengths)

      seqtype = 1      * 1: codons
    CodonFreq = {codon_freq}      * 0:1/61  1:F1X4  2:F3X4  3:codon table
        clock = 0      * no clock, unrooted tree
        ndata = {ndata}      * number of gene alignments in this run

        model = {model}      * 0: one omega  1: free-ratio  2: two or more ratios (branch)
      NSsites = {nssites}      * 0: one omega among sites
        icode = {icode}      * 0: universal genetic code

    fix_kappa = {fix_kappa}      * 0: estimate kappa
        kappa = {kappa}      * initial kappa
    fix_omega = {fix_omega}      * 0: estimate omega
        omega = {omega}      * initial omega

    fix_alpha = 1      * no gamma rate variation among sites (M0 / branch models)
        alpha = 0
       Malpha = 0
        ncatG = 3

    fix_blength = 0    * ignore any branch lengths on the input tree; estimate all
       method = 0      * 0: simultaneous optimisation

        getSE = 0
 RateAncestor = 0
    cleandata = {cleandata}      * 1: drop any codon column with a gap/ambiguity, once, for every model
"""


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------
def read_manifest_species(stage3_dir: Path) -> set[str]:
    man = stage3_dir / "stage3_manifest.tsv"
    if not man.exists():
        return set()
    # tip set is per-gene identical (Stage 2 invariant); read it from any
    # codon alignment header instead -- the manifest has no species column
    aln_dir = stage3_dir / "codon_aln"
    for pml in sorted(aln_dir.glob("*.pml"))[:1]:
        sp: set[str] = set()
        lines = pml.read_text(encoding="utf-8").splitlines()
        seqchars = set("ACGTNRYSWKMBDHVacgtn-?")
        for ln in lines[1:]:
            ln = ln.strip()
            if ln and not (set(ln) <= seqchars):
                sp.add(ln.split()[0])
        return sp
    return set()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", type=Path, required=True,
                    help="species tree (Newick); topology is all that is used")
    ap.add_argument("-o", "--out-dir", type=Path, default=Path("04_codeml"))
    ap.add_argument("--stage3-dir", type=Path,
                    help="Stage 3 out-dir; cross-checks tree tips vs the data")
    ap.add_argument("--focal", action="append", default=[], metavar="SPECIES",
                    help="species to get a 2-ratio run (repeatable; "
                         "default: every tip in the tree)")
    ap.add_argument("--analyses", default="m0,two_ratio",
                    help="comma list from {m0,two_ratio,free_ratio} "
                         "(default: m0,two_ratio -- free_ratio is opt-in, "
                         "see the module docstring)")
    ap.add_argument("--codon-freq", type=int, default=2, choices=[0, 1, 2, 3])
    ap.add_argument("--kappa", type=float, default=2.0)
    ap.add_argument("--fix-kappa", action="store_true")
    ap.add_argument("--omega", type=float, default=0.4)
    ap.add_argument("--icode", type=int, default=0)
    ap.add_argument("--cleandata", type=int, default=1, choices=[0, 1],
                    help="1 (default): drop gapped/ambiguous codon columns. "
                         "Set 0 only with a documented reason (e.g. many taxa, "
                         "poor site retention).")
    args = ap.parse_args(argv)

    t0 = time.perf_counter()
    analyses = [a.strip() for a in args.analyses.split(",") if a.strip()]
    bad = set(analyses) - set(ANALYSIS_MODELS)
    if bad:
        ap.error(f"unknown analysis type(s): {', '.join(sorted(bad))}")

    if not args.tree.exists():
        ap.error(f"tree not found: {args.tree}")
    topo = strip_to_topology(args.tree.read_text(encoding="utf-8"))
    tips = tip_labels(topo)
    if len(tips) < 3:
        ap.error(f"tree has {len(tips)} tips; need at least 3")
    if len(set(tips)) != len(tips):
        ap.error("tree has duplicate tip labels")

    if args.stage3_dir:
        data_sp = read_manifest_species(args.stage3_dir)
        if data_sp and data_sp != set(tips):
            only_tree = set(tips) - data_sp
            only_data = data_sp - set(tips)
            ap.error(
                "tree tips do not match the Stage 3 data.\n"
                f"  only in tree: {sorted(only_tree) or '-'}\n"
                f"  only in data: {sorted(only_data) or '-'}\n"
                "Fix the tree before Stage 5 -- a wrong topology biases every "
                "gene's omega.")

    focal = args.focal or tips
    unknown = [s for s in focal if s not in tips]
    if unknown:
        ap.error(f"--focal species not in tree: {', '.join(unknown)}")

    out = args.out_dir
    (out / "trees" / "two_ratio").mkdir(parents=True, exist_ok=True)
    (out / "ctl").mkdir(parents=True, exist_ok=True)
    log = out / "stage4.log"
    log.write_text("", encoding="utf-8")
    stage_log(log, "stage4", "start", tree=args.tree.name, tips=len(tips),
              analyses=",".join(analyses), focal=len(focal),
              cleandata=args.cleandata)

    if args.cleandata == 0:
        print("NOTE: cleandata=0 -- gapped/ambiguous codon columns are KEPT. "
              "Only do this with a recorded reason; M0 and 2-ratio must use "
              "the same setting to stay comparable.", file=sys.stderr)

    # ---- trees ----
    # Leading "<ntax> 1" header: codeml reads one tree per dataset unless
    # ntrees=1 tells it to reuse a single tree for every dataset in an
    # `ndata` batch. Without it, a batched run misreads (every dataset gets
    # tree #1's result) AND codeml exits non-zero. With it, batches are
    # correct and the exit code is meaningful. Verified on PAML 4.10.x.
    header = f"{len(tips)} 1\n"
    (out / "trees" / "m0.nwk").write_text(header + topo + "\n", encoding="utf-8")
    tree_path = {"m0": "trees/m0.nwk"}
    if "two_ratio" in analyses:
        for sp in focal:
            labelled = label_foreground(topo, sp)
            p = out / "trees" / "two_ratio" / f"{sp}.nwk"
            p.write_text(header + labelled + "\n", encoding="utf-8")
            assert labelled.count("#1") == 1
            assert len(tip_labels(labelled.replace(" #1", ""))) == len(tips)
        tree_path["two_ratio"] = "trees/two_ratio/<species>.nwk"

    # ---- control templates ----
    common = dict(cleandata=args.cleandata, codon_freq=args.codon_freq,
                  icode=args.icode, fix_kappa=int(args.fix_kappa),
                  kappa=args.kappa, fix_omega=0, omega=args.omega)
    for a in analyses:
        model, nssites = ANALYSIS_MODELS[a]
        ctl = write_control(seqfile=PH_SEQ, treefile=PH_TREE, outfile=PH_OUT,
                            ndata=PH_NDATA, model=model, nssites=nssites,
                            **common)
        (out / "ctl" / f"{a}.ctl").write_text(ctl, encoding="utf-8")

    # ---- plan table ----
    rows = ["analysis\ttype\tmodel\tnssites\ttreefile\tctl_template"]
    plan_n = 0
    for a in analyses:
        model, nssites = ANALYSIS_MODELS[a]
        if a == "two_ratio":
            for sp in focal:
                rows.append(f"two_ratio.{sp}\ttwo_ratio\t{model}\t{nssites}\t"
                            f"trees/two_ratio/{sp}.nwk\tctl/two_ratio.ctl")
                plan_n += 1
        else:
            rows.append(f"{a}\t{a}\t{model}\t{nssites}\t"
                        f"trees/m0.nwk\tctl/{a}.ctl")
            plan_n += 1
    (out / "stage4_analyses.tsv").write_text("\n".join(rows) + "\n",
                                             encoding="utf-8")

    # ---- params record ----
    (out / "stage4_params.tsv").write_text(
        "param\tvalue\n"
        f"cleandata\t{args.cleandata}\n"
        f"CodonFreq\t{args.codon_freq}\n"
        f"icode\t{args.icode}\n"
        f"fix_kappa\t{int(args.fix_kappa)}\n"
        f"kappa\t{args.kappa}\n"
        f"fix_omega\t0\n"
        f"omega\t{args.omega}\n"
        f"model_M0\t0\n"
        f"model_two_ratio\t2\n"
        f"NSsites\t0\n", encoding="utf-8")

    elapsed = time.perf_counter() - t0
    stage_log(log, "stage4", "done", analyses_planned=plan_n,
              elapsed_s=f"{elapsed:.2f}")
    print(f"{plan_n} analyses planned "
          f"({', '.join(analyses)}; {len(focal)} focal species) -> {out}/")
    print(f"  templates: {', '.join(a + '.ctl' for a in analyses)}")
    print(f"  trees: m0.nwk"
          + (f" + {len(focal)} two_ratio/<sp>.nwk" if 'two_ratio' in analyses
             else ""))
    print(f"  plan: stage4_analyses.tsv   done in {elapsed:.2f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
