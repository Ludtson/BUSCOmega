#!/usr/bin/env python3
"""
compare_omega_methods.py  --  dissertation robustness check, NOT part of the
BUSCOmega tool.

Estimates the per-species terminal-branch dN/dS four ways on the same set of
codon alignments, so the Ne proxy can be shown to be (or not to be) robust
to the estimation model:

  1. M0                     one omega tree-wide (the constrained baseline)
  2. 2-ratio, per-gene pooled     BUSCOmega's method (read from ne_proxy.tsv)
  3. free-ratio, per-gene pooled  codeml model=1 on each gene, pool the
                                  focal species' terminal branch across genes
  4. free-ratio on the concatenate   codeml model=1 on all genes joined into
                                  one alignment -> terminal omega directly
                                  (this is the eLife / Galtier-lab style)

Run inside the `buscomega` conda env (needs codeml on PATH).

    python compare_omega_methods.py \
        --aln-dir  ../examples/pilot_stage3_out/codon_aln \
        --tree     ../examples/pilot_stage4_out/trees/m0.nwk \
        --stage7-dir ../examples/pilot_stage7_out \
        -o out_compare/

Stdlib only.

Author: Adekolá Owoyemi (Protein Evolution Lab / Casola Lab, Texas A&M University)
"""
from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
from pathlib import Path

CORE = Path(__file__).resolve().parent.parent / "codes" / "core"


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, CORE / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


CTL = _load("codeml_control", "04_codeml_control.py")
PARSE = _load("parse_codeml", "06_parse_codeml_output.py")
NEPROXY = _load("ne_proxy", "07_ne_proxy.py")


# --------------------------------------------------------------------------
def read_pml(path: Path) -> tuple[list[str], dict[str, str]]:
    """PAML sequential -> (taxa in order, {taxon: sequence})."""
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    seqchars = set("ACGTNRYSWKMBDHVacgtn-?")
    taxa, seqs, cur = [], {}, None
    for ln in lines[1:]:
        s = ln.strip()
        if set(s) <= seqchars and cur is not None:
            seqs[cur] = seqs.get(cur, "") + s
        else:
            cur = s.split()[0]
            taxa.append(cur)
    return taxa, seqs


def concatenate(aln_dir: Path, dest: Path) -> tuple[list[str], int, int]:
    genes = sorted(aln_dir.glob("*.pml"))
    taxa0 = None
    cat: dict[str, list[str]] = {}
    used = 0
    for g in genes:
        taxa, seqs = read_pml(g)
        if taxa0 is None:
            taxa0 = taxa
            cat = {t: [] for t in taxa0}
        if sorted(taxa) != sorted(taxa0):
            continue
        w = len(next(iter(seqs.values())))
        if any(len(v) != w for v in seqs.values()):
            continue
        for t in taxa0:
            cat[t].append(seqs[t])
        used += 1
    total = sum(len(x) for x in cat[taxa0[0]])
    with dest.open("w") as fh:
        fh.write(f" {len(taxa0)} {total}\n")
        for t in taxa0:
            fh.write(f"{t}\n{''.join(cat[t])}\n")
    return taxa0, used, total


def run_codeml(workdir: Path, seqfile: str, treefile: str, ndata: int,
               model: int, out_name: str) -> Path:
    workdir.mkdir(parents=True, exist_ok=True)
    ctl = CTL.write_control(seqfile=seqfile, treefile=treefile,
                            outfile=out_name, ndata=ndata, model=model,
                            nssites="0", cleandata=1)
    (workdir / "codeml.ctl").write_text(ctl)
    subprocess.run(["codeml", "codeml.ctl"], cwd=workdir,
                   capture_output=True, text=True)
    return workdir / out_name


def pooled_by_species(records, species_list, ds_ceiling=1.5):
    """records: list of PARSE.BranchRecord. -> {species: pooled_omega}."""
    rows = [NEPROXY.Row(dict(gene_id=r.gene_id, tip=r.tip or "", t=r.t,
                             N=r.N, S=r.S, omega=r.omega, dN=r.dN, dS=r.dS,
                             qc_flag=r.qc_flag))
            for r in records]
    out = {}
    for sp in species_list:
        used, _, _, _ = NEPROXY.filter_focal(rows, sp, ds_ceiling)
        out[sp] = NEPROXY.pooled_omega(used)
    return out


# --------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln-dir", type=Path, required=True)
    ap.add_argument("--tree", type=Path, required=True,
                    help="topology-only Newick or a '<ntax> 1' header file")
    ap.add_argument("--stage7-dir", type=Path, required=True,
                    help="existing BUSCOmega Stage 7 out-dir (ne_proxy.tsv)")
    ap.add_argument("-o", "--out-dir", type=Path, default=Path("out_compare"))
    ap.add_argument("--batch-size", type=int, default=25)
    args = ap.parse_args(argv)

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    # tree with a "<ntax> 1" header for codeml batch reuse
    topo = CTL.strip_to_topology(args.tree.read_text())
    tips = CTL.tip_labels(topo)
    tree_hdr = out / "tree.nwk"
    tree_hdr.write_text(f"{len(tips)} 1\n{topo}\n")

    # ---- method 1 & 2 : from the existing BUSCOmega run -------------------
    npx = {}
    hdr, *rows = (args.stage7_dir / "ne_proxy.tsv").read_text().splitlines()
    hi = hdr.split("\t")
    for r in rows:
        d = dict(zip(hi, r.split("\t")))
        npx[d["species"]] = d
    m0_val = float(next(iter(npx.values()))["omega_M0"])
    two_ratio = {sp: float(npx[sp]["omega_pooled"]) for sp in npx}
    species = list(npx)

    # ---- method 3 : free-ratio (model=1) per gene, pooled ----------------
    genes = sorted(p.stem for p in args.aln_dir.glob("*.pml"))
    batches = [genes[i:i + args.batch_size]
               for i in range(0, len(genes), args.batch_size)]
    fr_records = []
    for bi, gsub in enumerate(batches):
        wd = out / "free_ratio_pergene" / f"batch_{bi:04d}"
        wd.mkdir(parents=True, exist_ok=True)
        with (wd / "batch.pml").open("w") as fh:
            for g in gsub:
                fh.write((args.aln_dir / f"{g}.pml").read_text().rstrip() + "\n\n")
        (wd / "genes.txt").write_text("\n".join(gsub) + "\n")
        (wd / "tree.nwk").write_text(tree_hdr.read_text())
        run_codeml(wd, "batch.pml", "tree.nwk", len(gsub), 1, "batch.mlc")
        for rec in PARSE.parse_mlc_file(wd / "batch.mlc"):
            idx = rec.dataset_index - 1
            if 0 <= idx < len(gsub):
                rec.gene_id = gsub[idx]
            fr_records.append(rec)
        print(f"  free-ratio batch {bi+1}/{len(batches)}", flush=True)
    free_ratio_pooled = pooled_by_species(fr_records, species)

    # ---- method 4 : free-ratio on the concatenate -----------------------
    taxa, n_used, n_sites = concatenate(args.aln_dir, out / "concat.pml")
    print(f"  concatenate: {n_used} genes, {n_sites} sites")
    wd = out / "free_ratio_concat"
    run_codeml(wd, str((out / "concat.pml").resolve()),
               str(tree_hdr.resolve()), 1, 1, "concat.mlc")
    concat_recs = PARSE.parse_mlc_file(wd / "concat.mlc")
    concat_fr = {r.tip: r.omega for r in concat_recs if r.tip in species}

    # ---- method 5 : M0 on the concatenate (genome-wide cross-check) ------
    wd = out / "m0_concat"
    run_codeml(wd, str((out / "concat.pml").resolve()),
               str(tree_hdr.resolve()), 1, 0, "concat.mlc")
    m0c = PARSE.parse_mlc_file(wd / "concat.mlc")
    m0_concat = m0c[0].omega if m0c else float("nan")

    # ---- table ----------------------------------------------------------
    tbl = ["method\t" + "\t".join(species)]
    tbl.append("M0 (per-gene, pooled -- tree-wide)\t"
               + "\t".join(f"{m0_val:.4f}" for _ in species))
    tbl.append("M0 (concatenate -- tree-wide)\t"
               + "\t".join(f"{m0_concat:.4f}" for _ in species))
    tbl.append("2-ratio (per-gene, pooled) -- BUSCOmega\t"
               + "\t".join(f"{two_ratio[s]:.4f}" for s in species))
    tbl.append("free-ratio (per-gene, pooled)\t"
               + "\t".join(f"{free_ratio_pooled[s]:.4f}" for s in species))
    tbl.append("free-ratio (concatenate) -- eLife-style\t"
               + "\t".join(f"{concat_fr.get(s, float('nan')):.4f}"
                           for s in species))
    (out / "method_comparison.tsv").write_text("\n".join(tbl) + "\n")
    print("\n" + "\n".join(tbl))
    print(f"\n-> {out}/method_comparison.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
