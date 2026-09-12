# species_tree.py — a codeml-ready species tree

**Optional.** Stage 4 (`04_codeml_control.py`) needs one Newick species
tree. If you already have a clean topology with tip names that match your
proteomes, you do not need this script — hand your tree straight to
Stage 4. Use this when you need to *clean* a tree you have, or *infer* one
you do not.

## What "codeml-ready" means

- **Topology only matters.** codeml re-estimates every branch length, so
  branch lengths, `)support` values, and `[&…]` annotations on your input
  are ignored — this script strips them.
- **Tip names must match** the species names in `common_scos.tsv` (which
  come from your BUSCO output directory names). A published tree that calls
  a tip `Arabidopsis_thaliana` while your data calls it `a_thaliana` will
  fail Stage 4's cross-check.
- **The `<ntax> 1` header line is not part of any tree.** A tree from
  IQ-TREE, from a paper, or from this script is a plain Newick with no such
  line. **Stage 4 adds `<ntax> 1`** (and the ` #1` foreground labels) when
  it writes the per-analysis tree files — that header tells codeml to reuse
  one tree across a batched `ndata` run. Do not add it yourself; this
  script deliberately does not.

## Dependencies

| subcommand | needs |
|---|---|
| `prepare` | Python 3.9+ only (stdlib) |
| `infer` (supermatrix) | Python 3.9+ only (stdlib) |
| `infer --run-iqtree` | IQ-TREE 2 on `PATH` (<https://iqtree.github.io>) |

IQ-TREE is **not** in the `buscomega` conda environment — it is only for
this optional step. Install it separately, or run the printed command on a
cluster.

## `prepare` — clean an existing tree

```bash
species_tree.py prepare published_tree.nwk -o species_tree.nwk \
    --rename name_map.tsv \
    --match-to 01_common_scos/common_scos.tsv
```

- reads any Newick (branch lengths, support, `[&…]`/`[&&NHX]` comments,
  `'quoted names'`, a `#NEXUS` wrapper) and reduces it to a bare topology
- `--rename` takes a tab-separated `old<TAB>new` map, one pair per line, and
  renames tips — this is how you reconcile `Arabidopsis_thaliana` →
  `a_thaliana`
- `--match-to` compares the resulting tips to your data and prints what is
  only in the tree vs only in the data; the command exits non-zero on a
  mismatch so it fails loudly in a script

Output is a plain topology Newick. Feed it to
`04_codeml_control.py --tree`.

### Getting a starting tree for a well-studied clade

For Brassicaceae and other resolved groups, take the topology from a
phylogenomic paper (e.g. Nikolov et al. 2019; Walden et al. 2020; Hendriks
et al. 2023), write or download the Newick, and clean it with `prepare`.
No inference needed, and a citation defends the topology.

## `match` — your tip names don't match your species names

Common with an **OrthoFinder** species tree
(`SpeciesTree_rooted[_node_labels].txt`): OrthoFinder's tips are your input
FASTA filenames with the extension stripped, so whatever you named the
files survives into the tree. If you ran OrthoFinder on
isoform-cleaned proteomes named e.g. `Athaliana_longest_isoforms.faa`, the
tip is `Athaliana_longest_isoforms` — not `Athaliana`, which is what your
`protein/`/`cds/` directories and `common_scos.tsv` use.

```bash
species_tree.py match SpeciesTree_rooted_node_labels.txt \
    --to protein/ -o rename_map.tsv
```

`--to` accepts a directory of FASTA files (species = filename up to the
first dot, same rule as Stage 1/2), a `common_scos.tsv`, or a one-per-line
list. For each tip it tries, in order: exact match, case-insensitive,
stripping a known suffix (`_longest_isoforms`, `_protein`, `_pep`, ... —
extend `COMMON_TIP_SUFFIXES` in the script for your own naming), then a
unique-prefix match if exactly one target is plausible.

**It never silently applies a guess.** Only exact / case-insensitive /
known-suffix matches are written as active `old<TAB>new` lines. Everything
else — a unique-prefix guess, an ambiguous tip with more than one plausible
target, or no match at all — is written as a `#`-commented line you must
read and either uncomment or fix by hand before it does anything:

```
SpeciesA_longest_isoforms	SpeciesA	suffix-stripped
# VERIFY  Foo_v1	Foo	unique-prefix -- only plausible match, not exact or a known suffix.
# NO MATCH   Totally_Unrelated	?
```

If the tree looks like OrthoFinder output (sequential `N0, N1, N2, ...`
internal-node labels, and/or the `SpeciesTree` filename), `match` prints a
one-line note saying so — informational only, it doesn't change what the
command does.

Feed the reviewed map straight to `prepare --rename` — a third column
(the rule, for your own reference) is fine; only the first two are used:

```bash
species_tree.py prepare SpeciesTree_rooted_node_labels.txt \
    -o species_tree.nwk --rename rename_map.tsv --match-to protein/
```

## `infer` — build a tree from your own alignments

For a clade with no published phylogeny.

```bash
# 1. concatenate the Stage 3 protein alignments into a supermatrix
species_tree.py infer 03_alignments/prot_aln -o tree_infer/

# 2a. run the tree search here (needs IQ-TREE 2)
species_tree.py infer 03_alignments/prot_aln -o tree_infer/ --run-iqtree --threads 8

# 2b. or run the printed iqtree2 command elsewhere, then:
species_tree.py prepare tree_infer/iqtree.treefile -o species_tree.nwk \
    --match-to 01_common_scos/common_scos.tsv
```

- `infer` writes `supermatrix.phy` (relaxed PHYLIP) and `partitions.txt`
  (one `AA` charset per gene, RAxML/IQ-TREE style)
- the single-copy-ortholog set is, by construction, clean phylogenetic
  data; a partitioned ML search recovers the backbone in minutes to about
  an hour for a few dozen taxa
- with `--run-iqtree`, the resulting `.treefile` is cleaned to a topology
  and written as `species_tree.nwk`
- if gene-tree discordance is a concern (rapid radiations), infer per-gene
  trees and summarise with ASTRAL instead — outside this script

## Pilot check

- `prepare examples/3species_pilot/species_tree.nwk` → strips the leading
  `3 1` and returns `(c_grandiflora,(a_halleri,a_thaliana));`
- `infer examples/pilot_stage3_out/prot_aln` → a 3-taxon × 187,766-aa
  supermatrix over 369 genes (a 3-taxon tree has only one unrooted
  topology, so the pilot never actually needs inference — this just
  exercises the concat)

## Test

`tests/test_species_tree.py` — the stdlib parts: Newick cleaning (lengths,
support, internal node labels like RAxML/IQ-TREE/OrthoFinder's `N0, N1,
...`, NEXUS wrapper, quoted names), `--rename`, the tip/data check, the
`match` matching rules and its refusal to auto-apply an uncertain guess,
and the supermatrix concatenation. The IQ-TREE call is not exercised.
