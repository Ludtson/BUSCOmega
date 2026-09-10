# Method comparison — 3-species pilot

**Question (raised by the PI):** the free-ratio model also produces ω values
per branch — why does BUSCOmega use the two-ratio branch model instead, and
does the choice change the answer?

**Short answer:** the branch model barely moves the per-species ω. What
moves it is gene filtering, and BUSCOmega's `--ds-ceiling` filter is doing
the same job the concatenation-based papers do by pre-filtering genes.

## The numbers

369 BUSCO genes, *A. halleri* / *A. thaliana* / *C. grandiflora*,
`viridiplantae_odb10`. Per-species terminal-branch ω, count-pooled:

| method | codeml model | layout | a_halleri | a_thaliana | c_grandiflora | genes used |
|---|---|---|---|---|---|---|
| M0, per-species (dS-filtered) | model=0 | per-gene, pooled | 0.157 | 0.159 | 0.158 | 367 / 367 / 366 |
| **2-ratio, per-gene pooled — BUSCOmega** | model=2 | per-gene, pooled | **0.163** | **0.170** | **0.151** | 368 / 368 / 306 |
| free-ratio, per-gene pooled | model=1 | per-gene, pooled | 0.163 | 0.167 | 0.153 | 366 / 368 / 328 |
| free-ratio, concatenate — no filter | model=1 | one alignment | 0.201 | 0.171 | 0.160 | 369 |
| free-ratio, concatenate — misalignment-filtered (≈ eLife) | model=1 | one alignment | 0.163 | 0.167 | 0.158 | 364 |
| free-ratio, concatenate — strict over-filter | model=1 | one alignment | 0.163 | 0.179 | 0.000 ⚠ | 187 |

**"Free model" = codeml `model=1`** — every branch gets its own ω. It is a
*model* choice, independent of per-gene vs concatenate. The actual eLife /
Galtier-lab method (bio++ substitution counting, not codeml) is closest in
*layout* to free-ratio on a concatenate, but it uses a counting method and
it **pre-filters genes**. The "no filter" concatenate row here is a
deliberately broken version, to show what skipping QC does.

`genes used` is per species where it varies (the per-gene methods filter
per species); a single number where all species share one alignment.

### Why M0 gives three (slightly) different per-species numbers

M0 fixes ω *per gene* — one value applied to every branch. But the reported
per-species number is a **branch-length-weighted average of the per-gene
ωɡ**: after `dNɡ,X = ωɡ · dSɡ,X`, the pooled formula reduces to a mean of
the ωɡ with weight ≈ `Sɡ · dSɡ,X` (synonymous substitutions on species X's
branch for gene g). Each species has a different divergence profile across
genes, so the weights differ, so the weighted mean comes out slightly
different: a_halleri 0.1571, a_thaliana 0.1590, c_grandiflora 0.1576
(weighted mean of ωɡ: 0.1567 / 0.1584 / 0.1568 — essentially the pooled
value). The spread is ~1 %, inside the bootstrap noise — M0 is correctly
reporting no per-species difference. (A second, tiny effect: the dS filter
drops one extra gene for c_grandiflora, whose longer branch pushes that
gene's dS past 1.5.)

### Why the gene counts differ

Only **5 genes are genuinely bad** — mis-aligned, dS > 1.5 on some branch
(`144680at33090`, `162794at33090`, `165103at33090`, `192456at33090`,
`5305at33090`). Everything else is a usable gene.

The rest of the drops are **model artifacts, not gene-quality problems**:

- The dS filters look at each gene's *estimated dS on that species' branch*,
  and each model estimates dS differently, so different genes cross the
  thresholds.
- **M0** fixes ω tree-wide → dS on a terminal branch ≈ its length →
  almost every gene usable (366–367 / 369).
- **2-ratio / free-ratio for *C. grandiflora*** keep only 306 / 328.
  On a 3-taxon tree the deep outgroup branch cannot be constrained for
  every gene, so ~60–170 genes collapse to dS ≈ 0 and are flagged
  `ds_floor`. **Small-sample artifact** — at 23 taxa every branch has
  neighbours to pin it and the collapse largely goes away.
- The **concatenate** uses one gene set for all three species, so its
  count is a single number. It handles the outgroup badly at 3 taxa in a
  different way: over-filter to only the "perfectly clean" 187 genes and
  *C. grandiflora*'s terminal branch collapses (a rooting confound — with
  only 3 branches, codeml can shift the outgroup's divergence onto the
  internal branch). The per-gene + pool approach does not have this
  problem because each gene's outgroup branch is estimated on its own.

## Reading it

1. **The branch model does not matter.** BUSCOmega's 2-ratio and the
   free-ratio model, both estimated per gene and pooled the same way, agree
   to within ~2 % for every species (0.163/0.170/0.151 vs
   0.163/0.167/0.153). The two-ratio constraint — one background ω shared
   across all non-focal branches — introduces no detectable bias. It is
   chosen because on a *single* ~400-codon gene free-ratio has
   `2·n_taxa − 3` ω parameters (≈ 43 for 23 taxa) and each is noise;
   two-ratio has two, and the background is well constrained by ~40
   branches' worth of substitution, which anchors the focal estimate.
   After pooling over thousands of genes the two converge — which this
   pilot confirms.

2. **Gene filtering matters, and the `--ds-ceiling` filter is load-bearing.**
   Without it: M0 for *C. grandiflora* falls to 0.109 and for *A. halleri*
   rises to 0.183; the free-ratio concatenate for *A. halleri* hits 0.201.
   In each case ~2 mis-aligned genes with dS ≈ 70 dominate the sums — they
   pull the pooled ω down where their branch is short (dS in the
   denominator collapses) and up where it is long (fake nonsynonymous
   differences in the numerator). With the filter, every route lands at
   0.15–0.17.

3. **The eLife / Galtier-lab concatenation approach needs the same QC,
   just earlier.** Bourguignon et al. 2024 (eLife 93629) and Weyna et al.
   (reviewed preprint 100574) concatenate BUSCO genes and estimate ω once
   with bio++/mapnh — but only after removing sequences with > 10 %
   insertions and genes with deviant topology. Skip that (as the
   "concatenate, not pre-filtered" row does) and the same mis-aligned
   genes inflate *A. halleri* to 0.201. BUSCOmega puts the QC at the
   pooling stage (`ds_floor` drop + `--ds-ceiling`) instead of before
   concatenation; the effect is equivalent.

## What this means for Chapter 3

- The reported Nₑ proxy is robust to the estimation model. State this with
  the comparison table as a supplementary item.
- The 3-taxon pilot is the **weakest** test — only three terminal branches,
  no internal branch, so free-ratio and two-ratio have little room to
  differ. Re-run `compare_omega_methods.py` on the 23-taxon output for the
  definitive version; there free-ratio has 43 ω per gene and the
  comparison is meaningful.
- One fix went into the pipeline from this: `07_ne_proxy.py` now applies
  the dS filter to `omega_M0` too, so the M0 baseline (0.158, not the old
  unfiltered 0.132) is computed consistently with the per-species numbers
  and sits sensibly among them.

## Reproduce

```bash
conda activate buscomega
cd analysis/
python compare_omega_methods.py \
    --aln-dir  ../examples/pilot_stage3_out/codon_aln \
    --tree     ../examples/pilot_stage4_out/trees/m0.nwk \
    --stage7-dir ../examples/pilot_stage7_out \
    -o out_compare/
```

~30 min (free-ratio is ~3× codeml's per-gene cost, and the degenerate
outgroup cases iterate hard).
