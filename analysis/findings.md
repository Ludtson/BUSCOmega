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

| method | a_halleri | a_thaliana | c_grandiflora |
|---|---|---|---|
| M0, per-species (dS-filtered) | 0.157 | 0.159 | 0.158 |
| M0, concatenate | 0.171 | 0.171 | 0.171 |
| **2-ratio, per-gene pooled — BUSCOmega** | **0.163** | **0.170** | **0.151** |
| free-ratio (`model=1`), per-gene pooled | 0.163 | 0.167 | 0.153 |
| free-ratio, concatenate (eLife-style, *not* pre-filtered) | 0.201 | 0.171 | 0.160 |

![method comparison](method_comparison.svg)

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
