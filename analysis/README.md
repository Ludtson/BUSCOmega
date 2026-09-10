# analysis/

Dissertation analyses that **use** BUSCOmega but are **not part of the
tool**. Kept here for reproducibility; not covered by the pipeline tests or
the stdlib-only guarantee for the core.

## compare_omega_methods.py

A robustness check: estimates the per-species terminal-branch dN/dS four
ways on the same codon alignments, so the Ne proxy can be shown to be (or
not to be) sensitive to the estimation model.

| method | what it is | why include it |
|---|---|---|
| **M0** | one ω tree-wide (per-gene pooled, and on the concatenate) | the constrained baseline — no lineage variation allowed |
| **2-ratio, per-gene pooled** | BUSCOmega's method: `model=2`, focal species foreground, count-pooled across genes | the number the pipeline reports |
| **free-ratio, per-gene pooled** | `model=1` on each gene (every branch its own ω), pool the focal terminal branch across genes | tests whether the 2-ratio constraint (one background ω) matters |
| **free-ratio on the concatenate** | `model=1` once on all genes joined into one alignment; read the terminal ω directly | the Galtier-lab / eLife approach (Bourguignon et al. 2024; Weyna et al. preprint), transcribed into codeml |

If the four agree to ~2 decimals, the Ne proxy is robust to the model
choice and the PI's "why not free-ratio" question is closed. The 3-taxon
pilot is the **weakest** case for free-ratio and 2-ratio to differ (only 3
terminal branches, no internal branch); the comparison is more informative
on the 23-taxon Chapter 3 set, where it should be re-run.

```bash
conda activate buscomega
python compare_omega_methods.py \
    --aln-dir  ../examples/pilot_stage3_out/codon_aln \
    --tree     ../examples/pilot_stage4_out/trees/m0.nwk \
    --stage7-dir ../examples/pilot_stage7_out \
    -o out_compare/
```

Output: `method_comparison.tsv` (method × species).
