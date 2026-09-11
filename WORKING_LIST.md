# BUSCO → dN/dS → Ne proxy: Working List

Status tracker for building this as a standalone, reusable pipeline (eventual
own git repo, README on GitHub). Nothing here is built yet except what's
explicitly marked done — this is the plan, not a log of finished work.

Legend: `[ ]` not started · `[~]` in progress · `[x]` done · `[?]` open decision, blocking

---

## Open decisions (resolve before building past this point)

- `[~]` **codeml branch model**: leaning **2-ratio branch-labeled**
  (`model=2`, one run per focal species per gene) as the tool's single
  engine — more general, more robust as taxa grow, degrades gracefully.
  Free-ratio (`model=1`) *not* built as a co-equal option unless a real
  user asks. Docs carry the decision rule (few taxa / no focal lineage →
  free-ratio territory; many taxa or a focal lineage → 2-ratio). For the
  Chapter 3 run specifically: 2-ratio, one run per ingroup species.
  **Confirm before Stage 4.**
- `[x]` **Repo location**: consolidated. `busco-dnds-ne\` folded into
  `working/` (its parser → `codes/core/06_parse_codeml_output.py`, its test
  and example data → `tests/` and `examples/`) and removed. Everything for
  this sub-project now lives under `ne_pipeline/working/`; spin out to its
  own git repo once it actually works end to end.
- `[x]` **`cleandata` setting**: standardized to `1` (strip ambiguous/gapped
  columns) across every codeml run, so stages stay comparable. Don't force
  `-nogap` **or `-nomismatch`** at PAL2NAL — let codeml's `cleandata` remove
  gaps once, consistently; a real translation mismatch drops the gene
  (Stage 3) rather than being silently trimmed. Full rationale in
  `docs/methods.md` Stage 3 ("How codeml handles gaps"). Caveat to check
  before Phase 2: at 23 taxa, `cleandata = 1` complete-deletion may discard
  many sites — measure site retention, fall back to `cleandata = 0` if poor.
- `[x]` **Dependencies**: Python side stays **stdlib-only**, matching
  `protein-preprocessing-isoform-pipeline`. No Biopython, no pandas, no
  gffread, no seqkit. Conda env is only for the external binaries.
- `[x]` **Sequence input**: consume the cleaned, ID-matched protein+CDS
  output of `protein-preprocessing-isoform-pipeline` (LINGUA Stage 1)
  rather than re-solving protein↔CDS matching in `busco_seq_extractor.py`.
- `[x]` **Parallelization model** (documented in `docs/methods.md`
  "Conventions"): `--jobs` (concurrent workers, default 1) everywhere;
  `--threads` (default 4) at Stage 3 only — MAFFT is threaded, codeml is
  not; `--batch-size` (genes per codeml call via `ndata`, default 1) at
  Stage 5 only. Batching is result-neutral (codeml resets between datasets;
  parser splits on `Dataset N`). Each concurrent codeml worker runs in its
  own dir (codeml scatters scratch files into cwd).
  Mechanism: stdlib `concurrent.futures.ProcessPoolExecutor`.

## Confirmed so far

- `[x]` **Original pipeline scripts recovered** from `E:\codes\` (not lost —
  were never in `E:\br_dng\`), then reorganized by necessity into
  `ne_pipeline/working/codes/`:
  `core/` — `01_common_scos.py`, `02_extract_sco_seqs.py`,
  `03_codon_align.py`, `04_codeml_control.py`, `05_run_codeml.py`,
  `06_parse_codeml_output.py`;
  `legacy/` — `merge_busco_table.py`, `busco_seq_extractor.py`,
  `run_mafft.sh`, `run_pal2nal.sh`, `generate_codeml_control.py`,
  `run_codeml.sh` (all superseded);
  `prep_optional/` — `run_busco.py`, `busco_summary.py`, `species_tree.py`
  (only needed for setup — BUSCO run, species tree — not part of the core
  tool);
  `not_used_here/` — `sort_phytozome_files.sh` (Phytozome-download-specific,
  belongs to the broader personal toolkit, not this tool). Not copied at
  all (available on the drive if wanted later): the longest-isoform family
  (`longest_isoforms.py` and friends — only relevant if a genome source
  needs isoform selection) and `run_orthofinder.sh` (that's for the
  separate DNG_protocol orthogroup step, not the Ne proxy).
  **This changes Stages 1-5 below from "write from scratch" to "review and
  adapt existing code" — reread those sections.**
- `[x]` BUSCO lineage: **`brassicales_odb10`** (order-level; Brassicaceae
  itself is not offered as its own family-level set). **4,596 markers** (verified from
  the tarball scores_cutoff; earlier "~665" was wrong); the common single-copy set
  across a set of Brassicaceae runs a few thousand, shrinking as more distantly
  related taxa (outgroups) are added.
- `[x]` Pilot data available as real test input: 3 species (*A. halleri*,
  *A. thaliana*, *C. grandiflora*), 368 common BUSCO genes, protein + CDS
  sequences already pulled per gene.

---

## Two phases — don't conflate them

**Phase 1 — build & validate BUSCOmega** on the **3-species pilot**
(*A. halleri*, *A. thaliana*, *C. grandiflora*) using the
**`viridiplantae_odb10`** lineage — the existing pilot BUSCO output already
uses it (`at33090` in the gene IDs), so there is **no BUSCO re-run**. Small
lineage (~425 markers, ~300–400 common), so the whole build/test loop is
fast. This is the deliverable: a working, reusable, documented pipeline.
Every stage below can reach `[x]` on this pilot.

**Phase 2 — apply the finished tool to the Chapter 3 dataset** (21
Brassicaceae + 2 outgroups) using **`brassicales_odb10`** (4,596 markers).
One downstream run of the completed pipeline, plus the 23-taxon BUSCO run
and species tree. Blocked on Phase 1. Produces the Ne-proxy values for
H3 / H5.

Stages 1–7 = Phase 1, on 3 species + viridiplantae. Phase 2 is the last
section.

---

## Stage 1 — BUSCO output → common single-copy gene set  *(Phase 1, viridiplantae)*

- `[x]` `codes/core/01_common_scos.py` written (supersedes
  `merge_busco_table.py`, now in `codes/legacy/`): `Complete`-only filter
  (single-copy orthologs), strict intersection, emits `common_scos.tsv` +
  `busco_status_summary.tsv`. Stdlib. Author block added, version 0.1.0.
- `[x]` Test: `tests/test_common_scos.py` + synthetic fixture — passes,
  every filter branch exercised (Complete-in-all / Duplicated / Missing /
  Fragmented). Demo output produced and shown.
- `[x]` `docs/methods.md` Stage 1 section written
- `[x]` Ran on real pilot data (`examples/pilot_busco/`, 3× protein-mode
  `full_table.tsv` from `E:\br_dng`): **369 common single-copy genes**
  (old pipeline: 368). Outputs in `examples/pilot_busco_stage1_out/`.
- `[x]` Summary table now splits non-Complete calls into **sole-blocker**
  (recoverable) vs. shared loss, per status. Pilot `recoverable`:
  a_thaliana 3, a_halleri 30, c_grandiflora 21 — a_halleri is the weak link.
- `[x]` Format verified current: pilot ran BUSCO v5.1.2/odb10, checked
  against v6.1.0/odb12.2 — `full_table.tsv` layout unchanged.
- `[x]` Also grabbed for Stage 2: 3× per-species proteome + CDS
  (`examples/pilot_proteomes/{prt,cds}/`). Drive can disconnect.
- `[ ]` *(deferred)* `--rescue-fragments <pct>` — pull near-complete
  Fragmented hits into the common set using the `Length` column vs. the
  lineage `lengths_cutoff`. Tradeoff: shorter alignment → noisier ω for
  that gene. Per-run choice, not default.
- `[ ]` *(deferred)* status-summary chart — stacked horizontal bar
  (BUSCO convention), emitted as hand-written SVG to stay stdlib-only.

**Stage 1 complete for the pilot** (core + sole-blocker enrichment).

## Stage 2 — Sequence extraction  *(Phase 1, 3 spp)*

- `[x]` `codes/core/02_extract_sco_seqs.py` written (stdlib, author block).
  Old `busco_seq_extractor.py` → `codes/legacy/` (rewrite, not adaptation:
  fuzzy substring matching replaced with exact + suffix-strip).
- `[x]` Protein: exact ID lookup. CDS: exact → built-in suffix strips
  (`.p .P -p -P .cds .CDS`) + `--cds-suffix-strip` → give up (no fuzzy).
  Match rule logged per species.
- `[x]` Length check `|cds − 3·prt| ≤ slop+3` → `length_flag` in manifest,
  not a drop. Whole-gene drop only if a species is missing protein or CDS.
- `[x]` Outputs: `prt/<id>.faa`, `cds/<id>.fna` (all species per file),
  `stage2_manifest.tsv` (long), `stage2_dropped.tsv`. High drop/flag rate
  prints an isoform-cleaning note.
- `[x]` Stage 1 got the matching heuristic: high Duplicated fraction prints
  the isoform-cleaning note there too.
- `[x]` Ran on real pilot: **369/369 written, 0 dropped, 0 length flags.**
  `c_grandiflora` matched via `strip:.p`, others exact. Output in
  `examples/pilot_stage2_out/`.
- `[x]` Test: `tests/test_extract_sco_seqs.py` + fixture (exact / strip /
  drop / length-flag). `docs/methods.md` Stage 2 section written.

**Stage 2 complete for the pilot.**

## Stage 3 — Codon alignment  *(Phase 1, 3 spp)*

- `[x]` `codes/core/03_codon_align.py` — one stdlib script replacing
  `run_mafft.sh` + `run_pal2nal.sh` (both → `codes/legacy/`)
- `[x]` Per gene: strip trailing stop → `mafft --auto --inputorder
  --thread T` → `pal2nal.pl -output paml`
- `[x]` PAL2NAL run with **no `-nogap` and no `-nomismatch`**: keep gap
  columns for codeml to remove once (`cleandata = 1`); a translation
  mismatch drops + logs the gene instead of silently deleting codons
- `[x]` `--jobs` concurrent genes (ProcessPoolExecutor, own temp dir each),
  `--threads` per MAFFT
- `[x]` Outputs: `prot_aln/*.aln`, `codon_aln/*.pml`, `stage3_manifest.tsv`
  (n_taxa, aln lengths, `gap_frac`, tool exit codes), `stage3_dropped.tsv`
  (busco_id, step ∈ input/mafft/pal2nal, reason); >5% drop rate warns
- `[x]` `tests/test_codon_align.py` (binaries faked) + `paml_stats` /
  `strip_stop` unit tests — all pass. `docs/methods.md` Stage 3 section
  written, incl. the codeml gap-handling / `cleandata` explanation and the
  23-taxon site-retention caveat to check before Phase 2
- `[x]` Real-binary run: `examples/pilot_stage2_out/` (369 genes) in WSL
  `buscomega` env (MAFFT 7.526 + PAL2NAL 14.1), `--jobs 8` → **369/369
  aligned, 0 dropped, all exit codes 0**, ~54 s. `gap_frac` median 0.7% /
  mean 2.2% / max 43% (23 genes >10%, 4 >20%). Output +manifests kept in
  `examples/pilot_stage3_out/` (~4 MB) for building Stage 4.

**Stage 3 complete for the pilot.**

## Stage 4 — codeml analysis set-up  *(Phase 1, 3 spp)*

- `[x]` `codes/core/04_codeml_control.py` — **rewrite** of
  `generate_codeml_control.py` (→ `codes/legacy/`), not an extension. Model
  decision resolved: **M0 (model=0) + 2-ratio (model=2), one 2-ratio run
  per focal species.** Not free-ratio (noisy), not branch-site (different
  question).
- `[x]` Template approach: one `.ctl` per analysis (`ctl/m0.ctl`,
  `ctl/two_ratio.ctl`) with `__SEQFILE__ __TREEFILE__ __OUTFILE__ __NDATA__`
  placeholders; `write_control()` is importable by Stage 5. Pilot = 4
  templates+trees, Phase 2 = 22. No per-gene `.ctl` files.
- `[x]` Labelled trees: `trees/m0.nwk` + `trees/two_ratio/<sp>.nwk` (that
  tip ` #1`, substring-safe). **Leading `<ntax> 1` header on every tree** —
  `ntax` from the tree, `ntrees=1` fixed on purpose (identical species set
  every gene). Without it, batched codeml misreads + exits nonzero.
- `[x]` **codeml exit code is unreliable** (returns 1 on clean runs without
  the tree header, and other cases) → Stage 5 must gate on a parseable
  `.mlc`, not `$?`. Documented.
- `[x]` `cleandata=1` baked in (override warns); `CodonFreq`/`kappa`/`omega`/
  `icode` as flags; choices recorded in `stage4_params.tsv`.
- `[x]` `--stage3-dir` cross-check: tree tips vs Stage 3 species → abort on
  mismatch (wrong topology biases every ω).
- `[x]` `tests/test_codeml_control.py` (5 test fns) + full suite green.
  End-to-end vs real codeml (`buscomega` env): `ndata=4` batches of pilot
  genes → exit 0, correct per-gene ω, M0 and 2-ratio. Output in
  `examples/pilot_stage4_out/`. `docs/methods.md` Stage 4 section written.
- `[x]` `prep_optional/species_tree.py` — `prepare` (clean any Newick:
  lengths/support/`[&…]`/NEXUS/quotes → topology; `--rename` reconciles tip
  names; `--match-to` diffs vs `common_scos.tsv`) and `infer` (concat Stage
  3 alignments → `supermatrix.phy` + `partitions.txt`; `--run-iqtree` runs
  the ML search — IQ-TREE 2, not in the env). Neither adds the `<ntax> 1`
  header — no tree file carries it; Stage 4 does. `tests/test_species_tree.py`
  + `species_tree.md`. Pilot: cleans the pilot tree, builds a 369-gene /
  187,766-aa supermatrix.
- `[x]` `docs/methods.md` — "Running BUSCOmega on your own clade" section
  (clade-agnostic; the 4 per-clade inputs).

**Stage 4 complete for the pilot.**

## Stage 5 — Run codeml  *(Phase 1, 3 spp)*

- `[x]` `codes/core/05_run_codeml.py` — Python rewrite of `run_codeml.sh`
  (→ `codes/legacy/`). Reads `stage4_analyses.tsv` + Stage 3 `codon_aln/`.
  Per (analysis × batch): concat `.pml`, copy the analysis tree, fill the 4
  placeholders (`fill_template()`), run codeml with `cwd` = the batch dir.
- `[x]` `--jobs` concurrent batches (ProcessPoolExecutor), `--batch-size`
  genes per call. **Success = one result block per gene in the `.mlc`, not
  the exit code.** Short batch → retry each gene solo into
  `retry/<id>.mlc`; still-failing → `stage5_failed.tsv`, rest kept.
- `[x]` `stage5_manifest.tsv` (analysis, batch, n_genes, n_ok, n_datasets,
  status, wall_s, mlc), `stage5_failed.tsv`, `stage5.log`.
- `[x]` **Real pilot run** (`buscomega`, PAML 4.10.10), 369 genes ×
  4 analyses, `--batch-size 20 --jobs 8`: **1476/1476 ok, 0 failed,
  ~470 s**. Per-analysis: M0 289 s, 2-ratio a_halleri 894 / a_thaliana
  827 / c_grandiflora 1583 s. Parsed end-to-end through Stage 6.
- `[x]` **Finding → feeds Stage 7:** 2-ratio foreground QC-flag rate scales
  with lineage divergence — a_halleri/a_thaliana ≈1%, **c_grandiflora
  ≈23%** (omega_boundary + ds_floor). 3 taxa + 1 short gene = little
  foreground signal. Stage 7 must dS-weight / exclude boundary cases, not
  average raw per-gene ω. Better powered at 23 taxa but pattern persists.
- Output kept: `examples/pilot_stage5_out/` (manifests + 4 parsed record
  tables + 1 sample batch dir, ~0.8 MB).

**Stage 5 complete for the pilot.**

## Stage 6 — Parse codeml output  *(Phase 1, 3 spp)*

- `[x]` `codes/core/06_parse_codeml_output.py` — stdlib, deterministic tip
  resolution (100%), QC flags (`omega_boundary` 0/999, `ds_floor` <1e-3).
- `[x]` `--stage5-dir` mode added: walks `<analysis>/batch_*/batch.mlc`,
  maps `Dataset N` → busco_id via the batch `genes.txt` (Stage 5 now writes
  it), picks up `retry/<id>.mlc`, → `06_parsed/<analysis>_records.tsv`.
- `[x]` Confirmed on real batched M0 + 2-ratio output (via the 5→6→7 chain).
- `[x]` Renamed 07→06 to match the stage number.

## Stage 7 — the Ne proxy  *(Phase 1, 3 spp)*

Design + justification: **`docs/primer.md` §8–9** (figures from real pilot
data). Pilot (5→6→7 chain): mean ω 5.6/3.0/226, median 0.14/0.15/0.18,
**count-pooled 0.163/0.170/0.151** (95% CIs mostly overlap). The
dS-ceiling filter caught genes with dS≈70 (mis-alignments) that otherwise
crushed the estimate to ~0.035.

- `[x]` `codes/core/07_ne_proxy.py` (stdlib). `--records-dir 06_parsed/`
  **or** `--stage5-dir 05_codeml_out/` (runs Stage 6 first).
- `[x]` Ne proxy per species = **count-pooled ω** on that species'
  foreground branch, `Σ(N·dN)/Σ(S·dS)` over genes. Not a flag; mean &
  median written alongside for transparency.
- `[x]` Filtering: keep `omega_boundary`; drop `ds_floor`; drop a gene for
  a lineage when its dS there > `--ds-ceiling` (default 1.5). Direction
  rule (deep clade → hold/lower, never raise) documented in methods.md +
  primer.
- `[x]` `--bootstrap` (default 1000, `--seed`): gene resample → 2.5/97.5 CI.
- `[x]` M0-only fallback: no `two_ratio.*` tables → emit just the
  genome-wide pooled M0 ω.
- `[x]` Outputs: `ne_proxy.tsv`, `per_gene_omega.tsv` (every gene, kept or
  dropped, with reason), `plots/{ne_proxy_forest,omega_per_gene_dist,
  omega_vs_divergence}.svg`, `stage7.log`.
- `[x]` `--png`: rasterise via `rsvg-convert` / `inkscape` / `cairosvg` if
  present; warn+skip otherwise. SVG is the source of truth. On for the
  dissertation runs.
- `[x]` `tests/test_ne_proxy.py` — pooling arithmetic, ds_floor/ceiling
  filtering + boundary retention, seeded bootstrap determinism, end-to-end.
- `[x]` Real pilot 5→6→7 chain run: 369 genes → per-species pooled ω
  0.163 / 0.170 / 0.151 (a_hal / a_tha / c_gra), CIs mostly overlapping.
  Output in `examples/pilot_stage7_out/`.
- `[x]` **Bug caught + fixed during the pilot run** (2026-09-10): pooled
  ω must be a ratio of pooled *rates*, `[Σ(N·dN)/ΣN] / [Σ(S·dS)/ΣS]`, not
  of pooled substitution counts `Σ(N·dN)/Σ(S·dS)` (which is ~3× high).
  **Verified correct now:** `test_ne_proxy.py::test_pooled_omega` asserts
  the exact formula on two hand-worked fixtures (equal and unequal gene
  sizes); the re-run chain gives pooled ω 0.163/0.170/0.151, sitting right
  on the per-gene medians (0.14/0.15/0.18) as expected for a well-behaved
  estimator. `primer.md` §8.4 (with the explicit "not `Σ(N·dN)/Σ(S·dS)`"
  note), §8.3/§9.3 tables, glossary, and all 4 figures updated to match.

## Scope & known gaps

- `[x]` "Scope and boundaries" section in `docs/methods.md` — what
  BUSCOmega is / is not (no site models, no free-ratio, no ancestral
  branches, no paralogs, no isoform selection — all by design, all
  documented).
- `[ ]` Stage 5 `--resume` — skip batches whose `.mlc` already parses.
  Matters for the multi-day 21+2 `brassicales_odb10` run.
- `[ ]` Stage 3 `--mafft-extra` passthrough for power users.
- `[ ]` Stage 1 `--rescue-fragments <pct>` (deferred earlier).
- `[x]` **`prep_optional/run_busco.py`** (2026-09-11) — rewrote the
  BUSCO-runner prerequisite in Python, picked from the two `run_busco.sh`
  variants recovered in `codes.zip` (sequential loop; fixed-size `&`/`wait`
  parallel batches). Real work queue (`ProcessPoolExecutor`, `--jobs` x
  `--threads`, same convention as the core stages) instead of wait-batches
  that idle on the slowest job; per-species drop-and-continue on a BUSCO
  failure; `--resume`; a manifest (`species status elapsed_s complete_pct`);
  calls `busco_summary.py` at the end. Old `.sh` variants → `codes/legacy/`.
  `tests/test_run_busco.py` (busco faked). `run_busco.md` rewritten.
- Note: the H5 regression (ω vs per-species DNG count) is Chapter 3 stats,
  not BUSCOmega — Stage 7 just emits the per-species table.

## Audit (2026-09-10, full pass)

- `[x]` **stdlib-only verified** — AST scan of all 7 core + 3 prep scripts,
  zero third-party imports.
- `[x]` **Determinism verified** — sorted inputs everywhere, seeded
  bootstrap, codeml deterministic, `results.sort()` before writing.
- `[x]` **Fixed:** `06_parse_codeml_output.py` `main()` now returns an int
  and is wrapped in `SystemExit` (was `main()` → always exit 0).
- `[x]` **Fixed:** `03_codon_align.py` `align_one` now has a top-level
  `except` → a malformed alignment drops the gene (`step=error`) instead of
  crashing the stage (Stage 5 already had this).
- `[x]` **Fixed:** `07_ne_proxy.py` plots guard against NaN / empty inputs
  (no M0 records, or all genes filtered) → a placeholder SVG, not a crash.
- `[x]` **Packaging blocker resolved:** `examples/pilot_proteomes/` (151 MB)
  is gitignored, confirmed not tracked (`git ls-files` → 0). Largest
  tracked blob on GitHub is 168 KB. Nowhere near a size problem.
- `[ ]` Minor (not blocking): Stage 1 `busco_status_summary.tsv`
  `duplicated` / `total` count *rows* not *genes* (Duplicated = 2+
  rows/gene); Stage 2 `n_flagged` can exceed the manifest's flag=1 count
  when a gene is flagged then later dropped. Cosmetic; document or tidy.
- `[x]` Pipeline illustration: `docs/figures/pipeline.svg` (regenerable via
  `_generate.py`), referenced from `methods.md`.

---

## Output structure (documented in `docs/methods.md`)

- `[x]` Decided: one run dir, numbered stage subdirs (`01_common_scos/` …
  `06_ne_proxy/`), matching the isoform pipeline convention.
- `[x]` Gene tracking: **per-stage manifests + a regenerable `audit` join**,
  not a mutated master table. Stage 2 already emits its manifests.
- `[x]` `run_summary.tsv` at the run root = concise per-stage counts + params
  + tool versions; the "start here if something's wrong" file.
- `[x]` **Top-level driver: `run_buscomega.py`** — runs stages 1→7 into
  the numbered subdirs as subprocesses (nothing re-implemented), stops on
  first failure, captures per-stage stdout. `--from-stage`/`--to-stage`
  slice, `--resume` (→ Stage 5), `--dry-run`. Tested end-to-end on the
  3-species pilot.
- `[x]` **`audit` join** (`build_audit` in the driver): `gene_tracking.tsv`
  = one row per gene, its fate per stage. Regenerated from manifests, not
  mutated. `run_summary.tsv` = every `stageN.log` line in one file.
- `[x]` `tests/test_orchestrator.py` — audit join + summary on a synthetic
  run dir (no subprocess/binaries).

## Packaging (Phase 1)

- `[x]` `docs/primer.md` — nearly-neutral theory from first principles, the
  estimator (§8–9) with figures, glossary.
- `[x]` `docs/methods.md` — every step + parameter + reasoning + citation;
  orchestrator + output structure + scope sections.
- `[x]` `README.md` — what it is, install, quick-start one-liner, I/O,
  the estimator in short, scope, tests, citation.
- `[x]` `environment.yml` (+ optional add-ons block), `docs/install.md`
  (conda / no-conda / WSL).
- `[x]` Stage 5 binary check + Stage 3/07 tool checks already report
  missing binaries; driver prints tools-on-PATH at start.
- `[x]` `LICENSE` (MIT), `CITATION.cff`, `pyproject.toml` (metadata +
  pytest config), `.gitignore` (excludes the 151 MB proteomes),
  `.github/workflows/tests.yml` (8 suites + stdlib-only check on 3.9/3.12).
- `[x]` `git init` + first commit + **pushed to
  github.com/Ludtson/BUSCOmega** (2026-09-11). 13 commits, `main` tracked.
- `[ ]` `--check-deps` top-level mode (report binaries + versions) — nice
  to have; the per-stage checks already cover the failure case.
- `[x]` Fancier pipeline figure: funnel layout + distinct glyph per stage,
  phylo trees hand-finished in Inkscape. `docs/figures/pipeline.svg`.
- `[ ]` Confirm CI (`tests.yml`) is green on GitHub now that it's public —
  check the Actions tab.

---

## Phase 2 — apply BUSCOmega to the Chapter 3 dataset (21 + 2)

- `[ ]` BUSCO run on all 23 taxa (`brassicales_odb10`)
- `[ ]` 23-taxon species tree (reuse Chapter 2's, or build)
- `[ ]` Run the finished pipeline end to end
- `[ ]` Output → the Ne-proxy table that feeds H3 / H5

**Blocked on Phase 1 being done.**
