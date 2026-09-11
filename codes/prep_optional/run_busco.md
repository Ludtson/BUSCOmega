# run_busco.py — run BUSCO across species (prerequisite)

**Optional.** BUSCOmega's Stage 1 starts from BUSCO's own `full_table.tsv`
per species. If you already have BUSCO output, skip this entirely. If you
don't, this runs it for you — in parallel, with a manifest and a resume
mode.

Rewritten in Python from the original toolkit's `run_busco.sh`, now in
`codes/legacy/`. The reasoning is the same as for the 7 core stages: a real
work queue instead of fixed-size `&`/`wait` batches, a manifest, one
species' failure doesn't kill the run, and it's testable.

## Dependencies

`busco` (>=5.5) on `PATH`. **Install it into its own conda env, not
`buscomega`:**

```bash
conda create -n busco -c bioconda -c conda-forge python=3.11 'busco>=5.5'
```

This script is stdlib-only, so it runs fine under any env's Python —
activate `busco` for this one step, then switch back:

```bash
conda activate busco
python run_busco.py <fasta_dir> -o <out_dir> --lineage <lineage_dir> --mode protein
conda activate buscomega   # back to this for the rest of the pipeline
```

**Why a separate env, not `-c conda-forge busco>=5.5` on top of
`buscomega`:** two different, real solve failures doing that —

1. `-c bioconda` alone: augustus (a BUSCO dependency) needs conda-forge's
   `boost-cpp`/`gsl`/`lp_solve` → "augustus ... no viable options".
2. `-c bioconda -c conda-forge` together, into an *existing* env: conda
   tries not to change packages already installed, which acts as an
   implicit pin on that env's Python version. If `buscomega`'s Python
   ended up newer than what BUSCO's `sepp` dependency has builds for
   (e.g. Python 3.14, released after the newest `sepp` build only goes up
   to <3.14), the solve fails outright — there is no BUSCO/sepp release
   for that Python yet.

Pinning `python=3.11` in a **fresh** env sidesteps both: nothing is
"already installed" to protect, and 3.11 has compatible builds all the
way through BUSCO's dependency tree. It also keeps BUSCO's large
dependency pull (hmmer, metaeuk, augustus, sepp, pandas, matplotlib, ...)
out of the env that runs mafft/pal2nal/codeml.

You also need a **pre-downloaded lineage dataset** — this script runs
BUSCO `--offline`:

```bash
curl -L https://busco-data.ezlab.org/v5/data/lineages/viridiplantae_odb10.2024-01-08.tar.gz \
    -o viridiplantae_odb10.tar.gz
tar -xzf viridiplantae_odb10.tar.gz
# -> pass the extracted directory as --lineage
```

## Usage

```bash
python run_busco.py <fasta_dir> -o <out_dir> \
    --lineage <lineage_dir> --mode protein \
    --jobs 4 --threads 6
```

| argument | meaning |
|---|---|
| `fasta_dir` | one `<species>.{faa,fa,fna,fasta}` per species. Species name = filename up to the first dot. |
| `--lineage` | the extracted lineage dataset directory (see above) |
| `--mode` | `genome` \| `protein` \| `transcriptome`. Use **`protein`** — that's what BUSCOmega's input contract expects (isoform-cleaned proteomes; see `docs/methods.md` "Input contract"). |
| `--jobs` | species run concurrently (default 1) |
| `--threads` | `--cpu` passed to each BUSCO call (default 4). Total cores used = `jobs x threads` — same convention as Stage 3. |
| `--resume` | skip a species that already has a `full_table.tsv` under `<out_dir>/<species>/` |

## What it does

Per species, in a worker pool of size `--jobs`:

```
busco --offline --in <fasta> --cpu <threads> \
      --lineage_dataset <lineage> --mode <mode> \
      --out <species> --out_path <out_dir> --force
```

captures stdout+stderr to `<out_dir>/logs/<species>.log`, and checks for a
resulting `full_table.tsv` to decide `ok` vs `failed`. A failed species is
logged and the run continues — it does not stop the other species.

## Outputs (into `--out-dir`)

- `<species>/run_<lineage>/full_table.tsv` — BUSCO's own output; this is
  what Stage 1 (`01_common_scos.py`) reads
- `logs/<species>.log` — that species' full BUSCO stdout+stderr
- `busco_run_manifest.tsv` — `species status elapsed_s complete_pct`
- `busco_summary.csv` — via `busco_summary.py`, run automatically at the
  end if at least one species succeeded
- `busco_run.log` — one timestamped start/done line (same `stageN.log`
  shape as the core stages)

## Test

`tests/test_run_busco.py` — species discovery, `--resume` skipping,
per-species drop-and-continue on a BUSCO failure, `complete_pct` parsing,
and the manifest — with `busco` replaced by a fake (the binary is not
assumed present in CI).
