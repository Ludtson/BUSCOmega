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

`busco` (>=5.5) on `PATH`. Deliberately **not** in the core `buscomega`
conda env — most users already have BUSCO output, and it pulls a large
dependency tree (hmmer, metaeuk, augustus/miniprot). Add it if you need it:

```bash
conda install -n buscomega -c bioconda 'busco>=5.5'
```

or keep BUSCO in its own environment and just put it on `PATH` when you
run this script.

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
