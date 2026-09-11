# Installing BUSCOmega

BUSCOmega's Python code is **standard-library only** — nothing to `pip install`.
All it needs is a Python interpreter and a few external bioinformatics
binaries on `PATH`.

## Option A — conda (recommended)

```bash
conda env create -f environment.yml     # or: mamba env create -f environment.yml
conda activate buscomega
```

This creates a named environment **`buscomega`** with:

| binary | stage | why |
|---|---|---|
| `python` >=3.9 | all | runs the pipeline scripts (stdlib only) |
| `mafft` >=7.5 | 3 | protein alignment |
| `pal2nal` =14.1 | 3 | protein alignment -> codon alignment |
| `codeml` (from `paml` >=4.10) | 4/5 | the dN/dS model fits |

To update after pulling a new BUSCOmega version:

```bash
conda env update -f environment.yml --prune
```

### Optional add-ons

None of these are installed by default — each covers one non-core feature.
Add what you need: `conda install -n buscomega -c <channel> <pkg>`.

| package | channel | needed for |
|---|---|---|
| `busco` >=5.5 | **its own env**, not `buscomega` | the prerequisite BUSCO run (`prep_optional/run_busco.py`) — only if you have no BUSCO output yet. `conda create -n busco -c bioconda -c conda-forge python=3.11 busco>=5.5`. Adding it to the existing `buscomega` env can fail to solve two ways: augustus needs conda-forge (bioconda alone → "no viable options"), and even with both channels, conda's reluctance to change an already-installed Python version can collide with BUSCO's `sepp` dependency having no build for a newer Python. A fresh env with `python=3.11` avoids both — see `run_busco.md` for the full story. |
| `iqtree` >=2.2 | bioconda | inferring a species tree from your alignments (`prep_optional/species_tree.py infer --run-iqtree`). Not needed if you supply a tree. |
| `librsvg` | conda-forge | PNG copies of the Stage 7 plots (`07_ne_proxy.py --png`) — provides `rsvg-convert`. `inkscape` or `cairosvg` also work. SVG is always written regardless. |

### On Windows: run it in WSL, not Git Bash

MAFFT, PAL2NAL, and codeml are Unix tools; bioconda has no Windows builds.
On a Windows host, do everything inside WSL2:

```bash
wsl                                   # drop into your Linux distro
conda env create -f /mnt/c/.../ne_pipeline/working/environment.yml
conda activate buscomega
```

The project on your `C:` drive is visible from WSL at
`/mnt/c/Users/<you>/...`. WSL2 defaults to about half the host's cores and
RAM; if a large run needs more, create `C:\Users\<you>\.wslconfig`:

```ini
[wsl2]
memory=28GB
# processors defaults to every logical core WSL sees; set it lower to cap
```

then `wsl --shutdown` and reopen. (`processors` cannot exceed the logical
core count WSL reports — check with `nproc` inside WSL.)

This is how the pilot was run and validated (WSL2 Ubuntu, `buscomega` env,
MAFFT 7.526 + PAL2NAL 14.1).

## Option B — no conda (HPC modules, restricted environments, or preference)

Put these on `PATH` however your system does it (`module load`, a package
manager, source build):

| tool | min version | source |
|---|---|---|
| Python | 3.9 | any — stdlib only |
| MAFFT (`mafft`) | 7.5 | <https://mafft.cbrc.jp/alignment/software/> |
| PAL2NAL (`pal2nal.pl`) | 14.1 | <http://www.bork.embl.de/pal2nal/> (Perl script; needs `perl`) |
| PAML (`codeml`) | 4.10 | <http://abacus.gene.ucl.ac.uk/software/paml.html> |
| BUSCO (`busco`) | 5.5 | <https://busco.ezlab.org/> — prerequisite step only |

Every stage script checks for the binaries it needs at start-up and prints
where to get anything missing, so a partial install fails loudly rather than
mid-run.

## Version notes

- The codeml output parser (`06_parse_codeml_output.py`) was validated
  against PAML **4.10.7** output. Re-run the parser test if you move to a
  different minor version.
- BUSCO's `full_table.tsv` format is stable across v5 and v6; the lineage
  dataset (`odb10` vs `odb12`) only changes gene-ID suffixes, not the
  columns BUSCOmega reads.
