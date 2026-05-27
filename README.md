# simulation-gaps

Tools for simulating phylogenetic data with realistic missing data patterns and signal conflict using IQ-TREE's alisim.

## Overview

This repository provides:

1. **`simulate_conflict.py`** — An automated pipeline for simulating phylogenetic data with conflicting signal using alisim (IQ-TREE). Reads a simple configuration file and handles all steps automatically.
2. **`introduce_gaps.py`** — Transfer per-site gap patterns from an empirical alignment to a simulated alignment.
3. **`introduce_gaps_random.py`** — Introduce random missing data on a per-individual basis using specified percentages.
4. **`five-parts-list.py`** — Utility to divide an alignment into equal partitions.

## Quick Start: Signal-Conflict Simulation Pipeline

The main script `simulate_conflict.py` automates the entire signal-conflict simulation workflow.

### Requirements

- Python 3.7+ (no external Python packages needed)
- [IQ-TREE 3](https://github.com/iqtree/iqtree3) (with alisim support)
- [AMAS](https://github.com/marekborowiec/AMAS) (`pip install amas`)

### Usage

```bash
# 1. Copy and edit the example configuration file
cp examples/params.cfg my_simulation.cfg
# Edit my_simulation.cfg with your alignment, trees, models, etc.

# 2. Preview what the script will do (recommended first!)
python simulate_conflict.py my_simulation.cfg --dry-run

# 3. Run the full pipeline
python simulate_conflict.py my_simulation.cfg

# 4. For SLURM clusters, set use_slurm = yes in the config file
#    The script will generate SLURM batch scripts instead of running directly
python simulate_conflict.py my_simulation.cfg
```

### What It Does

The pipeline automates these steps:

1. **Splits** the input alignment into two portions based on a ratio (e.g., 70:30)
2. **Sub-divides** each portion into multiple partitions (default: 5)
3. **Simulates** data on each partition using alisim with different evolutionary models and specified tree topologies
4. **Introduces gaps** from the empirical alignment into the simulated data
5. **Concatenates** all simulated partitions into a final alignment

### No-Alignment Mode (with Indel Model)

If you don't have a starting alignment, you can specify `alignment_length` instead and optionally provide an indel model. In this mode:

- AMAS.py is not required (no alignment splitting)
- alisim uses `--length` to set the root sequence length per partition
- Gaps are introduced naturally by the indel model during simulation
- The `--site-freq SAMPLING` and `--site-rate SAMPLING` flags are not used

```bash
# Use the indel example config
python simulate_conflict.py examples/params_indel_no_alignment.cfg --dry-run
```

See `examples/params_indel_no_alignment.cfg` for a fully-commented example.

### Configuration File

The configuration file (`params.cfg`) uses a simple INI format. See `examples/params.cfg` for a fully-commented template.

Key settings:
- `alignment` — Your empirical alignment file (PHYLIP or FASTA), or leave blank for no-alignment mode
- `alignment_length` — Total number of sites to simulate (required when `alignment` is blank)
- `ratio` — Signal conflict ratio (e.g., `70:30`)
- `tree1` / `tree2` — Two competing tree topologies
- `models` — Substitution models for each partition (e.g., `WAG+C10, 1.0, 1`)
- `indel` — Indel model rates (e.g., `0.03,0.10` for insertion/deletion rates)
- `indel_size` — Indel size distribution (e.g., `POW{1.7}`, `GEO{0.5}`)
- `use_slurm` — Generate SLURM scripts for HPC clusters

### Command-Line Options

| Flag | Description |
|------|-------------|
| `--dry-run` | Preview all commands without running anything |
| `--verbose` | Show detailed progress and commands |
| `--post-simulation` | Run only gap-introduction and concatenation (after SLURM jobs finish) |

## Motivation

When simulating phylogenetic data (e.g., amino acid alignments), the simulated sequences typically contain no missing data. However, real alignments often have substantial missing data with complex patterns — certain sites may be missing for many taxa while others are complete. Simply introducing random gaps doesn't capture this structure.

This tool solves that problem by transferring the **per-site missingness pattern** from a real (reference) alignment onto a simulated alignment with the same taxa. This preserves:

- Which taxa are missing at each alignment column
- The correlation of missingness across taxa at the same position
- The overall distribution of gaps across the alignment

## Requirements

- Python 3.7+
- No external Python packages (uses only the standard library)
- For the full pipeline: IQ-TREE 3 and AMAS

## Installation

Clone this repository:

```bash
git clone https://github.com/NathanWhelan/simulation-gaps.git
cd simulation-gaps
```

No additional installation is needed.

## Usage

```bash
python introduce_gaps.py reference.fasta simulated.fasta -o output.fasta
```

### Arguments

| Argument | Description |
|----------|-------------|
| `reference` | Empirical alignment with realistic gap pattern (FASTA or PHYLIP) |
| `simulated` | Simulated alignment to introduce gaps into (FASTA or PHYLIP) |
| `-o`, `--output` | Output file path (required) |
| `--output-format` | Force output format: `fasta` or `phylip` (default: auto-detect from extension) |
| `--gap-char` | Character to use for gaps (default: `-`) |
| `--proportional` | Use proportional column mapping when alignments differ in length |
| `--summary` | Print detailed missingness statistics |
| `--strict` | Exit with error if any reference taxon is missing from simulated alignment |

### Examples

Basic usage — apply gap pattern from an empirical alignment to simulated data:

```bash
python introduce_gaps.py empirical.fasta simulated.fasta -o simulated_with_gaps.fasta
```

Show missingness summary statistics:

```bash
python introduce_gaps.py empirical.fasta simulated.fasta -o output.fasta --summary
```

Use proportional mapping when alignments differ in length:

```bash
python introduce_gaps.py empirical.fasta simulated.fasta -o output.fasta --proportional
```

Use `?` instead of `-` for the gap character:

```bash
python introduce_gaps.py empirical.fasta simulated.fasta -o output.fasta --gap-char '?'
```

## How It Works

1. **Parse** the reference (empirical) alignment and the simulated alignment.
2. **Extract the gap mask**: For each column in the reference alignment, determine which taxa have missing data (characters: `-`, `?`, `X`, `x`).
3. **Match taxa**: Find taxa that are common between the reference and simulated alignments.
4. **Apply the gap mask**: For each column in the simulated alignment, replace the amino acid with the gap character for all taxa that were missing at that position in the reference.
5. **Write** the modified alignment.

### Column Mapping Modes

- **Direct mapping** (default): Column *i* in the simulated alignment gets the gap mask from column *i* in the reference. If alignments differ in length, extra columns in the longer alignment are unaffected.
- **Proportional mapping** (`--proportional`): Columns are mapped proportionally, so that the spatial distribution of gaps is stretched/compressed to fit the simulated alignment length.

## Supported Formats

- **FASTA** (`.fasta`, `.fa`, `.fas`, `.fna`, `.faa`)
- **PHYLIP** (`.phy`, `.phylip`) — sequential and interleaved

## Input Requirements

- Taxa names must match between the reference and simulated alignments.
- The reference alignment should contain realistic missing data (gaps).
- The simulated alignment typically has no gaps (complete sequences).
- Both alignments should be aligned (all sequences same length within each file).

## Example Workflow

1. **Obtain empirical data**: Download a real amino acid alignment for your taxa of interest.
2. **Simulate data**: Use a phylogenetic simulator (e.g., Seq-Gen, INDELible, pyvolve) to generate complete sequences for the same taxa.
3. **Introduce gaps**: Use this tool to apply the empirical gap pattern to the simulated data.

```bash
# Step 3: Apply the empirical gap pattern
python introduce_gaps.py real_alignment.fasta simulated_alignment.fasta \
    -o simulated_with_realistic_gaps.fasta --summary
```

## Random Per-Individual Missing Data

A second script, `introduce_gaps_random.py`, introduces missing data **randomly** on a per-individual basis according to percentages specified in a tab-delimited file. Unlike the original script (which copies an empirical gap pattern), this script inserts gaps at random positions independently for each individual, so the missing data pattern for one individual is unrelated to any other.

### Usage

```bash
python introduce_gaps_random.py simulated.phy percentages.tsv -o output.phy
```

### Percentages File Format

A tab-delimited file with two columns (no header required):

```
taxon1	25.0
taxon2	10.5
taxon3	50.0
```

Each line specifies the taxon name and the percentage of sites (0–100) to replace with missing data.

### Arguments

| Argument | Description |
|----------|-------------|
| `alignment` | Simulated alignment with no missing data (FASTA or relaxed PHYLIP) |
| `percentages` | Tab-delimited file: taxon_name, percent_missing |
| `-o`, `--output` | Output file path (required) |
| `--output-format` | Force output format: `fasta` or `phylip` (default: auto-detect from extension) |
| `--gap-char` | Character to use for gaps (default: `-`) |
| `--seed` | Random seed for reproducibility |
| `--summary` | Print detailed missingness statistics |

### Examples

```bash
# Basic usage with a relaxed PHYLIP file
python introduce_gaps_random.py simulated.phy percentages.tsv -o output.phy

# Use a random seed for reproducibility
python introduce_gaps_random.py simulated.fasta percentages.tsv -o output.fasta --seed 42

# Show missingness summary
python introduce_gaps_random.py simulated.phy percentages.tsv -o output.phy --summary
```

## Grouped Missing Data (Shared Patterns)

A third script, `introduce_gaps_grouped.py`, introduces missing data with **shared patterns** among specified groups of individuals. This is useful for simulating scenarios where certain taxa (e.g., closely related species or samples from the same sequencing run) tend to be missing data at the same alignment positions.

### Key Features

- **100% overlap mode** (`--overlap 1.0`): All group members receive gaps at exactly the same sites.
- **Partial overlap mode** (`--overlap 0.0–0.99`): A shared base set of gap sites is selected for the group, then each individual's pattern is varied. The overlap value controls what fraction of each individual's gaps come from the shared pool.
- **Group overlap statistics**: Use `--summary` to see Jaccard similarity and overlap coefficients between group members.

### Usage

```bash
python introduce_gaps_grouped.py simulated.fasta groups.tsv -o output.fasta
```

### Groups File Format

A tab-delimited file with three columns (no header required):

```
taxon1	groupA	25.0
taxon2	groupA	30.0
taxon3	groupB	10.0
taxon4	groupB	15.0
```

Each line specifies: taxon name, group name, and percentage of sites (0–100) to replace with missing data.

### Arguments

| Argument | Description |
|----------|-------------|
| `alignment` | Simulated alignment with no missing data (FASTA or relaxed PHYLIP) |
| `groups` | Tab-delimited file: taxon_name, group_name, percent_missing |
| `-o`, `--output` | Output file path (required) |
| `--overlap` | Fraction of gap sites shared within groups (0.0–1.0, default: 1.0) |
| `--output-format` | Force output format: `fasta` or `phylip` (default: auto-detect) |
| `--gap-char` | Character to use for gaps (default: `-`) |
| `--seed` | Random seed for reproducibility |
| `--summary` | Print missingness and group overlap statistics |

### Examples

```bash
# 100% shared gaps within groups (identical gap patterns)
python introduce_gaps_grouped.py simulated.fasta groups.tsv -o output.fasta --overlap 1.0

# 80% overlap within groups (some individual variation)
python introduce_gaps_grouped.py simulated.fasta groups.tsv -o output.fasta --overlap 0.8

# Fully independent gaps (no group sharing)
python introduce_gaps_grouped.py simulated.fasta groups.tsv -o output.fasta --overlap 0.0

# With reproducibility and summary statistics
python introduce_gaps_grouped.py simulated.fasta groups.tsv -o output.fasta --overlap 0.8 --seed 42 --summary
```

### How the Overlap Parameter Works

- `--overlap 1.0`: A shared set of gap positions is selected for the group. All members receive gaps at exactly the same sites (adjusted for their individual percentages).
- `--overlap 0.8`: 80% of each member's gap sites are drawn from a shared pool; the remaining 20% are selected independently for each individual.
- `--overlap 0.0`: All gap sites are selected independently (equivalent to `introduce_gaps_random.py` with per-individual randomness).

## License

MIT