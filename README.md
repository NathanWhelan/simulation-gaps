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

### Configuration File

The configuration file (`params.cfg`) uses a simple INI format. See `examples/params.cfg` for a fully-commented template.

Key settings:
- `alignment` — Your empirical alignment file (PHYLIP or FASTA)
- `ratio` — Signal conflict ratio (e.g., `70:30`)
- `tree1` / `tree2` — Two competing tree topologies
- `models` — Substitution models for each partition (e.g., `WAG+C10, 1.0, 1`)
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

## License

MIT