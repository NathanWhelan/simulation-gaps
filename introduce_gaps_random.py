#!/usr/bin/env python3
"""
Introduce random missing data into a simulated phylogenetic sequence alignment
on a per-individual basis, using percentages specified in a tab-delimited file.

The missing data is inserted at random positions independently for each
individual, so the pattern for one individual is unrelated to any other.

Usage:
    python introduce_gaps_random.py simulated.phy missing_percentages.tsv -o output.phy

Inputs:
    simulated alignment - A simulated alignment with no missing data
                          (relaxed PHYLIP or FASTA format)
    percentages file    - A tab-delimited file with two columns:
                          taxon_name<TAB>percent_missing
                          (percent_missing is a value between 0 and 100)

Output:
    A modified version of the simulated alignment where missing data has been
    introduced randomly at the specified percentage for each individual.

Supported formats: FASTA (.fasta, .fa, .fas) and relaxed PHYLIP (.phy, .phylip)
"""

import argparse
import random
import sys
from pathlib import Path


def parse_fasta(filepath):
    """Parse a FASTA format alignment file.

    Parameters
    ----------
    filepath : str or Path
        Path to the FASTA file.

    Returns
    -------
    dict
        Dictionary mapping taxon names to sequences (as strings).
    list
        Ordered list of taxon names (preserves input order).
    """
    sequences = {}
    taxa_order = []
    current_taxon = None
    current_seq = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_taxon is not None:
                    sequences[current_taxon] = "".join(current_seq)
                current_taxon = line[1:].split()[0]
                if current_taxon in sequences:
                    print(
                        f"Warning: Duplicate taxon '{current_taxon}' in FASTA file. "
                        f"Only the last occurrence will be used.",
                        file=sys.stderr,
                    )
                    taxa_order = [t for t in taxa_order if t != current_taxon]
                taxa_order.append(current_taxon)
                current_seq = []
            else:
                current_seq.append(line)
        if current_taxon is not None:
            sequences[current_taxon] = "".join(current_seq)

    return sequences, taxa_order


def parse_phylip(filepath):
    """Parse a relaxed, non-interleaved (sequential) PHYLIP format alignment file.

    Parameters
    ----------
    filepath : str or Path
        Path to the PHYLIP file.

    Returns
    -------
    dict
        Dictionary mapping taxon names to sequences (as strings).
    list
        Ordered list of taxon names.
    """
    sequences = {}
    taxa_order = []

    with open(filepath, "r") as f:
        lines = [line.rstrip() for line in f if line.strip()]

    header = lines[0].split()
    ntaxa = int(header[0])
    nchar = int(header[1])

    # Parse the first block: each line has "taxon_name sequence_data"
    idx = 1
    for i in range(ntaxa):
        if idx >= len(lines):
            print(
                f"Error: Expected {ntaxa} taxa but file ended after {i}.",
                file=sys.stderr,
            )
            break
        parts = lines[idx].split(None, 1)
        taxon = parts[0]
        seq = parts[1].replace(" ", "") if len(parts) > 1 else ""
        taxa_order.append(taxon)
        sequences[taxon] = seq
        idx += 1

    # Check if sequences are complete or if more lines are needed
    first_seq_len = len(sequences[taxa_order[0]]) if taxa_order else 0

    if first_seq_len >= nchar:
        # Sequential format: first block already has complete sequences.
        for taxon in taxa_order:
            sequences[taxon] = sequences[taxon][:nchar]
    elif idx < len(lines):
        # Interleaved format: subsequent blocks have sequence data only
        while idx < len(lines):
            for taxon in taxa_order:
                if idx < len(lines):
                    sequences[taxon] += lines[idx].replace(" ", "")
                    idx += 1

        # Truncate to nchar
        for taxon in taxa_order:
            sequences[taxon] = sequences[taxon][:nchar]

    return sequences, taxa_order


def parse_alignment(filepath):
    """Parse an alignment file, auto-detecting format from extension.

    Parameters
    ----------
    filepath : str or Path
        Path to the alignment file.

    Returns
    -------
    dict
        Dictionary mapping taxon names to sequences.
    list
        Ordered list of taxon names.
    """
    filepath = Path(filepath)
    ext = filepath.suffix.lower()

    if ext in (".fasta", ".fa", ".fas", ".fna", ".faa"):
        return parse_fasta(filepath)
    elif ext in (".phy", ".phylip"):
        return parse_phylip(filepath)
    else:
        # Default to FASTA
        return parse_fasta(filepath)


def parse_percentages(filepath):
    """Parse a tab-delimited file specifying per-individual missing data percentages.

    The file should have two columns separated by a tab:
        taxon_name<TAB>percent_missing

    Lines starting with '#' are treated as comments and skipped.
    A header line is skipped if the second column is not numeric.

    Parameters
    ----------
    filepath : str or Path
        Path to the tab-delimited percentages file.

    Returns
    -------
    dict
        Dictionary mapping taxon names to their missing data percentage (0-100).
    """
    percentages = {}

    with open(filepath, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                print(
                    f"Warning: Line {line_num} in percentages file does not have "
                    f"two tab-separated columns, skipping: '{line}'",
                    file=sys.stderr,
                )
                continue
            taxon = parts[0].strip()
            try:
                pct = float(parts[1].strip())
            except ValueError:
                # Skip header lines where second column is not numeric
                if line_num == 1:
                    continue
                print(
                    f"Warning: Line {line_num} has non-numeric percentage "
                    f"'{parts[1].strip()}', skipping.",
                    file=sys.stderr,
                )
                continue
            if pct < 0 or pct > 100:
                print(
                    f"Warning: Percentage for '{taxon}' is {pct}, "
                    f"which is outside [0, 100]. Clamping.",
                    file=sys.stderr,
                )
                pct = max(0, min(100, pct))
            percentages[taxon] = pct

    return percentages


def introduce_random_gaps(sequences, taxa_order, percentages, gap_char="-", seed=None):
    """Introduce random missing data into sequences based on per-taxon percentages.

    For each individual, randomly select positions to mask based on the
    specified percentage. The random selection is independent for each individual.

    Parameters
    ----------
    sequences : dict
        Dictionary mapping taxon names to sequences (no missing data).
    taxa_order : list
        Ordered list of taxon names.
    percentages : dict
        Dictionary mapping taxon names to their missing data percentage (0-100).
    gap_char : str
        Character to use for missing data (default: '-').
    seed : int or None
        Random seed for reproducibility (default: None).

    Returns
    -------
    dict
        Modified sequences with random gaps introduced.
    """
    if seed is not None:
        random.seed(seed)

    modified = {}
    alignment_length = len(sequences[taxa_order[0]]) if taxa_order else 0

    for taxon in taxa_order:
        seq_list = list(sequences[taxon])
        seq_len = len(seq_list)

        if taxon in percentages:
            pct = percentages[taxon]
            n_gaps = int(round(pct / 100.0 * seq_len))
            # Randomly select positions to mask
            if n_gaps > 0:
                positions = random.sample(range(seq_len), min(n_gaps, seq_len))
                for pos in positions:
                    seq_list[pos] = gap_char

        modified[taxon] = "".join(seq_list)

    return modified


def write_fasta(sequences, taxa_order, filepath, line_width=80):
    """Write sequences in FASTA format.

    Parameters
    ----------
    sequences : dict
        Dictionary mapping taxon names to sequences.
    taxa_order : list
        Ordered list of taxon names.
    filepath : str or Path
        Output file path.
    line_width : int
        Number of characters per sequence line.
    """
    with open(filepath, "w") as f:
        for taxon in taxa_order:
            f.write(f">{taxon}\n")
            seq = sequences[taxon]
            for i in range(0, len(seq), line_width):
                f.write(seq[i : i + line_width] + "\n")


def write_phylip(sequences, taxa_order, filepath):
    """Write sequences in relaxed PHYLIP format (non-interleaved).

    Parameters
    ----------
    sequences : dict
        Dictionary mapping taxon names to sequences.
    taxa_order : list
        Ordered list of taxon names.
    filepath : str or Path
        Output file path.
    """
    ntaxa = len(taxa_order)
    nchar = len(sequences[taxa_order[0]])

    with open(filepath, "w") as f:
        f.write(f"{ntaxa} {nchar}\n")
        for taxon in taxa_order:
            f.write(f"{taxon} {sequences[taxon]}\n")


def summarize_missingness(sequences, taxa_order, label=""):
    """Print a summary of missing data in an alignment.

    Parameters
    ----------
    sequences : dict
        Dictionary mapping taxon names to sequences.
    taxa_order : list
        Ordered list of taxon names.
    label : str
        Label for the summary output.
    """
    if not taxa_order:
        return

    alignment_length = len(sequences[taxa_order[0]])
    gap_chars = set("-?Xx")

    print(f"\n{'=' * 60}")
    print(f"Missingness Summary: {label}")
    print(f"{'=' * 60}")
    print(f"  Taxa: {len(taxa_order)}")
    print(f"  Alignment length: {alignment_length}")

    print(f"\n  Per-taxon missingness:")
    total_missing = 0
    for taxon in taxa_order:
        n_missing = sum(1 for c in sequences[taxon] if c in gap_chars)
        total_missing += n_missing
        pct = 100 * n_missing / alignment_length if alignment_length > 0 else 0
        print(f"    {taxon:<20} {n_missing:>6} / {alignment_length} ({pct:.1f}%)")

    overall_pct = (
        100 * total_missing / (len(taxa_order) * alignment_length)
        if alignment_length > 0
        else 0
    )
    print(f"\n  Overall missingness: {overall_pct:.1f}%")
    print(f"{'=' * 60}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Introduce random missing data into a simulated phylogenetic alignment "
            "based on per-individual percentages specified in a tab-delimited file."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python introduce_gaps_random.py simulated.phy percentages.tsv -o output.phy

  # Use a specific random seed for reproducibility
  python introduce_gaps_random.py simulated.fasta percentages.tsv -o output.fasta --seed 42

  # Show missingness summary statistics
  python introduce_gaps_random.py simulated.phy percentages.tsv -o output.phy --summary

  # Use '?' as the gap character instead of '-'
  python introduce_gaps_random.py simulated.phy percentages.tsv -o output.phy --gap-char '?'

Percentages file format (tab-delimited):
  taxon1\t25.0
  taxon2\t10.5
  taxon3\t50.0
""",
    )

    parser.add_argument(
        "alignment",
        help="Simulated alignment with no missing data (FASTA or relaxed PHYLIP)",
    )
    parser.add_argument(
        "percentages",
        help=(
            "Tab-delimited file with two columns: taxon_name and percent_missing "
            "(value between 0 and 100)"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output file path for the gapped alignment",
    )
    parser.add_argument(
        "--output-format",
        choices=["fasta", "phylip"],
        default=None,
        help="Output format (default: same as input alignment)",
    )
    parser.add_argument(
        "--gap-char",
        default="-",
        help="Character to use for missing data (default: '-')",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility (default: None)",
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help="Print missingness summary statistics",
    )

    args = parser.parse_args()

    # Parse input alignment
    sequences, taxa_order = parse_alignment(args.alignment)

    if not sequences:
        print("Error: Alignment is empty.", file=sys.stderr)
        sys.exit(1)

    # Validate that all sequences have the same length
    seq_lengths = {taxon: len(seq) for taxon, seq in sequences.items()}
    if len(set(seq_lengths.values())) > 1:
        print(
            "Error: Alignment has sequences of different lengths "
            "(not a valid alignment).",
            file=sys.stderr,
        )
        sys.exit(1)

    # Parse percentages file
    percentages = parse_percentages(args.percentages)

    if not percentages:
        print("Error: No valid percentages found in file.", file=sys.stderr)
        sys.exit(1)

    # Check for taxa in percentages file that are not in the alignment
    alignment_taxa_set = set(taxa_order)
    pct_taxa_set = set(percentages.keys())

    missing_from_alignment = pct_taxa_set - alignment_taxa_set
    missing_from_pct = alignment_taxa_set - pct_taxa_set

    if missing_from_alignment:
        print(
            f"Warning: {len(missing_from_alignment)} taxa in percentages file not "
            f"found in alignment: {sorted(missing_from_alignment)[:5]}",
            file=sys.stderr,
        )

    if missing_from_pct:
        print(
            f"Warning: {len(missing_from_pct)} taxa in alignment not found in "
            f"percentages file (no gaps will be introduced for these): "
            f"{sorted(missing_from_pct)[:5]}",
            file=sys.stderr,
        )

    alignment_length = len(sequences[taxa_order[0]])
    print(f"Alignment: {len(taxa_order)} taxa, {alignment_length} sites")
    print(f"Percentages specified for: {len(pct_taxa_set & alignment_taxa_set)} taxa")

    # Introduce random gaps
    modified_seqs = introduce_random_gaps(
        sequences, taxa_order, percentages, args.gap_char, args.seed
    )

    # Determine output format
    output_path = Path(args.output)
    if args.output_format:
        out_format = args.output_format
    else:
        out_ext = output_path.suffix.lower()
        if out_ext in (".phy", ".phylip"):
            out_format = "phylip"
        else:
            out_format = "fasta"

    # Write output
    if out_format == "phylip":
        write_phylip(modified_seqs, taxa_order, args.output)
    else:
        write_fasta(modified_seqs, taxa_order, args.output)

    print(f"\nOutput written to: {args.output}")

    # Print summaries if requested
    if args.summary:
        summarize_missingness(modified_seqs, taxa_order, "Output (with random gaps)")


if __name__ == "__main__":
    main()
