#!/usr/bin/env python3
"""
Introduce missing data into a simulated phylogenetic sequence alignment
to match the gap pattern observed in a reference (empirical) alignment.

The program preserves per-site missingness structure: for each alignment
column, the same set of taxa that have missing data in the reference
alignment will be masked in the simulated alignment.

Usage:
    python introduce_gaps.py reference.fasta simulated.fasta -o output.fasta

Inputs:
    reference.fasta  - Empirical alignment with realistic missing data pattern
    simulated.fasta  - Simulated alignment (no gaps) with the same taxa

Output:
    A modified version of the simulated alignment where missing data has been
    introduced to match the per-site pattern of the reference alignment.

Supported formats: FASTA (.fasta, .fa, .fas) and PHYLIP (.phy, .phylip)

Notes:
    - Missing data is represented as '-' or '?' characters in the reference.
    - Taxa must match between the two alignments (by name).
    - If alignments differ in length, the program can operate in two modes:
      (1) direct column mapping (truncate/pad) or (2) proportional mapping.
"""

import argparse
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
                    # Remove the earlier occurrence from taxa_order
                    taxa_order = [t for t in taxa_order if t != current_taxon]
                taxa_order.append(current_taxon)
                current_seq = []
            else:
                current_seq.append(line)
        if current_taxon is not None:
            sequences[current_taxon] = "".join(current_seq)

    return sequences, taxa_order


def parse_phylip(filepath):
    """Parse a PHYLIP format alignment file (sequential or interleaved).

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

    # Determine if this is sequential or interleaved based on nchar.
    # In sequential format, each taxon's full sequence follows its name
    # (possibly on multiple lines). In interleaved, all taxa appear in
    # blocks of equal length.
    first_seq_len = len(sequences[taxa_order[0]]) if taxa_order else 0

    if first_seq_len >= nchar:
        # Sequential format: first block already has complete sequences.
        # Truncate to nchar in case of trailing whitespace issues.
        for taxon in taxa_order:
            sequences[taxon] = sequences[taxon][:nchar]
    elif idx < len(lines):
        # Need more data. Determine if sequential-multiline or interleaved.
        # Heuristic: if after the first block we have exactly ntaxa taxa
        # and the next lines don't look like taxon names (no match to known
        # taxa), treat as sequential-multiline. Otherwise, interleaved.
        #
        # Standard approach: try interleaved (most common multi-block format).
        # In interleaved, subsequent blocks have sequence data only (no names).
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
    """Write sequences in PHYLIP format.

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
            # PHYLIP taxon names are typically padded to 10 characters
            f.write(f"{taxon:<10} {sequences[taxon]}\n")


def extract_gap_mask(sequences, taxa_order):
    """Extract the binary gap mask from an alignment.

    For each site (column), determine which taxa have missing data.
    Missing data characters are: '-', '?', 'X', 'x'

    Parameters
    ----------
    sequences : dict
        Dictionary mapping taxon names to sequences.
    taxa_order : list
        Ordered list of taxon names.

    Returns
    -------
    list of sets
        For each alignment column, a set of taxon names that have missing data.
    int
        Alignment length (number of columns).
    """
    if not taxa_order:
        return [], 0

    alignment_length = len(sequences[taxa_order[0]])
    gap_chars = set("-?Xx")

    gap_mask = []
    for col in range(alignment_length):
        missing_taxa = set()
        for taxon in taxa_order:
            if col < len(sequences[taxon]) and sequences[taxon][col] in gap_chars:
                missing_taxa.add(taxon)
        gap_mask.append(missing_taxa)

    return gap_mask, alignment_length


def apply_gap_mask(sequences, taxa_order, gap_mask, gap_char="-"):
    """Apply a gap mask to an alignment, introducing missing data.

    Parameters
    ----------
    sequences : dict
        Dictionary mapping taxon names to sequences (simulated, no gaps).
    taxa_order : list
        Ordered list of taxon names.
    gap_mask : list of sets
        For each column, a set of taxon names that should be masked.
    gap_char : str
        Character to use for missing data (default: '-').

    Returns
    -------
    dict
        Modified sequences with gaps introduced.
    """
    sim_length = len(sequences[taxa_order[0]])
    mask_length = len(gap_mask)

    # Build modified sequences
    modified = {}
    for taxon in taxa_order:
        seq_list = list(sequences[taxon])
        for col in range(min(sim_length, mask_length)):
            if taxon in gap_mask[col]:
                seq_list[col] = gap_char
        modified[taxon] = "".join(seq_list)

    return modified


def apply_gap_mask_proportional(sequences, taxa_order, gap_mask, gap_char="-"):
    """Apply a gap mask using proportional column mapping.

    When the reference and simulated alignments differ in length,
    map columns proportionally so the overall pattern of missingness
    is preserved across the full length of the simulated alignment.

    Parameters
    ----------
    sequences : dict
        Dictionary mapping taxon names to sequences.
    taxa_order : list
        Ordered list of taxon names.
    gap_mask : list of sets
        Gap mask from the reference alignment.
    gap_char : str
        Character to use for missing data.

    Returns
    -------
    dict
        Modified sequences with gaps introduced.
    """
    sim_length = len(sequences[taxa_order[0]])
    mask_length = len(gap_mask)

    # If the mask is empty, return sequences unmodified
    if mask_length == 0 or sim_length == 0:
        return {taxon: sequences[taxon] for taxon in taxa_order}

    modified = {}
    for taxon in taxa_order:
        seq_list = list(sequences[taxon])
        for col in range(sim_length):
            # Map simulated column to reference column proportionally
            ref_col = int(col * mask_length / sim_length)
            ref_col = min(ref_col, mask_length - 1)
            if taxon in gap_mask[ref_col]:
                seq_list[col] = gap_char
        modified[taxon] = "".join(seq_list)

    return modified


def summarize_missingness(sequences, taxa_order, label="", gap_char="-"):
    """Print a summary of missing data in an alignment.

    Parameters
    ----------
    sequences : dict
        Dictionary mapping taxon names to sequences.
    taxa_order : list
        Ordered list of taxon names.
    label : str
        Label for the summary output.
    gap_char : str
        The gap character used for missing data.
    """
    if not taxa_order:
        return

    alignment_length = len(sequences[taxa_order[0]])
    gap_chars = set("-?Xx")
    gap_chars.add(gap_char)

    print(f"\n{'=' * 60}")
    print(f"Missingness Summary: {label}")
    print(f"{'=' * 60}")
    print(f"  Taxa: {len(taxa_order)}")
    print(f"  Alignment length: {alignment_length}")

    # Per-taxon missingness
    print(f"\n  Per-taxon missingness:")
    total_missing = 0
    for taxon in taxa_order:
        n_missing = sum(1 for c in sequences[taxon] if c in gap_chars)
        total_missing += n_missing
        pct = 100 * n_missing / alignment_length if alignment_length > 0 else 0
        print(f"    {taxon:<20} {n_missing:>6} / {alignment_length} ({pct:.1f}%)")

    # Per-site missingness distribution
    site_missingness = []
    for col in range(alignment_length):
        n_missing = sum(
            1 for taxon in taxa_order if sequences[taxon][col] in gap_chars
        )
        site_missingness.append(n_missing)

    if site_missingness:
        avg_site = sum(site_missingness) / len(site_missingness)
        max_site = max(site_missingness)
        sites_with_gaps = sum(1 for x in site_missingness if x > 0)
        print(f"\n  Per-site missingness:")
        print(f"    Average taxa missing per site: {avg_site:.2f}")
        print(f"    Maximum taxa missing at a site: {max_site}")
        print(
            f"    Sites with any missing data: {sites_with_gaps} / "
            f"{alignment_length} ({100 * sites_with_gaps / alignment_length:.1f}%)"
        )

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
            "Introduce missing data into a simulated phylogenetic alignment "
            "to match the gap pattern of a reference alignment."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - apply gap pattern from empirical to simulated data
  python introduce_gaps.py empirical.fasta simulated.fasta -o output.fasta

  # Use proportional mapping when alignments differ in length
  python introduce_gaps.py empirical.fasta simulated.fasta -o output.fasta --proportional

  # Show missingness summary statistics
  python introduce_gaps.py empirical.fasta simulated.fasta -o output.fasta --summary

  # Use '?' as the gap character instead of '-'
  python introduce_gaps.py empirical.fasta simulated.fasta -o output.fasta --gap-char '?'
""",
    )

    parser.add_argument(
        "reference",
        help="Reference alignment with realistic gap pattern (FASTA or PHYLIP)",
    )
    parser.add_argument(
        "simulated",
        help="Simulated alignment to introduce gaps into (FASTA or PHYLIP)",
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
        help="Output format (default: same as simulated input)",
    )
    parser.add_argument(
        "--gap-char",
        default="-",
        help="Character to use for missing data (default: '-')",
    )
    parser.add_argument(
        "--proportional",
        action="store_true",
        help=(
            "Use proportional column mapping when alignments differ in length. "
            "By default, columns are mapped 1:1 (direct mapping)."
        ),
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help="Print missingness summary statistics",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Require all taxa in reference to be present in simulated alignment",
    )

    args = parser.parse_args()

    # Validate input files exist
    if not Path(args.reference).is_file():
        print(
            f"Error: Reference file not found: '{args.reference}'", file=sys.stderr
        )
        sys.exit(1)
    if not Path(args.simulated).is_file():
        print(
            f"Error: Simulated alignment file not found: '{args.simulated}'",
            file=sys.stderr,
        )
        sys.exit(1)

    # Parse input alignments
    ref_seqs, ref_taxa = parse_alignment(args.reference)
    sim_seqs, sim_taxa = parse_alignment(args.simulated)

    if not ref_seqs:
        print("Error: Reference alignment is empty.", file=sys.stderr)
        sys.exit(1)
    if not sim_seqs:
        print("Error: Simulated alignment is empty.", file=sys.stderr)
        sys.exit(1)

    # Validate that all sequences within each alignment have the same length
    ref_lengths = {taxon: len(seq) for taxon, seq in ref_seqs.items()}
    if len(set(ref_lengths.values())) > 1:
        print(
            "Error: Reference alignment has sequences of different lengths "
            "(not a valid alignment).",
            file=sys.stderr,
        )
        for taxon, length in ref_lengths.items():
            if length != len(ref_seqs[ref_taxa[0]]):
                print(
                    f"  {taxon}: {length} (expected {len(ref_seqs[ref_taxa[0]])})",
                    file=sys.stderr,
                )
        sys.exit(1)

    sim_lengths = {taxon: len(seq) for taxon, seq in sim_seqs.items()}
    if len(set(sim_lengths.values())) > 1:
        print(
            "Error: Simulated alignment has sequences of different lengths "
            "(not a valid alignment).",
            file=sys.stderr,
        )
        for taxon, length in sim_lengths.items():
            if length != len(sim_seqs[sim_taxa[0]]):
                print(
                    f"  {taxon}: {length} (expected {len(sim_seqs[sim_taxa[0]])})",
                    file=sys.stderr,
                )
        sys.exit(1)

    # Find common taxa
    ref_taxa_set = set(ref_taxa)
    sim_taxa_set = set(sim_taxa)
    common_taxa = ref_taxa_set & sim_taxa_set

    if not common_taxa:
        print(
            "Error: No matching taxa found between reference and simulated alignments.",
            file=sys.stderr,
        )
        print(f"  Reference taxa (first 5): {ref_taxa[:5]}", file=sys.stderr)
        print(f"  Simulated taxa (first 5): {sim_taxa[:5]}", file=sys.stderr)
        sys.exit(1)

    only_in_ref = ref_taxa_set - sim_taxa_set
    only_in_sim = sim_taxa_set - ref_taxa_set

    if only_in_ref:
        msg = (
            f"Warning: {len(only_in_ref)} taxa in reference not found in simulated: "
            f"{sorted(only_in_ref)[:5]}"
        )
        if args.strict:
            print(f"Error: {msg}", file=sys.stderr)
            sys.exit(1)
        else:
            print(msg, file=sys.stderr)

    if only_in_sim:
        print(
            f"Warning: {len(only_in_sim)} taxa in simulated not found in reference "
            f"(will not have gaps introduced): {sorted(only_in_sim)[:5]}",
            file=sys.stderr,
        )

    # Extract gap mask from reference
    gap_mask, ref_length = extract_gap_mask(ref_seqs, ref_taxa)
    sim_length = len(sim_seqs[sim_taxa[0]])

    print(f"Reference alignment: {len(ref_taxa)} taxa, {ref_length} sites")
    print(f"Simulated alignment: {len(sim_taxa)} taxa, {sim_length} sites")
    print(f"Common taxa: {len(common_taxa)}")

    if ref_length != sim_length:
        if args.proportional:
            print(
                f"Alignment lengths differ ({ref_length} vs {sim_length}). "
                f"Using proportional column mapping."
            )
        else:
            print(
                f"Alignment lengths differ ({ref_length} vs {sim_length}). "
                f"Using direct column mapping (columns beyond shorter alignment "
                f"will be unaffected)."
            )

    # Apply gap mask
    if args.proportional:
        modified_seqs = apply_gap_mask_proportional(
            sim_seqs, sim_taxa, gap_mask, args.gap_char
        )
    else:
        modified_seqs = apply_gap_mask(sim_seqs, sim_taxa, gap_mask, args.gap_char)

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
        write_phylip(modified_seqs, sim_taxa, args.output)
    else:
        write_fasta(modified_seqs, sim_taxa, args.output)

    print(f"\nOutput written to: {args.output}")

    # Print summaries if requested
    if args.summary:
        summarize_missingness(ref_seqs, ref_taxa, "Reference (empirical)", args.gap_char)
        summarize_missingness(
            modified_seqs, sim_taxa, "Output (simulated + gaps)", args.gap_char
        )


if __name__ == "__main__":
    main()
