#!/usr/bin/env python3
"""
Introduce missing data into a simulated phylogenetic sequence alignment
with shared (grouped) gap patterns among specified individuals.

Individuals can be assigned to groups so that members of the same group
share similar patterns of missing data (i.e., the same sites are masked).
Two modes control how similar the patterns are within a group:

  - identical mode (--overlap 1.0): All group members receive gaps at
    exactly the same sites (100% overlap).
  - partial overlap mode (--overlap 0.0-0.99): A shared base set of gap
    sites is selected for the group, then each individual's pattern is
    varied by swapping some gap sites for non-gap sites and vice versa.
    The --overlap value controls what fraction of each individual's gaps
    overlap with the shared group pattern.

Usage:
    python introduce_gaps_grouped.py simulated.fasta groups.tsv -o output.fasta

Inputs:
    simulated alignment - A simulated alignment with no missing data
                          (FASTA or relaxed PHYLIP format)
    groups file         - A tab-delimited file with columns:
                          taxon_name<TAB>group_name<TAB>percent_missing
                          (percent_missing is a value between 0 and 100)

Output:
    A modified version of the simulated alignment where missing data has been
    introduced with shared patterns within specified groups.

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

    first_seq_len = len(sequences[taxa_order[0]]) if taxa_order else 0

    if first_seq_len >= nchar:
        for taxon in taxa_order:
            sequences[taxon] = sequences[taxon][:nchar]
    elif idx < len(lines):
        while idx < len(lines):
            for taxon in taxa_order:
                if idx < len(lines):
                    sequences[taxon] += lines[idx].replace(" ", "")
                    idx += 1
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
        return parse_fasta(filepath)


def parse_groups(filepath):
    """Parse a tab-delimited file specifying group membership and missing data percentages.

    The file should have three columns separated by tabs:
        taxon_name<TAB>group_name<TAB>percent_missing

    Lines starting with '#' are treated as comments and skipped.
    A header line is skipped if the third column is not numeric.

    Parameters
    ----------
    filepath : str or Path
        Path to the tab-delimited groups file.

    Returns
    -------
    dict
        Dictionary mapping group names to lists of (taxon_name, percent_missing) tuples.
    """
    groups = {}

    with open(filepath, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                print(
                    f"Warning: Line {line_num} in groups file does not have "
                    f"three tab-separated columns, skipping: '{line}'",
                    file=sys.stderr,
                )
                continue
            taxon = parts[0].strip()
            group = parts[1].strip()
            try:
                pct = float(parts[2].strip())
            except ValueError:
                # Skip header lines where third column is not numeric
                if line_num == 1:
                    continue
                print(
                    f"Warning: Line {line_num} has non-numeric percentage "
                    f"'{parts[2].strip()}', skipping.",
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

            if group not in groups:
                groups[group] = []
            groups[group].append((taxon, pct))

    return groups


def introduce_grouped_gaps(
    sequences, taxa_order, groups, overlap=1.0, gap_char="-", seed=None
):
    """Introduce missing data with shared patterns within groups.

    For each group, a shared set of candidate gap sites is selected. The
    overlap parameter controls how much each individual's final gap pattern
    matches the shared group pattern.

    Parameters
    ----------
    sequences : dict
        Dictionary mapping taxon names to sequences (no missing data).
    taxa_order : list
        Ordered list of taxon names.
    groups : dict
        Dictionary mapping group names to lists of (taxon_name, percent_missing).
    overlap : float
        Fraction of gap sites shared among group members (0.0 to 1.0).
        1.0 means identical gap patterns within the group.
    gap_char : str
        Character to use for missing data (default: '-').
    seed : int or None
        Random seed for reproducibility (default: None).

    Returns
    -------
    dict
        Modified sequences with grouped gaps introduced.
    """
    if seed is not None:
        random.seed(seed)

    seq_len = len(sequences[taxa_order[0]])
    all_positions = list(range(seq_len))

    # Start with original sequences
    modified = {taxon: list(sequences[taxon]) for taxon in taxa_order}

    # Track which taxa have been processed via groups
    processed_taxa = set()

    for group_name, members in groups.items():
        # Filter to members that exist in the alignment
        valid_members = [
            (taxon, pct) for taxon, pct in members if taxon in sequences
        ]
        if not valid_members:
            continue

        # Determine the maximum number of gap sites needed in this group
        # (used to select the shared pool of candidate sites)
        max_n_gaps = max(
            int(round(pct / 100.0 * seq_len)) for _, pct in valid_members
        )
        max_n_gaps = min(max_n_gaps, seq_len)

        if max_n_gaps == 0:
            for taxon, _ in valid_members:
                processed_taxa.add(taxon)
            continue

        # Number of shared sites based on overlap fraction
        # shared_n is the number of sites ALL members will have in common
        # For each member, their gaps = shared_sites + individual_sites
        n_shared = int(round(overlap * max_n_gaps))
        n_shared = min(n_shared, seq_len)

        # Select the shared gap sites for this group
        shared_positions = set(random.sample(all_positions, n_shared))

        # Remaining positions available for individual variation
        remaining_positions = [p for p in all_positions if p not in shared_positions]

        for taxon, pct in valid_members:
            n_gaps = int(round(pct / 100.0 * seq_len))
            n_gaps = min(n_gaps, seq_len)

            if n_gaps == 0:
                processed_taxa.add(taxon)
                continue

            if overlap >= 1.0:
                # 100% overlap: all members use exactly the shared positions
                # truncated/extended to match their individual percentage
                if n_gaps <= len(shared_positions):
                    # Use a subset of shared positions
                    taxon_positions = set(random.sample(
                        sorted(shared_positions), n_gaps
                    ))
                else:
                    # Use all shared + draw extra from remaining
                    extra_needed = n_gaps - len(shared_positions)
                    extra_needed = min(extra_needed, len(remaining_positions))
                    extra = set(random.sample(remaining_positions, extra_needed))
                    taxon_positions = shared_positions | extra
            else:
                # Partial overlap: some sites are shared, rest are individual
                # n_from_shared = how many of this taxon's gaps come from shared pool
                n_from_shared = int(round(overlap * n_gaps))
                n_from_shared = min(n_from_shared, len(shared_positions), n_gaps)
                n_individual = n_gaps - n_from_shared
                n_individual = min(n_individual, len(remaining_positions))

                # Select from shared pool
                from_shared = set(random.sample(
                    sorted(shared_positions), n_from_shared
                ))
                # Select individual-specific positions
                from_individual = set(random.sample(
                    remaining_positions, n_individual
                ))

                taxon_positions = from_shared | from_individual

            # Apply gaps
            for pos in taxon_positions:
                modified[taxon][pos] = gap_char

            processed_taxa.add(taxon)

    # Convert lists back to strings
    result = {}
    for taxon in taxa_order:
        result[taxon] = "".join(modified[taxon])

    return result


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


def summarize_group_overlap(sequences, taxa_order, groups, gap_char="-"):
    """Print overlap statistics for each group.

    Parameters
    ----------
    sequences : dict
        Dictionary mapping taxon names to sequences.
    taxa_order : list
        Ordered list of taxon names.
    groups : dict
        Dictionary mapping group names to lists of (taxon_name, percent_missing).
    gap_char : str
        The gap character used for missing data.
    """
    gap_chars = set("-?Xx")
    gap_chars.add(gap_char)
    alignment_length = len(sequences[taxa_order[0]])

    print(f"\n{'=' * 60}")
    print("Group Overlap Summary")
    print(f"{'=' * 60}")

    for group_name, members in groups.items():
        valid_members = [
            taxon for taxon, _ in members if taxon in sequences
        ]
        if len(valid_members) < 2:
            continue

        print(f"\n  Group: {group_name} ({len(valid_members)} members)")

        # Get gap positions for each member
        gap_sets = {}
        for taxon in valid_members:
            gap_sets[taxon] = set(
                i for i in range(alignment_length)
                if sequences[taxon][i] in gap_chars
            )

        # Compute pairwise overlap (Jaccard similarity)
        overlaps = []
        for i in range(len(valid_members)):
            for j in range(i + 1, len(valid_members)):
                t1, t2 = valid_members[i], valid_members[j]
                s1, s2 = gap_sets[t1], gap_sets[t2]
                if s1 or s2:
                    union_size = len(s1 | s2)
                    intersect_size = len(s1 & s2)
                    jaccard = intersect_size / union_size if union_size > 0 else 0
                    # Also compute overlap coefficient (intersection / min)
                    min_size = min(len(s1), len(s2))
                    overlap_coef = (
                        intersect_size / min_size if min_size > 0 else 0
                    )
                    overlaps.append((t1, t2, jaccard, overlap_coef, intersect_size))

        if overlaps:
            avg_jaccard = sum(o[2] for o in overlaps) / len(overlaps)
            avg_overlap = sum(o[3] for o in overlaps) / len(overlaps)
            print(f"    Average Jaccard similarity: {avg_jaccard:.3f}")
            print(f"    Average overlap coefficient: {avg_overlap:.3f}")
            for t1, t2, jac, ovl, n_shared in overlaps:
                print(
                    f"    {t1} vs {t2}: "
                    f"shared={n_shared}, Jaccard={jac:.3f}, overlap={ovl:.3f}"
                )

    print(f"{'=' * 60}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Introduce missing data into a simulated phylogenetic alignment "
            "with shared (grouped) gap patterns among specified individuals."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 100%% shared gaps within groups (default)
  python introduce_gaps_grouped.py simulated.fasta groups.tsv -o output.fasta

  # 80%% overlap within groups (some individual variation)
  python introduce_gaps_grouped.py simulated.fasta groups.tsv -o output.fasta --overlap 0.8

  # Fully independent gaps (no group sharing, equivalent to introduce_gaps_random)
  python introduce_gaps_grouped.py simulated.fasta groups.tsv -o output.fasta --overlap 0.0

  # Use a specific random seed for reproducibility
  python introduce_gaps_grouped.py simulated.fasta groups.tsv -o output.fasta --seed 42

  # Show overlap statistics among group members
  python introduce_gaps_grouped.py simulated.fasta groups.tsv -o output.fasta --summary

Groups file format (tab-delimited):
  taxon1\tgroupA\t25.0
  taxon2\tgroupA\t30.0
  taxon3\tgroupB\t10.0
  taxon4\tgroupB\t15.0

The --overlap parameter controls how similar gap patterns are within a group:
  1.0 = identical sites are masked for all group members (100%% overlap)
  0.8 = 80%% of each member's gap sites come from a shared pool
  0.0 = fully independent random gaps (no shared pattern)
""",
    )

    parser.add_argument(
        "alignment",
        help="Simulated alignment with no missing data (FASTA or relaxed PHYLIP)",
    )
    parser.add_argument(
        "groups",
        help=(
            "Tab-delimited file with three columns: "
            "taxon_name, group_name, and percent_missing (0-100)"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output file path for the gapped alignment",
    )
    parser.add_argument(
        "--overlap",
        type=float,
        default=1.0,
        help=(
            "Fraction of gap sites shared among group members (0.0 to 1.0). "
            "1.0 = identical patterns (100%% overlap), "
            "0.0 = fully independent random patterns. (default: 1.0)"
        ),
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
        help="Print missingness and group overlap summary statistics",
    )

    args = parser.parse_args()

    # Validate overlap parameter
    if args.overlap < 0.0 or args.overlap > 1.0:
        print(
            f"Error: --overlap must be between 0.0 and 1.0, got {args.overlap}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Validate input files exist
    if not Path(args.alignment).is_file():
        print(
            f"Error: Alignment file not found: '{args.alignment}'", file=sys.stderr
        )
        sys.exit(1)
    if not Path(args.groups).is_file():
        print(
            f"Error: Groups file not found: '{args.groups}'", file=sys.stderr
        )
        sys.exit(1)

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

    # Parse groups file
    groups = parse_groups(args.groups)

    if not groups:
        print("Error: No valid group entries found in file.", file=sys.stderr)
        sys.exit(1)

    # Check for taxa in groups file that are not in the alignment
    alignment_taxa_set = set(taxa_order)
    group_taxa = set()
    for members in groups.values():
        for taxon, _ in members:
            group_taxa.add(taxon)

    missing_from_alignment = group_taxa - alignment_taxa_set
    missing_from_groups = alignment_taxa_set - group_taxa

    if missing_from_alignment:
        print(
            f"Warning: {len(missing_from_alignment)} taxa in groups file not "
            f"found in alignment: {sorted(missing_from_alignment)[:5]}",
            file=sys.stderr,
        )

    if missing_from_groups:
        print(
            f"Note: {len(missing_from_groups)} taxa in alignment not assigned to "
            f"any group (no gaps will be introduced for these): "
            f"{sorted(missing_from_groups)[:5]}",
            file=sys.stderr,
        )

    alignment_length = len(sequences[taxa_order[0]])
    print(f"Alignment: {len(taxa_order)} taxa, {alignment_length} sites")
    print(f"Groups: {len(groups)} ({', '.join(groups.keys())})")
    print(f"Taxa with group assignments: {len(group_taxa & alignment_taxa_set)}")
    print(f"Overlap: {args.overlap:.1%}")

    # Introduce grouped gaps
    modified_seqs = introduce_grouped_gaps(
        sequences, taxa_order, groups, args.overlap, args.gap_char, args.seed
    )

    # Determine output format
    output_path = Path(args.output)
    if args.output_format:
        out_format = args.output_format
    else:
        out_ext = output_path.suffix.lower()
        if out_ext in (".phy", ".phylip"):
            out_format = "phylip"
        elif out_ext in (".fasta", ".fa", ".fas", ".fna", ".faa"):
            out_format = "fasta"
        else:
            in_ext = Path(args.alignment).suffix.lower()
            if in_ext in (".phy", ".phylip"):
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
        summarize_missingness(
            modified_seqs, taxa_order, "Output (with grouped gaps)", args.gap_char
        )
        summarize_group_overlap(modified_seqs, taxa_order, groups, args.gap_char)


if __name__ == "__main__":
    main()
