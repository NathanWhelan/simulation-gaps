#!/usr/bin/env python3
"""
simulate_conflict.py - Automated signal-conflict simulation pipeline using alisim.

This script automates the process of simulating phylogenetic data with
conflicting phylogenetic signal using IQ-TREE's alisim. It reads simulation
parameters from a configuration file and performs all steps automatically.

Workflow:
    1. Split the input alignment into two portions based on a ratio
    2. Further split each portion into sub-partitions
    3. Simulate data on each sub-partition using alisim with specified models
    4. Optionally introduce empirical gap patterns into simulated data
    5. Concatenate all simulated partitions into a final alignment

Usage:
    python simulate_conflict.py params.cfg
    python simulate_conflict.py params.cfg --dry-run
    python simulate_conflict.py params.cfg --verbose

For biologists: This script replaces the manual steps of splitting alignments,
running alisim on each partition, introducing gaps, and concatenating results.
All you need to do is fill in the params.cfg file and run this script.

Requirements:
    - Python 3.7+
    - iqtree3 (IQ-TREE with alisim support)
    - AMAS.py (for alignment splitting and manipulation)
    - No Python packages beyond the standard library
"""

import argparse
import configparser
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path


# =============================================================================
# ANSI Color codes for terminal output (disabled if not a TTY)
# =============================================================================

class Colors:
    """Terminal colors for user-friendly output."""

    def __init__(self):
        use_color = hasattr(sys.stdout, "isatty") and sys.stdout.isatty()
        if use_color:
            self.RED = "\033[91m"
            self.GREEN = "\033[92m"
            self.YELLOW = "\033[93m"
            self.BLUE = "\033[94m"
            self.BOLD = "\033[1m"
            self.RESET = "\033[0m"
        else:
            self.RED = ""
            self.GREEN = ""
            self.YELLOW = ""
            self.BLUE = ""
            self.BOLD = ""
            self.RESET = ""


COLORS = Colors()


# =============================================================================
# Utility functions
# =============================================================================


def print_error(msg):
    """Print an error message in red and exit."""
    print(f"\n{COLORS.RED}{COLORS.BOLD}ERROR:{COLORS.RESET} {COLORS.RED}{msg}{COLORS.RESET}", file=sys.stderr)


def print_warning(msg):
    """Print a warning message in yellow."""
    print(f"{COLORS.YELLOW}WARNING:{COLORS.RESET} {msg}", file=sys.stderr)


def print_step(step_num, total, msg):
    """Print a step progress message."""
    print(f"\n{COLORS.BLUE}{COLORS.BOLD}[Step {step_num}/{total}]{COLORS.RESET} {msg}")


def print_success(msg):
    """Print a success message in green."""
    print(f"{COLORS.GREEN}✓{COLORS.RESET} {msg}")


def print_info(msg):
    """Print an informational message."""
    print(f"  {msg}")


def run_command(cmd, dry_run=False, verbose=False, description="", cwd=None):
    """Run a shell command with error handling.

    Parameters
    ----------
    cmd : list
        Command and arguments as a list.
    dry_run : bool
        If True, print command but don't execute.
    verbose : bool
        If True, print command before executing.
    description : str
        Human-readable description of what the command does.
    cwd : str or Path or None
        Working directory for the command (default: current directory).

    Returns
    -------
    subprocess.CompletedProcess or None
        The completed process result, or None in dry-run mode.

    Raises
    ------
    SystemExit
        If the command fails.
    """
    cmd_str = " ".join(str(c) for c in cmd)

    if description:
        print_info(f"{description}")

    if dry_run:
        cwd_msg = f" (in {cwd})" if cwd else ""
        print(f"  {COLORS.YELLOW}[DRY RUN]{COLORS.RESET} {cmd_str}{cwd_msg}")
        return None

    if verbose:
        cwd_msg = f" (in {cwd})" if cwd else ""
        print(f"  $ {cmd_str}{cwd_msg}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=86400,  # 24 hour timeout
            cwd=cwd,
        )
    except FileNotFoundError:
        print_error(
            f"Command not found: '{cmd[0]}'\n"
            f"Please ensure '{cmd[0]}' is installed and in your PATH.\n"
            f"Full command: {cmd_str}"
        )
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print_error(f"Command timed out after 24 hours:\n  {cmd_str}")
        sys.exit(1)
    except OSError as e:
        print_error(f"Failed to run command: {e}\n  Command: {cmd_str}")
        sys.exit(1)

    if result.returncode != 0:
        print_error(
            f"Command failed with exit code {result.returncode}:\n"
            f"  Command: {cmd_str}\n"
            f"  Stderr: {result.stderr.strip()[:500] if result.stderr else '(none)'}\n"
            f"  Stdout: {result.stdout.strip()[:500] if result.stdout else '(none)'}"
        )
        sys.exit(1)

    return result


def check_tool_available(tool_name, tool_path=None):
    """Check if a command-line tool is available.

    Parameters
    ----------
    tool_name : str
        Human-readable name of the tool.
    tool_path : str or None
        Path or command name to check.

    Returns
    -------
    str
        The resolved path/command for the tool.
    """
    if tool_path is None:
        tool_path = tool_name

    resolved = shutil.which(tool_path)
    if resolved is None:
        print_error(
            f"Required tool '{tool_name}' not found.\n"
            f"  Looked for: '{tool_path}'\n"
            f"  Please install it or provide the full path in your params file.\n"
            f"\n"
            f"  For iqtree3: https://github.com/iqtree/iqtree3\n"
            f"  For AMAS.py: pip install amas"
        )
        sys.exit(1)

    return resolved


def read_alignment_dimensions(filepath):
    """Read the number of taxa and sites from an alignment file.

    Supports PHYLIP and FASTA formats.

    Parameters
    ----------
    filepath : str or Path
        Path to the alignment file.

    Returns
    -------
    tuple of (int, int)
        Number of taxa and number of sites.
    """
    filepath = Path(filepath)
    ext = filepath.suffix.lower()

    if ext in (".phy", ".phylip"):
        with open(filepath, "r") as f:
            header = f.readline().strip().split()
            if len(header) < 2:
                print_error(
                    f"Invalid PHYLIP file: '{filepath}'\n"
                    f"  First line should contain 'ntaxa nsites' but got: '{' '.join(header)}'"
                )
                sys.exit(1)
            try:
                ntaxa = int(header[0])
                nsites = int(header[1])
            except ValueError:
                print_error(
                    f"Invalid PHYLIP header in '{filepath}'\n"
                    f"  Expected two integers (ntaxa nsites), got: '{' '.join(header)}'"
                )
                sys.exit(1)
        return ntaxa, nsites
    elif ext in (".fasta", ".fa", ".fas", ".fna", ".faa"):
        # Parse FASTA to count taxa and determine alignment length
        ntaxa = 0
        first_seq_len = None
        current_seq_len = 0
        in_seq = False

        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if in_seq and first_seq_len is None:
                        first_seq_len = current_seq_len
                    ntaxa += 1
                    current_seq_len = 0
                    in_seq = True
                else:
                    current_seq_len += len(line)

        if first_seq_len is None:
            first_seq_len = current_seq_len

        return ntaxa, first_seq_len
    else:
        print_error(
            f"Unrecognized alignment file format: '{filepath}'\n"
            f"  Supported extensions: .phy, .phylip, .fasta, .fa, .fas"
        )
        sys.exit(1)


def validate_tree_file(filepath):
    """Validate that a tree file exists and contains a Newick tree.

    Parameters
    ----------
    filepath : str or Path
        Path to the tree file.
    """
    filepath = Path(filepath)
    if not filepath.is_file():
        print_error(f"Tree file not found: '{filepath}'")
        sys.exit(1)

    with open(filepath, "r") as f:
        content = f.read().strip()

    if not content:
        print_error(f"Tree file is empty: '{filepath}'")
        sys.exit(1)

    # Basic Newick validation: should contain parentheses and end with semicolon
    if "(" not in content or ")" not in content:
        print_error(
            f"Tree file does not appear to contain a valid Newick tree: '{filepath}'\n"
            f"  A Newick tree should contain parentheses, e.g., ((A,B),(C,D));"
        )
        sys.exit(1)

    if not content.endswith(";"):
        print_warning(
            f"Tree file '{filepath}' does not end with a semicolon. "
            f"This may cause issues with iqtree3."
        )


def divide_into_parts(total, num_parts):
    """Divide a number into approximately equal parts.

    Parameters
    ----------
    total : int
        Total number to divide.
    num_parts : int
        Number of parts to divide into.

    Returns
    -------
    list of tuple
        List of (start, end) 1-indexed ranges.
    """
    if total < num_parts:
        print_error(
            f"Cannot divide {total} sites into {num_parts} partitions.\n"
            f"  The alignment must have at least {num_parts} sites."
        )
        sys.exit(1)

    base = total // num_parts
    remainder = total % num_parts

    ranges = []
    start = 1
    for i in range(num_parts):
        size = base + (1 if i < remainder else 0)
        end = start + size - 1
        ranges.append((start, end))
        start = end + 1

    return ranges


# =============================================================================
# Configuration parsing
# =============================================================================


def parse_config(config_path):
    """Parse the simulation parameters configuration file.

    Parameters
    ----------
    config_path : str or Path
        Path to the configuration file.

    Returns
    -------
    dict
        Dictionary containing all simulation parameters.
    """
    config_path = Path(config_path)
    if not config_path.is_file():
        print_error(
            f"Configuration file not found: '{config_path}'\n"
            f"  Please provide a valid params.cfg file.\n"
            f"  See examples/params.cfg for a template."
        )
        sys.exit(1)

    config = configparser.ConfigParser(
        interpolation=None,  # Disable interpolation to allow % in model strings
        comment_prefixes=("#",),
        inline_comment_prefixes=("#",),
    )

    try:
        config.read(config_path)
    except configparser.Error as e:
        print_error(
            f"Failed to parse configuration file '{config_path}':\n"
            f"  {e}\n"
            f"  Please check the file format."
        )
        sys.exit(1)

    params = {}

    # --- [general] section ---
    if "general" not in config:
        print_error(
            f"Configuration file is missing the [general] section.\n"
            f"  Please check your params file format.\n"
            f"  See examples/params.cfg for a template."
        )
        sys.exit(1)

    general = config["general"]

    # Required fields
    required_general = ["output_prefix", "ratio"]
    for field in required_general:
        if field not in general or not general[field].strip():
            print_error(
                f"Missing required field '{field}' in [general] section.\n"
                f"  Please set this value in your params file."
            )
            sys.exit(1)

    # alignment is optional if alignment_length is provided
    params["alignment"] = general.get("alignment", "").strip()
    params["alignment_length"] = general.get("alignment_length", "").strip()
    params["output_prefix"] = general["output_prefix"].strip()
    params["datatype"] = general.get("datatype", "aa").strip().lower()
    params["threads"] = general.get("threads", "1").strip()
    params["iqtree"] = general.get("iqtree", "iqtree3").strip()
    params["amas"] = general.get("amas", "AMAS.py").strip()
    params["output_dir"] = general.get("output_dir", "sim_output").strip()
    params["seed_base"] = general.get("seed_base", "").strip()
    params["introduce_gaps"] = general.get("introduce_gaps", "yes").strip().lower() in ("yes", "true", "1")
    params["gap_method"] = general.get("gap_method", "direct").strip().lower()
    params["concatenate"] = general.get("concatenate", "yes").strip().lower() in ("yes", "true", "1")

    # Indel model parameters (optional)
    params["indel"] = general.get("indel", "").strip()
    params["indel_size"] = general.get("indel_size", "").strip()

    # Parse num_partitions
    try:
        params["num_partitions"] = int(general.get("num_partitions", "5").strip())
    except ValueError:
        print_error(
            f"Invalid value for 'num_partitions' in [general] section.\n"
            f"  Expected a positive integer, got: '{general.get('num_partitions', '')}'"
        )
        sys.exit(1)

    if params["num_partitions"] < 1:
        print_error("num_partitions must be at least 1.")
        sys.exit(1)

    # Parse ratio
    ratio_str = general["ratio"].strip()
    try:
        parts = ratio_str.replace("-", ":").split(":")
        if len(parts) != 2:
            raise ValueError("Expected format like 70:30")
        ratio_a = int(parts[0])
        ratio_b = int(parts[1])
        if ratio_a <= 0 or ratio_b <= 0:
            raise ValueError("Both ratio values must be positive")
        params["ratio"] = (ratio_a, ratio_b)
    except ValueError as e:
        print_error(
            f"Invalid ratio format: '{ratio_str}'\n"
            f"  Expected format: 'majority:minority' (e.g., 70:30 or 60:40)\n"
            f"  {e}"
        )
        sys.exit(1)

    # Validate datatype
    if params["datatype"] not in ("aa", "dna"):
        print_error(
            f"Invalid datatype: '{params['datatype']}'\n"
            f"  Must be 'aa' (amino acid) or 'dna'"
        )
        sys.exit(1)

    # Validate threads
    try:
        threads_int = int(params["threads"])
        if threads_int < 1:
            raise ValueError()
    except ValueError:
        print_error(
            f"Invalid value for 'threads': '{params['threads']}'\n"
            f"  Must be a positive integer."
        )
        sys.exit(1)

    # Validate alignment vs alignment_length
    if not params["alignment"] and not params["alignment_length"]:
        print_error(
            "You must provide either 'alignment' (path to an existing alignment)\n"
            "  or 'alignment_length' (number of sites to simulate) in [general] section."
        )
        sys.exit(1)

    if params["alignment_length"]:
        try:
            params["alignment_length"] = int(params["alignment_length"])
            if params["alignment_length"] < 1:
                raise ValueError()
        except ValueError:
            print_error(
                f"Invalid value for 'alignment_length': '{params['alignment_length']}'\n"
                f"  Must be a positive integer."
            )
            sys.exit(1)
    else:
        params["alignment_length"] = None

    # Validate indel parameters if provided
    if params["indel"]:
        indel_parts = params["indel"].split(",")
        if len(indel_parts) != 2:
            print_error(
                f"Invalid indel format: '{params['indel']}'\n"
                f"  Expected format: insertion_rate,deletion_rate\n"
                f"  Example: 0.03,0.10"
            )
            sys.exit(1)
        try:
            ins_rate = float(indel_parts[0])
            del_rate = float(indel_parts[1])
            if ins_rate < 0 or del_rate < 0:
                raise ValueError()
        except ValueError:
            print_error(
                f"Invalid indel rates: '{params['indel']}'\n"
                f"  Both rates must be non-negative numbers.\n"
                f"  Example: 0.03,0.10"
            )
            sys.exit(1)

    # --- [tree1] and [tree2] sections ---
    for tree_key in ("tree1", "tree2"):
        if tree_key not in config:
            print_error(
                f"Configuration file is missing the [{tree_key}] section.\n"
                f"  You must specify two tree files (conflicting topologies)."
            )
            sys.exit(1)

        tree_section = config[tree_key]
        if "file" not in tree_section or not tree_section["file"].strip():
            print_error(f"Missing 'file' in [{tree_key}] section.")
            sys.exit(1)
        if "label" not in tree_section or not tree_section["label"].strip():
            print_error(f"Missing 'label' in [{tree_key}] section.")
            sys.exit(1)

        params[tree_key] = {
            "file": tree_section["file"].strip(),
            "label": tree_section["label"].strip(),
        }

    # Validate that labels don't contain problematic characters
    for tree_key in ("tree1", "tree2"):
        label = params[tree_key]["label"]
        if any(c in label for c in " /\\:*?\"<>|"):
            print_error(
                f"Tree label '{label}' contains invalid characters.\n"
                f"  Labels are used in filenames and should not contain:\n"
                f"  spaces, /, \\, :, *, ?, \", <, >, |"
            )
            sys.exit(1)

    # --- [models] section ---
    if "models" not in config:
        print_error(
            f"Configuration file is missing the [models] section.\n"
            f"  You must specify a model for each partition."
        )
        sys.exit(1)

    models = {}
    for key, value in config["models"].items():
        try:
            part_num = int(key)
        except ValueError:
            print_warning(f"Ignoring non-numeric key '{key}' in [models] section.")
            continue

        parts = [p.strip() for p in value.split(",")]
        if len(parts) != 3:
            print_error(
                f"Invalid model specification for partition {part_num}: '{value}'\n"
                f"  Expected format: model, gamma_shape, seed\n"
                f"  Example: WAG+C10, 1.0, 1"
            )
            sys.exit(1)

        model_name = parts[0]
        try:
            gamma_shape = float(parts[1])
        except ValueError:
            print_error(
                f"Invalid gamma shape for partition {part_num}: '{parts[1]}'\n"
                f"  Must be a positive number (e.g., 0.5, 1.0, 3.0)"
            )
            sys.exit(1)

        if gamma_shape <= 0:
            print_error(
                f"Gamma shape for partition {part_num} must be positive, got: {gamma_shape}"
            )
            sys.exit(1)

        try:
            seed = int(parts[2])
        except ValueError:
            print_error(
                f"Invalid seed for partition {part_num}: '{parts[2]}'\n"
                f"  Must be an integer."
            )
            sys.exit(1)

        models[part_num] = {
            "model": model_name,
            "gamma": gamma_shape,
            "seed": seed,
        }

    # Validate we have models for all partitions
    for i in range(1, params["num_partitions"] + 1):
        if i not in models:
            print_error(
                f"Missing model specification for partition {i} in [models] section.\n"
                f"  You specified num_partitions = {params['num_partitions']}, "
                f"so you need entries for partitions 1 through {params['num_partitions']}."
            )
            sys.exit(1)

    params["models"] = models

    # --- [slurm] section (optional) ---
    params["slurm"] = {"use_slurm": False}
    if "slurm" in config:
        slurm_section = config["slurm"]
        use_slurm = slurm_section.get("use_slurm", "no").strip().lower() in ("yes", "true", "1")
        params["slurm"] = {
            "use_slurm": use_slurm,
            "partition": slurm_section.get("partition", "").strip(),
            "time": slurm_section.get("time", "24:00:00").strip(),
            "memory": slurm_section.get("memory", "32GB").strip(),
            "nodes": slurm_section.get("nodes", "1").strip(),
            "email": slurm_section.get("email", "").strip(),
            "mail_type": slurm_section.get("mail_type", "ALL").strip(),
        }

    return params


# =============================================================================
# Core pipeline steps
# =============================================================================


def write_split_file(filepath, partitions):
    """Write an AMAS-compatible partition file.

    Parameters
    ----------
    filepath : str or Path
        Output file path.
    partitions : list of tuple
        List of (name, start, end) tuples.
    """
    with open(filepath, "w") as f:
        for name, start, end in partitions:
            f.write(f"{name} = {start}-{end}\n")


def generate_slurm_script(script_path, commands, params, tree_label):
    """Generate a SLURM batch script for alisim commands.

    Parameters
    ----------
    script_path : str or Path
        Output path for the SLURM script.
    commands : list of list
        List of commands (each command is a list of strings).
    params : dict
        Simulation parameters.
    tree_label : str
        Label for the tree topology.
    """
    slurm = params["slurm"]
    job_name = f"alisim_{tree_label}_{params['ratio'][0]}-{params['ratio'][1]}"

    with open(script_path, "w") as f:
        f.write("#!/usr/bin/env bash\n\n")
        f.write(f"#SBATCH --time={slurm['time']}\n")
        f.write(f"#SBATCH --nodes={slurm['nodes']}\n")
        f.write("#SBATCH --ntasks-per-core=2\n")
        f.write("#SBATCH --hint=multithread\n")
        if slurm["partition"]:
            f.write(f"#SBATCH --partition={slurm['partition']}\n")
        f.write(f"#SBATCH --ntasks={params['threads']}\n")
        f.write(f"#SBATCH --mem={slurm['memory']}\n")
        f.write(f"#SBATCH --job-name={job_name}\n")
        if slurm["email"]:
            f.write(f"#SBATCH --mail-type={slurm['mail_type']}\n")
            f.write(f"#SBATCH --mail-user={slurm['email']}\n")
        f.write("\n# cd to dir from which I submit job\n")
        f.write("cd $SLURM_SUBMIT_DIR\n\n")
        f.write("# alisim commands\n\n")

        for cmd in commands:
            f.write(" ".join(str(c) for c in cmd) + "\n\n")

    os.chmod(script_path, 0o755)


def build_alisim_command(params, partition_num, source_alignment, tree_file, tree_label,
                         partition_length=None):
    """Build an alisim command for a single partition.

    Parameters
    ----------
    params : dict
        Simulation parameters.
    partition_num : int
        Partition number (1-indexed).
    source_alignment : str or None
        Path to the source alignment for this partition. None if using --length mode.
    tree_file : str
        Path to the tree file.
    tree_label : str
        Label for the tree topology.
    partition_length : int or None
        Number of sites for this partition (used when source_alignment is None).

    Returns
    -------
    tuple of (list, str)
        The command as a list and the output alignment name.
    """
    model_info = params["models"][partition_num]
    ratio_str = f"{params['ratio'][0]}-{params['ratio'][1]}"

    # Build model string with gamma
    model_with_gamma = f"{model_info['model']}+G{{{model_info['gamma']}}}"

    # Build a safe model label for filenames (replace problematic chars)
    # Format gamma: remove trailing zeros after decimal, then remove the dot
    # 1.0 -> "1", 0.5 -> "05", 0.2 -> "02", 3.0 -> "3", 0.8 -> "08"
    gamma_val = model_info["gamma"]
    if gamma_val == int(gamma_val):
        gamma_str = str(int(gamma_val))
    else:
        # Format as minimal decimal, then remove the dot
        gamma_str = f"{gamma_val:g}".replace(".", "")
    model_label = model_info["model"].replace(".", "") + "+G" + gamma_str

    # Output name
    output_name = f"{params['output_prefix']}_{tree_label}_{ratio_str}_{model_label}"
    prefix = output_name

    cmd = [
        params["iqtree"],
        "--alisim", output_name,
    ]

    if source_alignment is not None:
        # Use source alignment for site frequencies and rates
        cmd.extend(["-s", str(source_alignment)])
        cmd.extend(["-te", str(tree_file)])
        cmd.extend(["-m", model_with_gamma])
        cmd.extend(["-pre", prefix])
        cmd.extend(["--site-freq", "SAMPLING"])
        cmd.extend(["--site-rate", "SAMPLING"])
    else:
        # No source alignment: use --length and tree
        cmd.extend(["-t", str(tree_file)])
        cmd.extend(["-m", model_with_gamma])
        cmd.extend(["-pre", prefix])
        if partition_length is not None:
            cmd.extend(["--length", str(partition_length)])

    # Add indel model if specified
    if params.get("indel"):
        cmd.extend(["--indel", params["indel"]])
    if params.get("indel_size"):
        cmd.extend(["--indel-size", params["indel_size"]])

    cmd.extend(["-nt", str(params["threads"])])
    cmd.extend(["--seed", str(model_info["seed"])])

    if source_alignment is not None:
        cmd.append("-blfix")

    return cmd, output_name


def introduce_gaps_into_alignment(reference_path, simulated_path, output_path, method="direct"):
    """Introduce gaps from a reference alignment into a simulated alignment.

    This is a self-contained implementation (not calling the external script)
    to avoid dependency on script location.

    Parameters
    ----------
    reference_path : str or Path
        Path to the reference alignment with gap pattern.
    simulated_path : str or Path
        Path to the simulated alignment (no gaps).
    output_path : str or Path
        Path for the output alignment with gaps.
    method : str
        'direct' or 'proportional' column mapping.

    Returns
    -------
    bool
        True if successful.
    """
    # Parse alignments
    ref_seqs, ref_taxa = _parse_alignment_file(reference_path)
    sim_seqs, sim_taxa = _parse_alignment_file(simulated_path)

    if not ref_seqs or not sim_seqs:
        print_warning(f"Could not introduce gaps: empty alignment(s)")
        return False

    # Find common taxa
    common = set(ref_taxa) & set(sim_taxa)
    if not common:
        print_warning(
            f"No common taxa between reference and simulated alignment. "
            f"Skipping gap introduction."
        )
        return False

    ref_length = len(ref_seqs[ref_taxa[0]])
    sim_length = len(sim_seqs[sim_taxa[0]])

    gap_chars = set("-?Xx")

    # Extract gap mask from reference
    gap_mask = []
    for col in range(ref_length):
        missing_taxa = set()
        for taxon in ref_taxa:
            if col < len(ref_seqs[taxon]) and ref_seqs[taxon][col] in gap_chars:
                missing_taxa.add(taxon)
        gap_mask.append(missing_taxa)

    # Apply gap mask to simulated
    modified = {}
    for taxon in sim_taxa:
        seq_list = list(sim_seqs[taxon])
        if method == "proportional":
            for col in range(sim_length):
                ref_col = int(col * ref_length / sim_length) if sim_length > 0 else 0
                ref_col = min(ref_col, ref_length - 1)
                if taxon in gap_mask[ref_col]:
                    seq_list[col] = "-"
        else:
            for col in range(min(sim_length, ref_length)):
                if taxon in gap_mask[col]:
                    seq_list[col] = "-"
        modified[taxon] = "".join(seq_list)

    # Write output in the same format as the simulated input
    _write_alignment_file(modified, sim_taxa, output_path)
    return True


def _parse_alignment_file(filepath):
    """Parse an alignment file (PHYLIP or FASTA).

    Parameters
    ----------
    filepath : str or Path
        Path to the alignment file.

    Returns
    -------
    tuple of (dict, list)
        Sequences dict and taxa order list.
    """
    filepath = Path(filepath)
    ext = filepath.suffix.lower()

    if ext in (".phy", ".phylip"):
        return _parse_phylip(filepath)
    else:
        return _parse_fasta(filepath)


def _parse_phylip(filepath):
    """Parse a PHYLIP alignment file."""
    sequences = {}
    taxa_order = []

    with open(filepath, "r") as f:
        lines = [line.rstrip() for line in f if line.strip()]

    if not lines:
        return {}, []

    header = lines[0].split()
    ntaxa = int(header[0])
    nchar = int(header[1])

    idx = 1
    for i in range(ntaxa):
        if idx >= len(lines):
            break
        parts = lines[idx].split(None, 1)
        taxon = parts[0]
        seq = parts[1].replace(" ", "") if len(parts) > 1 else ""
        taxa_order.append(taxon)
        sequences[taxon] = seq
        idx += 1

    # Handle interleaved format
    first_seq_len = len(sequences[taxa_order[0]]) if taxa_order else 0
    if first_seq_len < nchar and idx < len(lines):
        while idx < len(lines):
            for taxon in taxa_order:
                if idx < len(lines):
                    sequences[taxon] += lines[idx].replace(" ", "")
                    idx += 1

    for taxon in taxa_order:
        sequences[taxon] = sequences[taxon][:nchar]

    return sequences, taxa_order


def _parse_fasta(filepath):
    """Parse a FASTA alignment file."""
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
                taxa_order.append(current_taxon)
                current_seq = []
            else:
                current_seq.append(line)
        if current_taxon is not None:
            sequences[current_taxon] = "".join(current_seq)

    return sequences, taxa_order


def _write_alignment_file(sequences, taxa_order, filepath):
    """Write an alignment file (format detected from extension).

    Parameters
    ----------
    sequences : dict
        Sequences dict.
    taxa_order : list
        Ordered taxa list.
    filepath : str or Path
        Output file path.
    """
    filepath = Path(filepath)
    ext = filepath.suffix.lower()

    if ext in (".phy", ".phylip"):
        ntaxa = len(taxa_order)
        nchar = len(sequences[taxa_order[0]]) if taxa_order else 0
        with open(filepath, "w") as f:
            f.write(f"{ntaxa} {nchar}\n")
            for taxon in taxa_order:
                f.write(f"{taxon} {sequences[taxon]}\n")
    else:
        with open(filepath, "w") as f:
            for taxon in taxa_order:
                f.write(f">{taxon}\n")
                seq = sequences[taxon]
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i + 80] + "\n")


def concatenate_phylip_files(input_files, output_file):
    """Concatenate multiple PHYLIP alignment files.

    All files must have the same taxa (in the same order).
    This is a simple concatenation without requiring FASconCAT.

    Parameters
    ----------
    input_files : list of str/Path
        Paths to the PHYLIP files to concatenate.
    output_file : str or Path
        Path for the concatenated output.

    Returns
    -------
    bool
        True if successful.
    """
    if not input_files:
        print_warning("No files to concatenate.")
        return False

    all_sequences = {}
    taxa_order = None
    total_length = 0

    for fpath in input_files:
        seqs, taxa = _parse_alignment_file(fpath)
        if not seqs:
            print_warning(f"Skipping empty file during concatenation: {fpath}")
            continue

        if taxa_order is None:
            taxa_order = taxa
            all_sequences = {t: "" for t in taxa}
        else:
            # Verify taxa match
            if set(taxa) != set(taxa_order):
                print_warning(
                    f"Taxa mismatch in file {fpath}. "
                    f"Expected {len(taxa_order)} taxa, found {len(taxa)}. "
                    f"Skipping this file."
                )
                continue

        seq_len = len(seqs[taxa[0]])
        total_length += seq_len
        for taxon in taxa_order:
            if taxon in seqs:
                all_sequences[taxon] += seqs[taxon]
            else:
                # Fill with gaps if taxon is missing
                all_sequences[taxon] += "-" * seq_len

    if taxa_order is None:
        print_warning("No valid files found for concatenation.")
        return False

    _write_alignment_file(all_sequences, taxa_order, output_file)
    return True


# =============================================================================
# Main pipeline
# =============================================================================


def run_pipeline(params, dry_run=False, verbose=False):
    """Execute the full simulation pipeline.

    Parameters
    ----------
    params : dict
        Parsed simulation parameters.
    dry_run : bool
        If True, print commands without executing.
    verbose : bool
        If True, print additional detail.
    """
    has_alignment = bool(params["alignment"])

    total_steps = 7
    if not has_alignment:
        # Skip AMAS split steps (3 and 4) when no alignment
        total_steps -= 2
    if not params["introduce_gaps"] or not has_alignment:
        # Can't introduce empirical gaps without a source alignment
        if has_alignment and not params["introduce_gaps"]:
            total_steps -= 1
        elif not has_alignment:
            total_steps -= 1
    if not params["concatenate"]:
        total_steps -= 1

    ratio_a, ratio_b = params["ratio"]
    ratio_str = f"{ratio_a}-{ratio_b}"
    ratio_total = ratio_a + ratio_b

    # =========================================================================
    # Step 1: Validate inputs
    # =========================================================================
    print_step(1, total_steps, "Validating inputs")

    alignment_path = None
    nsites = None

    if has_alignment:
        # Check alignment file
        alignment_path = Path(params["alignment"])
        if not alignment_path.is_file():
            print_error(
                f"Alignment file not found: '{alignment_path}'\n"
                f"  Please check the path in your params file."
            )
            sys.exit(1)

        ntaxa, nsites = read_alignment_dimensions(alignment_path)
        print_success(f"Alignment: {ntaxa} taxa, {nsites} sites ({alignment_path})")
    else:
        nsites = params["alignment_length"]
        print_success(f"No source alignment; using alignment_length = {nsites} sites")
        if params["indel"]:
            print_success(f"Indel model: --indel {params['indel']}")
        if params["indel_size"]:
            print_success(f"Indel size distribution: --indel-size {params['indel_size']}")

    # Check tree files
    for tree_key in ("tree1", "tree2"):
        tree_path = Path(params[tree_key]["file"])
        validate_tree_file(tree_path)
        print_success(f"{tree_key}: {params[tree_key]['label']} ({tree_path})")

    # Check tools (skip if dry run and tools might not be on this machine)
    if not dry_run:
        check_tool_available("iqtree3", params["iqtree"])
        if has_alignment:
            check_tool_available("AMAS.py", params["amas"])
            print_success("Required tools found (iqtree3, AMAS.py)")
        else:
            print_success("Required tools found (iqtree3)")
    else:
        print_info("Skipping tool availability check in dry-run mode")

    # Validate ratio makes sense with alignment length
    sites_tree1 = int(math.floor(nsites * ratio_a / ratio_total))
    sites_tree2 = nsites - sites_tree1

    if sites_tree1 < params["num_partitions"] or sites_tree2 < params["num_partitions"]:
        print_error(
            f"Alignment too short ({nsites} sites) for the specified ratio ({ratio_str}) "
            f"and num_partitions ({params['num_partitions']}).\n"
            f"  tree1 would get {sites_tree1} sites, tree2 would get {sites_tree2} sites.\n"
            f"  Each must have at least {params['num_partitions']} sites."
        )
        sys.exit(1)

    print_success(
        f"Ratio {ratio_str}: {sites_tree1} sites for {params['tree1']['label']}, "
        f"{sites_tree2} sites for {params['tree2']['label']}"
    )

    # =========================================================================
    # Step 2: Create output directory structure
    # =========================================================================
    current_step = 2
    print_step(current_step, total_steps, "Creating output directory structure")

    output_dir = Path(params["output_dir"])
    tree1_dir = output_dir / params["tree1"]["label"]
    tree2_dir = output_dir / params["tree2"]["label"]
    combined_dir = output_dir / "combined"

    for d in [output_dir, tree1_dir, tree2_dir, combined_dir]:
        if not dry_run:
            d.mkdir(parents=True, exist_ok=True)
        print_info(f"Directory: {d}")

    print_success("Output directories ready")

    # Use absolute paths to avoid confusion with working directories
    if alignment_path is not None:
        alignment_path = alignment_path.resolve()
    output_dir = output_dir.resolve()
    tree1_dir = tree1_dir.resolve()
    tree2_dir = tree2_dir.resolve()
    combined_dir = combined_dir.resolve()

    # =========================================================================
    # Steps 3-4: Split alignment (only when source alignment is provided)
    # =========================================================================
    if has_alignment:
        current_step += 1
        print_step(current_step, total_steps, f"Splitting alignment by ratio ({ratio_str})")

        # Create split file for AMAS
        split_file = output_dir / "ratio_split.txt"
        partitions_split = [
            (params["tree1"]["label"], 1, sites_tree1),
            (params["tree2"]["label"], sites_tree1 + 1, nsites),
        ]

        if not dry_run:
            write_split_file(split_file, partitions_split)
        print_info(f"Split file: {split_file}")

        # Run AMAS split (run from output_dir so output lands there)
        amas_cmd = [
            params["amas"],
            "split",
            "-l", str(split_file),
            "-u", "phylip",
            "-i", str(alignment_path),
            "-f", "phylip",
            "-d", params["datatype"],
        ]
        run_command(amas_cmd, dry_run=dry_run, verbose=verbose,
                    description="Running AMAS split...", cwd=str(output_dir))

        # AMAS outputs files with suffixes based on partition names
        alignment_stem = alignment_path.stem
        tree1_alignment_name = f"{alignment_stem}_{params['tree1']['label']}-out.phy"
        tree2_alignment_name = f"{alignment_stem}_{params['tree2']['label']}-out.phy"
        tree1_alignment_in_outdir = output_dir / tree1_alignment_name
        tree2_alignment_in_outdir = output_dir / tree2_alignment_name

        # Move the split files to their directories
        if not dry_run:
            if tree1_alignment_in_outdir.exists():
                shutil.move(str(tree1_alignment_in_outdir), str(tree1_dir / tree1_alignment_name))
            else:
                print_error(
                    f"Expected AMAS output file not found: '{tree1_alignment_in_outdir}'\n"
                    f"  AMAS may have used a different naming convention.\n"
                    f"  Please check the output directory for split files."
                )
                sys.exit(1)

            if tree2_alignment_in_outdir.exists():
                shutil.move(str(tree2_alignment_in_outdir), str(tree2_dir / tree2_alignment_name))
            else:
                print_error(f"Expected AMAS output file not found: '{tree2_alignment_in_outdir}'")
                sys.exit(1)

        tree1_alignment_path = tree1_dir / tree1_alignment_name
        tree2_alignment_path = tree2_dir / tree2_alignment_name

        print_success(f"Split complete: {tree1_alignment_name} ({sites_tree1} sites), {tree2_alignment_name} ({sites_tree2} sites)")

        # Split each portion into sub-partitions
        current_step += 1
        print_step(current_step, total_steps, f"Splitting each portion into {params['num_partitions']} sub-partitions")

        # Split tree1 portion
        tree1_ranges = divide_into_parts(sites_tree1, params["num_partitions"])
        tree1_split_file = tree1_dir / "sub_split.txt"
        tree1_partitions = [
            (f"gene{i + 1}", start, end) for i, (start, end) in enumerate(tree1_ranges)
        ]
        if not dry_run:
            write_split_file(tree1_split_file, tree1_partitions)

        amas_cmd_tree1 = [
            params["amas"],
            "split",
            "-l", str(tree1_split_file),
            "-u", "phylip",
            "-i", str(tree1_alignment_path),
            "-f", "phylip",
            "-d", params["datatype"],
        ]
        run_command(amas_cmd_tree1, dry_run=dry_run, verbose=verbose,
                    description=f"Splitting {params['tree1']['label']} portion into {params['num_partitions']} parts...",
                    cwd=str(tree1_dir))

        # Verify AMAS output files exist in tree1 directory
        if not dry_run:
            tree1_al_stem = Path(tree1_alignment_name).stem
            for i in range(1, params["num_partitions"] + 1):
                gene_file = tree1_dir / f"{tree1_al_stem}_gene{i}-out.phy"
                if not gene_file.exists():
                    print_error(
                        f"Expected sub-partition file not found: '{gene_file}'\n"
                        f"  AMAS split may have failed."
                    )
                    sys.exit(1)

        # Split tree2 portion
        tree2_ranges = divide_into_parts(sites_tree2, params["num_partitions"])
        tree2_split_file = tree2_dir / "sub_split.txt"
        tree2_partitions = [
            (f"gene{i + 1}", start, end) for i, (start, end) in enumerate(tree2_ranges)
        ]
        if not dry_run:
            write_split_file(tree2_split_file, tree2_partitions)

        amas_cmd_tree2 = [
            params["amas"],
            "split",
            "-l", str(tree2_split_file),
            "-u", "phylip",
            "-i", str(tree2_alignment_path),
            "-f", "phylip",
            "-d", params["datatype"],
        ]
        run_command(amas_cmd_tree2, dry_run=dry_run, verbose=verbose,
                    description=f"Splitting {params['tree2']['label']} portion into {params['num_partitions']} parts...",
                    cwd=str(tree2_dir))

        # Verify AMAS output files exist in tree2 directory
        if not dry_run:
            tree2_al_stem = Path(tree2_alignment_name).stem
            for i in range(1, params["num_partitions"] + 1):
                gene_file = tree2_dir / f"{tree2_al_stem}_gene{i}-out.phy"
                if not gene_file.exists():
                    print_error(
                        f"Expected sub-partition file not found: '{gene_file}'\n"
                        f"  AMAS split may have failed."
                    )
                    sys.exit(1)

        print_success(f"Sub-partitions created for both topologies")

    # =========================================================================
    # Run alisim on each sub-partition (or directly with --length)
    # =========================================================================
    current_step += 1
    print_step(current_step, total_steps, "Running alisim simulations")

    simulated_files = {"tree1": [], "tree2": []}

    if has_alignment:
        # With source alignment: use sub-partition files
        tree1_al_stem = Path(tree1_alignment_name).stem
        tree2_al_stem = Path(tree2_alignment_name).stem

        for tree_key, tree_dir_path, al_stem in [
            ("tree1", tree1_dir, tree1_al_stem),
            ("tree2", tree2_dir, tree2_al_stem),
        ]:
            tree_label = params[tree_key]["label"]
            tree_file = Path(params[tree_key]["file"]).resolve()
            alisim_commands = []

            for i in range(1, params["num_partitions"] + 1):
                source_alignment = tree_dir_path / f"{al_stem}_gene{i}-out.phy"
                cmd, output_name = build_alisim_command(
                    params, i, source_alignment, tree_file, tree_label
                )
                alisim_commands.append(cmd)
                simulated_files[tree_key].append(output_name + ".phy")

            if params["slurm"]["use_slurm"]:
                slurm_script = tree_dir_path / f"alisim_{tree_label}_{ratio_str}.sh"
                if not dry_run:
                    generate_slurm_script(slurm_script, alisim_commands, params, tree_label)
                print_info(f"SLURM script generated: {slurm_script}")
                print_info(f"  Submit with: sbatch {slurm_script}")
            else:
                for i, cmd in enumerate(alisim_commands, 1):
                    run_command(
                        cmd, dry_run=dry_run, verbose=verbose,
                        description=f"Simulating {tree_label} partition {i}/{params['num_partitions']}..."
                    )
    else:
        # Without source alignment: use --length for each partition
        tree1_part_ranges = divide_into_parts(sites_tree1, params["num_partitions"])
        tree2_part_ranges = divide_into_parts(sites_tree2, params["num_partitions"])

        for tree_key, tree_dir_path, part_ranges in [
            ("tree1", tree1_dir, tree1_part_ranges),
            ("tree2", tree2_dir, tree2_part_ranges),
        ]:
            tree_label = params[tree_key]["label"]
            tree_file = Path(params[tree_key]["file"]).resolve()
            alisim_commands = []

            for i in range(1, params["num_partitions"] + 1):
                part_start, part_end = part_ranges[i - 1]
                part_length = part_end - part_start + 1
                cmd, output_name = build_alisim_command(
                    params, i, None, tree_file, tree_label,
                    partition_length=part_length
                )
                alisim_commands.append(cmd)
                simulated_files[tree_key].append(output_name + ".phy")

            if params["slurm"]["use_slurm"]:
                slurm_script = tree_dir_path / f"alisim_{tree_label}_{ratio_str}.sh"
                if not dry_run:
                    generate_slurm_script(slurm_script, alisim_commands, params, tree_label)
                print_info(f"SLURM script generated: {slurm_script}")
                print_info(f"  Submit with: sbatch {slurm_script}")
            else:
                for i, cmd in enumerate(alisim_commands, 1):
                    run_command(
                        cmd, dry_run=dry_run, verbose=verbose,
                        description=f"Simulating {tree_label} partition {i}/{params['num_partitions']}..."
                    )

    if params["slurm"]["use_slurm"]:
        print_success("SLURM scripts generated. Submit them to run simulations.")
        print_warning(
            "When SLURM jobs are complete, re-run this script with --post-simulation\n"
            "  to introduce gaps and concatenate results."
        )
        return
    else:
        print_success("All alisim simulations complete")

    # =========================================================================
    # Introduce gaps (optional, only when source alignment is available)
    # =========================================================================
    if params["introduce_gaps"] and has_alignment:
        current_step += 1
        print_step(current_step, total_steps, "Introducing empirical gap patterns")

        for tree_key in ("tree1", "tree2"):
            tree_label = params[tree_key]["label"]
            tree_dir_path = tree1_dir if tree_key == "tree1" else tree2_dir

            for sim_file in simulated_files[tree_key]:
                sim_path = Path(sim_file)
                if not sim_path.exists():
                    # Try in tree directory
                    sim_path = tree_dir_path / sim_file
                if not sim_path.exists() and not dry_run:
                    print_warning(f"Simulated file not found: {sim_file}. Skipping gap introduction.")
                    continue

                output_gap_file = combined_dir / sim_path.name
                if not dry_run:
                    success = introduce_gaps_into_alignment(
                        str(alignment_path),
                        str(sim_path),
                        str(output_gap_file),
                        method=params["gap_method"],
                    )
                    if success:
                        print_info(f"Gaps introduced: {output_gap_file.name}")
                    else:
                        # Copy without gaps as fallback
                        shutil.copy2(str(sim_path), str(output_gap_file))
                        print_info(f"Copied without gaps: {output_gap_file.name}")
                else:
                    print_info(f"[DRY RUN] Would introduce gaps: {sim_file} -> {combined_dir / Path(sim_file).name}")

        print_success("Gap introduction complete")
    else:
        # Copy simulated files directly to combined directory
        for tree_key in ("tree1", "tree2"):
            tree_dir_path = tree1_dir if tree_key == "tree1" else tree2_dir
            for sim_file in simulated_files[tree_key]:
                sim_path = Path(sim_file)
                if not sim_path.exists():
                    sim_path = tree_dir_path / sim_file
                if sim_path.exists() and not dry_run:
                    shutil.copy2(str(sim_path), str(combined_dir / sim_path.name))

    # =========================================================================
    # Step 7: Concatenate final results (optional)
    # =========================================================================
    if params["concatenate"]:
        current_step += 1
        print_step(current_step, total_steps, "Concatenating simulated partitions")

        # Collect all files in combined directory
        combined_files = sorted(combined_dir.glob("*.phy")) if not dry_run else []

        if dry_run:
            print_info("[DRY RUN] Would concatenate all .phy files in combined/ directory")
        elif combined_files:
            concat_output = output_dir / f"{params['output_prefix']}_{ratio_str}_concatenated.phy"
            success = concatenate_phylip_files(combined_files, str(concat_output))
            if success:
                print_success(f"Concatenated alignment: {concat_output}")
            else:
                print_warning("Concatenation failed. Individual partition files are still available.")
        else:
            print_warning("No files found to concatenate in combined/ directory.")

    # =========================================================================
    # Summary
    # =========================================================================
    print(f"\n{'=' * 60}")
    print(f"{COLORS.GREEN}{COLORS.BOLD}SIMULATION COMPLETE{COLORS.RESET}")
    print(f"{'=' * 60}")
    print(f"  Output directory: {output_dir}")
    print(f"  Ratio: {ratio_str} ({params['tree1']['label']} : {params['tree2']['label']})")
    print(f"  Sites per topology: {sites_tree1} ({params['tree1']['label']}), {sites_tree2} ({params['tree2']['label']})")
    print(f"  Partitions per topology: {params['num_partitions']}")
    print(f"  Models used: {', '.join(params['models'][i]['model'] for i in range(1, params['num_partitions'] + 1))}")
    if params.get("indel"):
        print(f"  Indel model: {params['indel']}")
    if params.get("indel_size"):
        print(f"  Indel size distribution: {params['indel_size']}")
    if params["introduce_gaps"] and has_alignment:
        print(f"  Gap introduction: {params['gap_method']} mapping")
    if params["concatenate"]:
        concat_file = output_dir / f"{params['output_prefix']}_{ratio_str}_concatenated.phy"
        print(f"  Final concatenated file: {concat_file}")
    print(f"{'=' * 60}")


# =============================================================================
# Entry point
# =============================================================================


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Automated signal-conflict simulation pipeline using alisim (IQ-TREE).\n\n"
            "This script automates the creation of simulated phylogenetic data with\n"
            "conflicting phylogenetic signal, using a configuration file to specify\n"
            "all simulation parameters."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run the full simulation pipeline
  python simulate_conflict.py params.cfg

  # Preview all commands without running them (recommended first!)
  python simulate_conflict.py params.cfg --dry-run

  # Run with verbose output (shows all commands)
  python simulate_conflict.py params.cfg --verbose

  # Continue after SLURM jobs finish (gap introduction + concatenation)
  python simulate_conflict.py params.cfg --post-simulation

For help creating a params file, see: examples/params.cfg
""",
    )

    parser.add_argument(
        "config",
        help="Path to the simulation parameters file (e.g., params.cfg)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print all commands without executing them (use this to check your setup!)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed progress information",
    )
    parser.add_argument(
        "--post-simulation",
        action="store_true",
        help=(
            "Run only the post-simulation steps (gap introduction and concatenation). "
            "Use this after SLURM simulation jobs have completed."
        ),
    )

    args = parser.parse_args()

    # Print header
    print(f"\n{COLORS.BOLD}{'=' * 60}{COLORS.RESET}")
    print(f"{COLORS.BOLD}  Signal-Conflict Simulation Pipeline (alisim){COLORS.RESET}")
    print(f"{COLORS.BOLD}{'=' * 60}{COLORS.RESET}")

    if args.dry_run:
        print(f"\n{COLORS.YELLOW}{COLORS.BOLD}  *** DRY RUN MODE - No commands will be executed ***{COLORS.RESET}")

    # Parse configuration
    print(f"\n  Reading configuration: {args.config}")
    params = parse_config(args.config)
    print_success("Configuration parsed successfully")

    # Print summary of parameters
    print(f"\n  {'—' * 40}")
    if params["alignment"]:
        print(f"  Alignment:    {params['alignment']}")
    else:
        print(f"  Alignment:    (none - using alignment_length = {params['alignment_length']})")
    print(f"  Tree 1:       {params['tree1']['label']} ({params['tree1']['file']})")
    print(f"  Tree 2:       {params['tree2']['label']} ({params['tree2']['file']})")
    print(f"  Ratio:        {params['ratio'][0]}:{params['ratio'][1]}")
    print(f"  Partitions:   {params['num_partitions']}")
    if params["indel"]:
        print(f"  Indel:        {params['indel']}")
    if params["indel_size"]:
        print(f"  Indel size:   {params['indel_size']}")
    print(f"  Output:       {params['output_dir']}/")
    print(f"  {'—' * 40}")

    # Run the pipeline
    if args.post_simulation:
        print_warning("--post-simulation mode not yet fully implemented. Running full pipeline.")

    run_pipeline(params, dry_run=args.dry_run, verbose=args.verbose)


if __name__ == "__main__":
    main()
