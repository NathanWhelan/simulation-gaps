#!/usr/bin/env python3
import sys


def divide_into_fifths(total_number, pad_zeros=False):
    base = total_number // 5
    remainder = total_number % 5

    # Determine the character width of the largest raw number string
    max_digit_len = len(str(total_number))

    output_lines = []
    start = 1
    for i in range(1, 6):
        size = base + (1 if (i - 1) < remainder else 0)
        end = start + size - 1

        if pad_zeros:
            # Pad both numbers to the exact maximum digit width
            start_str = f"{start:0{max_digit_len}d}"
            end_str = f"{end:0{max_digit_len}d}"
        else:
            start_str = str(start)
            end_str = str(end)

        line = f"gene{i} = {start_str}-{end_str}"
        output_lines.append(line)
        start = end + 1

    # Fix: If padding is ON, all lines are structurally guaranteed to be identical length.
    # If padding is OFF, we spaces-pad the right side of shorter rows so the overall text file
    # maintains completely uniform row character counts.
    max_line_len = max(len(line) for line in output_lines)
    final_output = [line.ljust(max_line_len) for line in output_lines]

    # Write the output to a file in the working directory
    output_filename = "five-split.txt"
    try:
        with open(output_filename, "w") as f:
            f.write("\n".join(final_output) + "\n")
        print(f"Success: Output written to '{output_filename}'")
    except IOError as e:
        print(f"Error: Could not write to file '{output_filename}': {e}")
        sys.exit(1)


if __name__ == "__main__":
    # Clean argument parsing that ignores flag order
    args = sys.argv[1:]

    use_padding = "--padded" in args
    if use_padding:
        args.remove("--padded")

    # Check if we have our primary number argument left
    if not args:
        print("Error: Missing input number.")
        print("Usage: python script.py <number> [--padded]")
        sys.exit(1)

    # Validate that the remaining argument is a clean integer
    input_val = args[0]
    if not input_val.isdigit():
        print(
            f"Error: '{input_val}' is not a valid positive integer. Please enter a number greater than or equal to 5."
        )
        sys.exit(1)

    target_number = int(input_val)

    if target_number < 5:
        print("Error: The number must be 5 or greater to divide into fifths.")
        sys.exit(1)

    divide_into_fifths(target_number, pad_zeros=use_padding)
