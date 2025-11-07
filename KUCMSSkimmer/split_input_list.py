#!/usr/bin/env python3
import os
import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: python split_condor_lists.py <input_file> <item_to_add>")
        sys.exit(1)

    input_file = sys.argv[1]
    #new_item1 = bla

    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)

    with open(input_file, "r") as infile:
        for line in infile:
            line = line.strip()

            # Skip blank lines, comments, or lines starting with space/# 
            if not line or line.startswith("#"):
                continue

            parts = line.split()

            # Expect at least two entries per line
            if len(parts) < 9:
                print(f"Skipping malformed line: {line}")
                continue

            output_file = parts[1]
            if output_file.endswith(".txt"):
                output_file = output_file.replace(".txt", "_condor.txt")
            else:
                output_file = output_file + "_condor.txt"

            # Add the new item to the end of the list
            #parts.append(new_item)

            # Write modified line to new file
            with open(output_file, "w") as outfile:
                outfile.write(" ".join(parts) + "\n")

            print(f"[OK] Created {output_file}")

if __name__ == "__main__":
    main()

