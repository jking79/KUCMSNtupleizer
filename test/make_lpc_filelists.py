#!/usr/bin/env python3

import os
import sys

directory = "sms_fullpath_filelist"
prefix_to_remove = "root://cmseos.fnal.gov/"


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <match_string>")
        print(f"Example: {sys.argv[0]} filelist_T5Zg")
        sys.exit(1)

    match_string = sys.argv[1]

    if not os.path.isdir(directory):
        print(f"ERROR: directory does not exist: {directory}")
        sys.exit(1)

    # Find matching input files.
    input_files = sorted(
        filename
        for filename in os.listdir(directory)
        if filename.startswith(match_string)
        and os.path.isfile(os.path.join(directory, filename))
    )

    if not input_files:
        print(
            f'No files starting with "{match_string}" '
            f"found in {directory}/"
        )
        sys.exit(0)

    print(f'Found {len(input_files)} file(s) matching "{match_string}"')
    print()

    total_files = 0
    total_lines = 0
    total_corrected = 0
    total_unexpected = 0

    for filename in input_files:

        input_path = os.path.join(directory, filename)

        # Convert:
        #
        #   filelist_T5Zg_...
        #
        # to:
        #
        #   filelist_lpc_T5Zg_...
        #
        if filename.startswith("filelist_"):
            output_filename = (
                "filelist_lpc_" + filename[len("filelist_"):]
            )
        else:
            output_filename = "lpc_" + filename

        output_path = os.path.join(directory, output_filename)

        n_lines = 0
        n_corrected = 0
        n_unexpected = 0

        with open(input_path, "r") as infile, \
             open(output_path, "w") as outfile:

            for line in infile:
                line = line.strip()

                if not line:
                    continue

                n_lines += 1

                if line.startswith(prefix_to_remove):
                    line = line[len(prefix_to_remove):]
                    n_corrected += 1
                else:
                    print(
                        f"WARNING [{filename}]: "
                        f"unexpected line format:\n  {line}"
                    )
                    n_unexpected += 1

                outfile.write(line + "\n")

        print(f"Input:   {input_path}")
        print(f"Output:  {output_path}")
        print(f"  Lines:       {n_lines}")
        print(f"  Corrected:   {n_corrected}")
        print(f"  Unexpected:  {n_unexpected}")
        print()

        total_files += 1
        total_lines += n_lines
        total_corrected += n_corrected
        total_unexpected += n_unexpected

    print("--------------------------------------------------")
    print(f"Files processed:    {total_files}")
    print(f"Total lines:        {total_lines}")
    print(f"Total corrected:    {total_corrected}")
    print(f"Total unexpected:   {total_unexpected}")


if __name__ == "__main__":
    main()
