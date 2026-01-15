#!/usr/bin/env python3
import sys
from pathlib import Path

EOS_PREFIX = "root://cmseos.fnal.gov/"

def add_prefix(line: str) -> str:
    """
    Strip whitespace/newlines, and if the line is non-empty and not already prefixed,
    add EOS_PREFIX to the front.
    """
    s = line.strip()
    if not s:
        return ""
    if s.startswith("root://"):
        # already a full xrootd URI
        return s
    # If input begins with /store/... then this makes:
    # root://cmseos.fnal.gov//store/...
    return EOS_PREFIX + s

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} INPUT.txt OUTPUT.txt")
        sys.exit(1)

    infile = Path(sys.argv[1])
    outfile = Path(sys.argv[2])

    if not infile.exists():
        print(f"ERROR: Input file not found: {infile}")
        sys.exit(2)

    n_in = 0
    n_out = 0

    with infile.open("r", encoding="utf-8") as fin, outfile.open("w", encoding="utf-8") as fout:
        for line in fin:
            n_in += 1
            new_line = add_prefix(line)
            if new_line:
                fout.write(new_line + "\n")
                n_out += 1

    print(f"Done.")
    print(f"  Read lines:    {n_in}")
    print(f"  Wrote entries: {n_out}")
    print(f"  Output file:   {outfile}")

if __name__ == "__main__":
    main()

