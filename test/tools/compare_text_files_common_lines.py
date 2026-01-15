#!/usr/bin/env python3
import sys
from pathlib import Path
from collections import Counter


def normalize_line(line: str, strip: bool = True) -> str:
    """
    Normalize a line for comparison.
    Default behavior:
      - strip whitespace/newline at both ends
    """
    return line.strip() if strip else line.rstrip("\n")


def load_lines(path: Path, ignore_blank: bool = True) -> list[str]:
    """
    Read lines from file and normalize them.
    Returns a list (keeps duplicates).
    """
    lines = []
    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            s = normalize_line(raw)
            if ignore_blank and not s:
                continue
            lines.append(s)
    return lines


def main():
    if len(sys.argv) not in (3, 4):
        print(f"Usage: {sys.argv[0]} file1.txt file2.txt [shared_lines_output.txt]")
        sys.exit(1)

    file1 = Path(sys.argv[1])
    file2 = Path(sys.argv[2])
    out_shared = Path(sys.argv[3]) if len(sys.argv) == 4 else None

    if not file1.exists():
        print(f"ERROR: file not found: {file1}")
        sys.exit(2)
    if not file2.exists():
        print(f"ERROR: file not found: {file2}")
        sys.exit(2)

    # Load lines (normalized)
    lines1 = load_lines(file1)
    lines2 = load_lines(file2)

    c1 = Counter(lines1)
    c2 = Counter(lines2)

    # Intersection: items that appear in BOTH files
    shared = sorted(set(c1.keys()) & set(c2.keys()))

    print("========================================")
    print("Compare 2 text files: shared line check")
    print("========================================")
    print(f"File 1: {file1}")
    print(f"File 2: {file2}")
    print("")
    print(f"Normalized non-blank lines in file1: {len(lines1)}  (unique: {len(c1)})")
    print(f"Normalized non-blank lines in file2: {len(lines2)}  (unique: {len(c2)})")
    print("")
    print(f"Shared unique lines: {len(shared)}")

    if len(shared) == 0:
        print("\n✅ No identical lines found between the two files.")
        return

    print("\n⚠️ Identical lines found in BOTH files:\n")

    # Print shared lines + how many times they appear in each file
    for i, line in enumerate(shared, start=1):
        print(f"[{i:4d}] (file1 x{c1[line]}) (file2 x{c2[line]})  {line}")

    # Optional output file
    if out_shared:
        with out_shared.open("w", encoding="utf-8") as fout:
            for line in shared:
                fout.write(line + "\n")
        print(f"\nWrote shared unique lines to: {out_shared}")


if __name__ == "__main__":
    main()

