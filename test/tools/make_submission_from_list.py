#!/usr/bin/env python3
import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 prepare_submission_list.py <missing_files.txt> [output.txt]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "submission_list.txt"

    prefix = "root://cmseos.fnal.gov/"

    with open(input_file, "r") as f:
        missing_files = [line.strip() for line in f if line.strip()]

    with open(output_file, "w") as out:
        for path in missing_files:
            # Skip if it already has a prefix (avoid double prefixing)
            if path.startswith("root://"):
                out.write(f"{path}\n")
            else:
                out.write(f"{prefix}{path}\n")

    print(f"✅ Wrote {len(missing_files)} entries to {output_file}")

if __name__ == "__main__":
    main()

