#!/usr/bin/env python3
import json
import sys

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 json_to_filelist.py <input.json> <output.txt> [--xrootd]")
        print("Example:")
        print("  python3 json_to_filelist.py dataset.json filelist.txt --xrootd")
        sys.exit(1)

    input_json = sys.argv[1]
    output_txt = sys.argv[2]
    use_xrootd = ("--xrootd" in sys.argv) or ("-x" in sys.argv)

    XROOTD_PREFIX = "root://cmseos.fnal.gov/"

    # --- Read JSON file ---
    with open(input_json, "r") as f:
        try:
            data = json.load(f)
        except json.JSONDecodeError as e:
            print(f"❌ Error parsing JSON: {e}")
            sys.exit(1)

    if not isinstance(data, dict):
        print("❌ Error: JSON root must be an object (dict of lists).")
        sys.exit(1)

    # --- Collect all file paths ---
    all_files = []
    for key, value in data.items():
        if isinstance(value, list):
            all_files.extend(value)
        else:
            print(f"⚠️ Warning: key '{key}' does not contain a list, skipping.")

    if not all_files:
        print("❌ No file paths found in JSON.")
        sys.exit(1)

    # --- Write to output ---
    with open(output_txt, "w") as out:
        for path in all_files:
            if use_xrootd:
                out.write(f"{XROOTD_PREFIX}{path}\n")
            else:
                out.write(f"{path}\n")

    print(f"✅ Wrote {len(all_files)} file paths to '{output_txt}'")
    if use_xrootd:
        print(f"   (with prefix '{XROOTD_PREFIX}')")

if __name__ == "__main__":
    main()

