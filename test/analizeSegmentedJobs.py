#!/usr/bin/env python3
import json
import sys
from collections import defaultdict, Counter

def read_json_file(filename):
    """Read CMS-style JSON, returning both a flat file list and a map of file → keys."""
    with open(filename, "r") as f:
        data = json.load(f)

    if not isinstance(data, dict):
        print(f"❌ Error: {filename} is not a JSON object.")
        return [], {}

    files = []
    file_to_keys = defaultdict(list)

    for key, value in data.items():
        if isinstance(value, list):
            for path in value:
                files.append(path)
                file_to_keys[path].append(key)
        else:
            print(f"⚠️ Warning: key '{key}' in {filename} does not contain a list.")

    return files, file_to_keys

def read_filelist(filename):
    """Reads a plain text file list."""
    with open(filename, "r") as f:
        return [line.strip() for line in f if line.strip()]

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 compare_json_filelists_detailed.py <all_files.txt> <json1> <json2> <json3> <json4> ...")
        sys.exit(1)

    nEntries = len(sys.argv)
    json_files = sys.argv[2:nEntries]
    ref_file = sys.argv[1]

    # --- Collect info across all JSONs ---
    all_entries = []
    file_locations = defaultdict(list)  # file_path → list of (json_file, key)
    
    for jf in json_files:
        print(f"📖 Reading {jf} ...")
        file_list, file_to_keys = read_json_file(jf)
        for path, keys in file_to_keys.items():
            for k in keys:
                file_locations[path].append((jf, k))
        all_entries.extend(file_list)

    print(f"✅ Extracted {len(all_entries)} total entries from {len(json_files)} JSON files.")

    # --- Find duplicates (files appearing in more than one place) ---
    counts = Counter(all_entries)
    duplicates = {f: file_locations[f] for f, n in counts.items() if n > 1}

    # --- Compare with reference list ---
    all_possible = read_filelist(ref_file)
    print(f"📚 Reference list has {len(all_possible)} total files.")

    found_set = set(all_entries)
    ref_set = set(all_possible)

    missing = sorted(list(ref_set - found_set))
    extra = sorted(list(found_set - ref_set))

    # --- Summary ---
    print("\n================ Summary ================")
    print(f"Total files (unique):     {len(found_set)}")
    print(f"Total reference files:    {len(ref_set)}")
    print(f"Missing files:            {len(missing)}")
    print(f"Extra/unexpected files:   {len(extra)}")
    print(f"Duplicate files:          {len(duplicates)}")
    print("=========================================\n")

    # --- Detailed duplicate report ---
    if duplicates:
        print("⚠️  Duplicates found:")
        for fpath, locs in duplicates.items():
            print(f"\n→ {fpath}")
            for jf, key in locs:
                print(f"     • in {jf}  (key: {key})")

    # --- Missing and extras ---
    if missing:
        print("\n❌ Missing files:")
        for f in missing:
            print("   ", f)
    if extra:
        print("\n⚠️ Extra/unexpected files:")
        for f in extra:
            print("   ", f)

    # --- Write reports ---
    with open("missing_files.txt", "w") as out:
        out.writelines(f"{f}\n" for f in missing)
    with open("duplicates_detailed.txt", "w") as out:
        for fpath, locs in duplicates.items():
            out.write(f"{fpath}\n")
            for jf, key in locs:
                out.write(f"    {jf}  (key: {key})\n")
            out.write("\n")

    print("\n✅ Reports written to:")
    print("   missing_files.txt")
    print("   duplicates_detailed.txt")

if __name__ == "__main__":
    main()

