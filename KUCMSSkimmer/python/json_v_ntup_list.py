#!/usr/bin/env python3
import json
import sys
import re
import os
from datetime import datetime

def main(json_file, txt_file):
    # Load JSON
    with open(json_file, "r") as jf:
        data = json.load(jf)

    # Load text file of file paths
    with open(txt_file, "r") as tf:
        filenames = [line.strip() for line in tf if line.strip()]

    # Extract numeric indices from file names (e.g. kucmsntuple_123.root → "123")
    pattern = re.compile(r'_(\d+)\.root$')
    file_indices = {match.group(1) for name in filenames if (match := pattern.search(name))}

    # JSON keys as string set
    json_indices = set(data.keys())

    # Find missing indices (present in JSON but not in filenames)
    missing_indices = [i for i in json_indices if i not in file_indices]

    # Find skipped indices (missing numbers in JSON index sequence)
    json_ints = sorted(map(int, json_indices))
    skipped = []
    if len(json_ints) > 1:
        expected = set(range(json_ints[0], json_ints[-1] + 1))
        skipped = sorted(expected - set(json_ints))

    # Build timestamped report filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    txt_base = os.path.splitext(os.path.basename(txt_file))[0]
    report_name = f"missing_report_{txt_base}_{timestamp}.txt"

    # Write missing and skipped to report
    with open(report_name, "w") as out:
        out.write(f"Report generated using file list: {os.path.basename(txt_file)}\n")
        out.write(f"Timestamp: {timestamp}\n")
        out.write("=" * 60 + "\n\n")

        # Missing files section
        out.write("===== Missing Files (in JSON but not in file list) =====\n")
        if missing_indices:
            for idx in sorted(map(int, missing_indices)):
                paths = data.get(str(idx), [])
                for p in paths:
                    out.write(p + "\n")
        else:
            out.write("None\n")

        # Skipped indices section
        out.write("\n===== Skipped Indices (gaps in JSON sequence) =====\n")
        if skipped:
            for idx in skipped:
                str_idx = str(idx)
                paths = data.get(str_idx)
                if paths:
                    for p in paths:
                        out.write(p + "\n")
                else:
                    out.write(f"[No entry for index {idx} in JSON]\n")
        else:
            out.write("None\n")

    # Console summary
    print(f"✅ Report written to {report_name}")
    print(f"❌ Missing indices: {len(missing_indices)}")
    print(f"⚠️ Skipped indices: {len(skipped)}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python check_json_missing_indices.py <json_file> <file_list.txt>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])

