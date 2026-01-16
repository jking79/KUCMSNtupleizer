#!/usr/bin/env python3

import subprocess
import sys
import os
from collections import defaultdict


EOS_MGM = "root://cmseos.fnal.gov"

# Letters expected in Run era; include all possible
RUN_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def eos_ls(path):
    cmd = f"eos {EOS_MGM} ls {path}"
    result = subprocess.run(
        cmd,
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def eos_find_rootfiles(path):
    cmd = f"eos {EOS_MGM} find -f {path}"
    result = subprocess.run(
        cmd,
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    files = []
    for line in result.stdout.splitlines():
        p = line.strip()
        if p.endswith(".root"):
            files.append(p)
    return files


def normalize_store_path(p: str) -> str:
    """
    Normalize EOS physical mount paths to canonical /store/... form.
    Also collapses duplicate slashes.
    """
    p = p.strip()

    for prefix in ("/eos/uscms/store", "/eos/cms/store"):
        if p.startswith(prefix):
            p = "/store" + p[len(prefix):]
            break

    return os.path.normpath(p)


def relative_to_start(full_path: str, start_dir: str) -> str:
    """
    Make full_path relative to start_dir, robust to /eos/uscms/store vs /store.
    """
    full_n = normalize_store_path(full_path)
    start_n = normalize_store_path(start_dir)

    start_prefix = start_n.rstrip("/") + "/"
    if full_n.startswith(start_prefix):
        return full_n[len(start_prefix):]
    return full_n


def extract_run_era_by_string(path_str: str, year: str = "2024") -> str:
    """
    String-only era finding:
      finds 'Run2024B', 'Run2024C', etc as substrings.
    Returns 'Run2024X' or 'UNKNOWN'.
    """
    tag_prefix = "Run" + year
    for letter in RUN_LETTERS:
        token = tag_prefix + letter
        if token in path_str:
            return token
    return "UNKNOWN"


def uniq_preserve_order(items):
    seen = set()
    out = []
    for x in items:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def main():
    # Inputs
    start_dir = sys.argv[1] if len(sys.argv) > 1 else "/store/user/lpcsusylep/jaking/KUCMSNtuple/"
    dir_pattern = sys.argv[2] if len(sys.argv) > 2 else ""
    year = sys.argv[3] if len(sys.argv) > 3 else "2024"  # optional year override

    start_dir = start_dir.rstrip("/")

    print(f"[INFO] Scanning first-level subdirectories under: {start_dir}")
    if dir_pattern:
        print(f"[INFO] Filtering first-level dirs by pattern: '{dir_pattern}'")
    print(f"[INFO] Grouping by run era tokens: Run{year}A..Z (string match)")

    # List first-level entries and filter to matching subdirs
    entries = eos_ls(start_dir)
    subdirs = [d for d in entries if (dir_pattern in d)] if dir_pattern else entries

    if not subdirs:
        print("[WARN] No matching subdirectories found.")
        return

    grouped = defaultdict(list)

    for sub in subdirs:
        full_subdir = f"{start_dir}/{sub}"
        print(f"[INFO] Finding .root files under: {full_subdir}")

        root_files = eos_find_rootfiles(full_subdir)
        if not root_files:
            print(f"[INFO]   No .root files found under {full_subdir}, skipping.")
            continue

        print(f"[INFO]   Found {len(root_files)} .root files under {full_subdir}")

        # DEBUG: show sample paths
        print("[DEBUG] Sample paths (relative):")
        for samp in root_files[:3]:
            rel_samp = relative_to_start(samp, start_dir)
            print("   ", rel_samp)

        # Group
        for rf in root_files:
            rel = relative_to_start(rf, start_dir)
            era = extract_run_era_by_string(rel, year=year)
            grouped[era].append(rel)

    if not grouped:
        print("[WARN] No .root files found anywhere.")
        return

    out_prefix = dir_pattern if dir_pattern else "filelist"

    print("\n[INFO] Writing grouped output files:")
    for era in sorted(grouped.keys()):
        rel_paths = uniq_preserve_order(grouped[era])
        outname = f"{out_prefix}_{era}.txt"
        print(f"  - {outname}: {len(rel_paths)} files")

        with open(outname, "w", encoding="utf-8") as f:
            for p in rel_paths:
                f.write(p + "\n")

    print("\n[SUMMARY] (pre-unique counts)")
    for era in sorted(grouped.keys()):
        print(f"  {era}: {len(grouped[era])} files")


if __name__ == "__main__":
    main()

