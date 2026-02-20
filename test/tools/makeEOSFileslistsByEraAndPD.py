#!/usr/bin/env python3

import subprocess
import sys
import os
import re
from collections import defaultdict

EOS_MGM = "root://cmseos.fnal.gov"
RUN_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def eos_ls(path):
    cmd = f"eos {EOS_MGM} ls {path}"
    result = subprocess.run(
        cmd, shell=True, check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def eos_find_rootfiles(path):
    cmd = f"eos {EOS_MGM} find -f {path}"
    result = subprocess.run(
        cmd, shell=True, check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip().endswith(".root")]


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
    Finds 'Run2024B', 'Run2024C', etc as substrings.
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


def safe_name(s: str) -> str:
    """Make a filesystem-safe token."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s)


def main():
    # Inputs
    start_dir = sys.argv[1] if len(sys.argv) > 1 else "/store/user/lpcsusylep/jaking/KUCMSNtuple/"
    dir_pattern = sys.argv[2] if len(sys.argv) > 2 else ""     # filter first-level dirs (PD shards)
    year = sys.argv[3] if len(sys.argv) > 3 else "2024"

    start_dir = start_dir.rstrip("/")

    # output prefix: use last path element if no explicit pattern, otherwise pattern
    out_prefix = safe_name(dir_pattern) if dir_pattern else safe_name(os.path.basename(start_dir))

    print(f"[INFO] Scanning first-level subdirectories under: {start_dir}")
    if dir_pattern:
        print(f"[INFO] Filtering first-level dirs by pattern: '{dir_pattern}'")
    print(f"[INFO] Grouping by PD shard (first-level dir) AND run era tokens: Run{year}A..Z")

    # List first-level entries and filter to matching subdirs
    entries = eos_ls(start_dir)
    subdirs = [d for d in entries if (dir_pattern in d)] if dir_pattern else entries

    if not subdirs:
        print("[WARN] No matching subdirectories found.")
        return

    # grouped[pd_shard][era] = [rel_paths...]
    grouped = defaultdict(lambda: defaultdict(list))

    for sub in subdirs:
        pd_shard = sub  # e.g. EGamma0
        full_subdir = f"{start_dir}/{sub}"
        print(f"[INFO] Finding .root files under PD shard: {pd_shard}  ({full_subdir})")

        root_files = eos_find_rootfiles(full_subdir)
        if not root_files:
            print(f"[INFO]   No .root files found under {full_subdir}, skipping.")
            continue

        print(f"[INFO]   Found {len(root_files)} .root files under {full_subdir}")

        # DEBUG sample
        print("[DEBUG] Sample paths (relative):")
        for samp in root_files[:3]:
            rel_samp = relative_to_start(samp, start_dir)
            print("   ", rel_samp)

        for rf in root_files:
            rel = relative_to_start(rf, start_dir)
            era = extract_run_era_by_string(rel, year=year)
            grouped[pd_shard][era].append(rel)

    if not grouped:
        print("[WARN] No .root files found anywhere.")
        return

    print("\n[INFO] Writing grouped output files:")
    total_files_written = 0
    for pd_shard in sorted(grouped.keys()):
        for era in sorted(grouped[pd_shard].keys()):
            rel_paths = uniq_preserve_order(grouped[pd_shard][era])
            outname = f"{out_prefix}_{safe_name(pd_shard)}_{era}.txt"
            print(f"  - {outname}: {len(rel_paths)} files")
            with open(outname, "w", encoding="utf-8") as f:
                for p in rel_paths:
                    f.write(p + "\n")
            total_files_written += len(rel_paths)

    print("\n[SUMMARY] (pre-unique counts)")
    for pd_shard in sorted(grouped.keys()):
        for era in sorted(grouped[pd_shard].keys()):
            print(f"  {pd_shard:12s}  {era:9s}: {len(grouped[pd_shard][era])} files")

    print(f"\n[SUMMARY] Total unique-written (sum of per-filelist uniques): {total_files_written}")


if __name__ == "__main__":
    main()
