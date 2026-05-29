#!/usr/bin/env python3

import subprocess
import sys
import os
import re

EOS_MGM = "root://cmseos.fnal.gov"


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
    """
    Recursively find all .root files under 'path' on EOS.
    Returns a list of absolute EOS paths as returned by eos find.
    Often /eos/uscms/store/... or /store/...
    """
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
        if p and p.endswith(".root"):
            files.append(p)

    return files


def normalize_store_path(p: str) -> str:
    """
    Normalize EOS physical mount paths to a canonical /store/... form.
    Also collapses duplicate slashes.

    Examples:
      /eos/uscms/store/user/...  -> /store/user/...
      /store//user/...           -> /store/user/...
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


def safe_filename(s: str) -> str:
    """
    Make a directory name safe to use in an output filename.
    Mostly keeps normal dataset-name characters, but replaces anything weird.
    """
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", s)


def main():
    if len(sys.argv) > 1:
        start_dir = sys.argv[1]
    else:
        start_dir = "/store/user/lpcsusylep/jaking/KUCMSNtuple/"

    if len(sys.argv) > 2:
        campaign_pattern = sys.argv[2]
    else:
        campaign_pattern = ""

    if len(sys.argv) > 3:
        ht_pattern = sys.argv[3]
    else:
        ht_pattern = ""

    start_dir = start_dir.rstrip("/")

    print(f"[INFO] Scanning first-level subdirectories under: {start_dir}")

    campaigns = eos_ls(start_dir)

    matching_campaigns = [d for d in campaigns if campaign_pattern in d]

    if not matching_campaigns:
        print("[WARN] No matching top-level subdirectories found.")
        return

    for campaign in matching_campaigns:
        campaign_dir = f"{start_dir}/{campaign}"

        print(f"[INFO] Scanning HT-slice subdirectories under: {campaign_dir}")

        ht_slices = eos_ls(campaign_dir)

        matching_ht_slices = [d for d in ht_slices if ht_pattern in d]

        if not matching_ht_slices:
            print(f"[WARN]   No matching HT-slice subdirectories found under {campaign_dir}")
            continue

        for ht_slice in matching_ht_slices:
            ht_dir = f"{campaign_dir}/{ht_slice}"

            print(f"[INFO] Finding .root files under: {ht_dir}")

            root_files = eos_find_rootfiles(ht_dir)

            if not root_files:
                print(f"[INFO]   No .root files found under {ht_dir}, skipping.")
                continue

            outfile = f"filelist_{safe_filename(campaign)}__{safe_filename(ht_slice)}.txt"

            print(f"[INFO]   Writing {len(root_files)} files to {outfile}")

            with open(outfile, "w") as f:
                for rf in root_files:
                    f.write(relative_to_start(rf, start_dir) + "\n")


if __name__ == "__main__":
    main()
