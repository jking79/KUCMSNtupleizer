#!/usr/bin/env python3

import subprocess
import sys

EOS_MGM = "root://cmseos.fnal.gov"
EOS_PREFIX = f"{EOS_MGM}/"


def eos_ls(path):
    """
    List entries (files + dirs) in an EOS directory.
    Returns a list of names (no leading path).
    """
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
    Returns a list of full xrootd paths (root://cmseos.fnal.gov/…).
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
        if not p:
            continue
        if p.endswith(".root"):
            # make sure there’s exactly one leading slash after EOS_PREFIX
            files.append(EOS_PREFIX + '/' + p.replace("/eos/uscms/", "/").lstrip("/"))

    return files


def main():
    # ------------------------------------------------------------------
    # Start directory: either from command line or hard-coded default
    # ------------------------------------------------------------------
    if len(sys.argv) > 1:
        start_dir = sys.argv[1]
    else:
        # default example – change this to whatever you like
        start_dir = "/store/user/lpcsusylep/jaking/KUCMSNtuple/"
    if len(sys.argv) > 2:
        dir_pattern = sys.argv[2]
    else:
        dir_pattern = ""

    start_dir = start_dir.rstrip("/")

    print(f"[INFO] Scanning first-level subdirectories under: {start_dir}")

    # Get first-level subdirectories
    entries = eos_ls(start_dir)

    # Optional: filter entries if you only want certain names
    # e.g. dir_pattern = "2500_mN2-1500_mN1-1000"
    subdirs = [
        d for d in entries
        if (dir_pattern in d)  # keep this if you want pattern matching
    ]

    if not subdirs:
        print("[WARN] No matching subdirectories found.")
        return

    for sub in subdirs:
        full_subdir = f"{start_dir}/{sub}"
        print(f"[INFO] Finding .root files under: {full_subdir}")

        root_files = eos_find_rootfiles(full_subdir)

        if not root_files:
            print(f"[INFO]   No .root files found under {full_subdir}, skipping.")
            continue

        # One output file per first-level subdirectory
        outfile = f"filelist_{sub}.txt"
        print(f"[INFO]   Writing {len(root_files)} files to {outfile}")

        with open(outfile, "w") as f:
            for rf in root_files:
                f.write(rf + "\n")


if __name__ == "__main__":
    main()

