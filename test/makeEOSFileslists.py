#!/usr/bin/env python3

import subprocess
import sys
import os

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
    (Often /eos/uscms/store/... or /store/...)
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
      /store//user/...          -> /store/user/...
    """
    p = p.strip()

    # Convert common EOS mount prefixes to /store
    for prefix in ("/eos/uscms/store", "/eos/cms/store"):
        if p.startswith(prefix):
            p = "/store" + p[len(prefix):]
            break

    # Collapse duplicate slashes
    p = os.path.normpath(p)

    # normpath removes trailing slash; that's fine for file paths
    return p


def relative_to_start(full_path: str, start_dir: str) -> str:
    """
    Make full_path relative to start_dir, robust to /eos/uscms/store vs /store.
    """
    full_n = normalize_store_path(full_path)
    start_n = normalize_store_path(start_dir)

    # Ensure start_n is treated like a directory prefix
    start_prefix = start_n.rstrip("/") + "/"

    if full_n.startswith(start_prefix):
        return full_n[len(start_prefix):]

    # Fallback: if something doesn't match, return normalized full path
    return full_n


def main():
    if len(sys.argv) > 1:
        start_dir = sys.argv[1]
    else:
        start_dir = "/store/user/lpcsusylep/jaking/KUCMSNtuple/"

    start_dir = start_dir.rstrip("/")

    print(f"[INFO] Scanning first-level subdirectories under: {start_dir}")

    entries = eos_ls(start_dir)

    dir_pattern = ""
    subdirs = [d for d in entries if dir_pattern in d]

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

        outfile = f"filelist_{sub}.txt"
        print(f"[INFO]   Writing {len(root_files)} files to {outfile}")

        with open(outfile, "w") as f:
            for rf in root_files:
                f.write(relative_to_start(rf, start_dir) + "\n")


if __name__ == "__main__":
    main()

