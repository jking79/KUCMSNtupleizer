#!/usr/bin/env python3

import argparse
import os
import re
import shlex
import subprocess
import sys
from typing import List


def run_cmd(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command and return the CompletedProcess."""
    return subprocess.run(cmd, text=True, capture_output=True, check=check)


def eos_find_files(eos_path: str, recursive: bool) -> List[str]:
    """
    Return a list of file paths on EOS.
    Uses `eos <root> find -f <path>` when recursive,
    otherwise uses `eos <root> ls <path>` and filters to files by convention.
    """
    eos_path = eos_path.rstrip("/")

    if recursive:
        cmd = ["eos", "root://cmseos.fnal.gov", "find", "-f", eos_path]
        proc = run_cmd(cmd)
        files = [line.strip() for line in proc.stdout.splitlines() if line.strip()]
        return files

    cmd = ["eos", "root://cmseos.fnal.gov", "ls", eos_path]
    proc = run_cmd(cmd)

    files = []
    for entry in proc.stdout.splitlines():
        entry = entry.strip()
        if not entry:
            continue
        full = f"{eos_path}/{entry}"
        files.append(full)
    return files


def rename_path(old_path: str, search: str, replace: str, use_regex: bool) -> str:
    """Return the renamed path, changing only the basename."""
    directory = os.path.dirname(old_path)
    basename = os.path.basename(old_path)

    if use_regex:
        new_basename = re.sub(search, replace, basename)
    else:
        new_basename = basename.replace(search, replace)

    return os.path.join(directory, new_basename)


def eos_exists(path: str) -> bool:
    """Check whether a path exists on EOS."""
    cmd = ["eos", "root://cmseos.fnal.gov", "ls", path]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    return proc.returncode == 0


def eos_mv(old_path: str, new_path: str, dry_run: bool) -> bool:
    """Rename a file on EOS using eos mv."""
    cmd = ["eos", "root://cmseos.fnal.gov", "mv", old_path, new_path]

    if dry_run:
        print("[DRY-RUN]", " ".join(shlex.quote(x) for x in cmd))
        return True

    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.returncode != 0:
        print(f"[ERROR] Failed to rename:\n  {old_path}\n  -> {new_path}", file=sys.stderr)
        if proc.stderr:
            print(proc.stderr.strip(), file=sys.stderr)
        return False

    print(f"[RENAMED] {old_path} -> {new_path}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Rename EOS files using string or regex search/replace."
    )
    parser.add_argument(
        "eos_path",
        help="EOS directory to scan, e.g. /store/user/kingj/mydir"
    )
    parser.add_argument(
        "search",
        help="Search string or regex pattern"
    )
    parser.add_argument(
        "replace",
        help="Replacement string"
    )
    parser.add_argument(
        "--regex",
        action="store_true",
        help="Treat search as a regex pattern"
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search recursively with `eos find -f`"
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually perform the rename. Without this, script does a dry run."
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip rename if destination already exists"
    )
    parser.add_argument(
        "--match-full-path",
        action="store_true",
        help="Apply replacement to the full path instead of only the basename"
    )

    args = parser.parse_args()
    dry_run = not args.execute

    try:
        files = eos_find_files(args.eos_path, args.recursive)
    except subprocess.CalledProcessError as e:
        print("[ERROR] Failed to list EOS files.", file=sys.stderr)
        if e.stderr:
            print(e.stderr.strip(), file=sys.stderr)
        sys.exit(1)

    if not files:
        print("No files found.")
        sys.exit(0)

    rename_pairs = []

    for old_path in files:
        if args.match_full_path:
            if args.regex:
                new_path = re.sub(args.search, args.replace, old_path)
            else:
                new_path = old_path.replace(args.search, args.replace)
        else:
            new_path = rename_path(old_path, args.search, args.replace, args.regex)

        if new_path == old_path:
            continue

        rename_pairs.append((old_path, new_path))

    if not rename_pairs:
        print("No matching files found for rename.")
        sys.exit(0)

    print(f"Found {len(rename_pairs)} file(s) to rename:\n")
    for old_path, new_path in rename_pairs:
        print(f"  {old_path}")
        print(f"    -> {new_path}")

    print()

    success = 0
    skipped = 0
    failed = 0

    for old_path, new_path in rename_pairs:
        if args.skip_existing and eos_exists(new_path):
            print(f"[SKIP] Destination exists: {new_path}")
            skipped += 1
            continue

        ok = eos_mv(old_path, new_path, dry_run=dry_run)
        if ok:
            success += 1
        else:
            failed += 1

    print("\nSummary:")
    print(f"  Successful: {success}")
    print(f"  Skipped:    {skipped}")
    print(f"  Failed:     {failed}")
    print(f"  Dry-run:    {dry_run}")


if __name__ == "__main__":
    main()
