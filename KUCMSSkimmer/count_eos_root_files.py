#!/usr/bin/env python3

import argparse
import subprocess
import sys


def run_command(cmd):
    """Run a shell command and return stdout lines."""
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except Exception as err:
        print(f"ERROR: failed to run command: {' '.join(cmd)}", file=sys.stderr)
        print(f"ERROR: {err}", file=sys.stderr)
        return []

    if result.returncode != 0:
        print(f"WARNING: command returned nonzero status: {' '.join(cmd)}", file=sys.stderr)
        print(result.stderr.strip(), file=sys.stderr)

    lines = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if line:
            lines.append(line)

    return lines


def count_root_files(eos_dir, eos_host="cmseos.fnal.gov", match_string=""):
    """Count ROOT files under eos_dir using eos find -f."""

    cmd = [
        "eos",
        f"root://{eos_host}",
        "find",
        "-f",
        eos_dir,
    ]

    print("Running command:")
    print("  " + " ".join(cmd))

    lines = run_command(cmd)

    n_lines = 0
    n_root = 0
    n_match_fail = 0
    root_files = []

    for path in lines:
        n_lines += 1

        if not path.endswith(".root"):
            continue

        n_root += 1

        if match_string and match_string not in path:
            n_match_fail += 1
            continue

        root_files.append(path)

    print()
    print("EOS ROOT file count summary")
    print("---------------------------")
    print(f"eos_host     : {eos_host}")
    print(f"eos_dir      : {eos_dir}")
    print(f"match_string : {match_string}")
    print(f"raw lines    : {n_lines}")
    print(f"root files   : {n_root}")
    print(f"match fail   : {n_match_fail}")
    print(f"accepted     : {len(root_files)}")

    return root_files


def main():
    parser = argparse.ArgumentParser(
        description="Count ROOT files recursively under an EOS directory."
    )

    parser.add_argument(
        "eos_dir",
        help="EOS directory to scan, for example /store/user/lpcsusylep/jaking/KUCMSNtuple/...",
    )

    parser.add_argument(
        "--eos-host",
        default="cmseos.fnal.gov",
        help="EOS host. Default: cmseos.fnal.gov",
    )

    parser.add_argument(
        "--match",
        default="",
        help="Optional substring required in the full EOS path.",
    )

    parser.add_argument(
        "--write-list",
        default="",
        help="Optional output text file to write accepted ROOT file paths.",
    )

    parser.add_argument(
        "--xrootd",
        action="store_true",
        help="Write/list paths as root:// URLs instead of /store/... paths.",
    )

    args = parser.parse_args()

    root_files = count_root_files(
        eos_dir=args.eos_dir,
        eos_host=args.eos_host,
        match_string=args.match,
    )

    if args.xrootd:
        root_files = [
            f"root://{args.eos_host}/{path}" if path.startswith("/store/") else path
            for path in root_files
        ]

    if args.write_list:
        with open(args.write_list, "w") as fout:
            for path in root_files:
                fout.write(path + "\n")

        print()
        print(f"Wrote accepted ROOT file list to: {args.write_list}")


if __name__ == "__main__":
    main()
