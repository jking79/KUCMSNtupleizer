#!/usr/bin/env python3
"""Build one EOS ROOT-file list per T5Zg MiniAOD mass point.

Expected EOS layout:

  <start-dir>/T5Zg_FullSim_Mini_<year>/
    T5Zg_mGo-..._<year>/MINI/<production timestamp>/<slice>/*.root

All ROOT files below each mass point's MINI directory are placed in a
single list.  Other tiers beside MINI (AOD, GEN, HLT, SIM, digi, ...) are
never scanned.
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path


EOS_MGM = "root://cmseos.fnal.gov"
DEFAULT_START_DIR = "/store/user/lpcsusylep/jaking/KUCMSNtuple"


def run_eos(*arguments):
    """Run an EOS command and return its non-empty output lines."""
    command = ["eos", EOS_MGM, *arguments]
    try:
        result = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except FileNotFoundError:
        raise RuntimeError(
            "The 'eos' command was not found. Run this on cmslpc after "
            "setting up the EOS client environment."
        ) from None
    except subprocess.CalledProcessError as error:
        detail = error.stderr.strip() or error.stdout.strip()
        raise RuntimeError(
            f"EOS command failed ({' '.join(command)}): {detail}"
        ) from error

    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def eos_ls(path):
    return run_eos("ls", path)


def eos_find_rootfiles(path):
    """Recursively find ROOT files below one EOS directory."""
    return sorted(
        {
            line
            for line in run_eos("find", "-f", path)
            if line.endswith(".root")
        }
    )


def normalize_store_path(path):
    """Normalize EOS mount paths to the canonical /store/... form."""
    path = path.strip()
    for prefix in ("/eos/uscms/store", "/eos/cms/store"):
        if path.startswith(prefix):
            path = "/store" + path[len(prefix):]
            break
    return os.path.normpath(path)


def xrootd_url(path):
    """Convert an EOS path to a root://cmseos.fnal.gov//store/... URL."""
    normalized = normalize_store_path(path)
    return f"{EOS_MGM}/{normalized}"


def safe_filename(value):
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", value)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "start_dir",
        nargs="?",
        default=DEFAULT_START_DIR,
        help=f"EOS directory containing the campaigns (default: {DEFAULT_START_DIR})",
    )
    parser.add_argument(
        "--campaign-pattern",
        default="T5Zg_FullSim_Mini_",
        help="substring required in a campaign directory name",
    )
    parser.add_argument(
        "--mass-point-prefix",
        "--sample-pattern",
        dest="mass_point_prefix",
        default="T5Zg_mGo-",
        help="prefix required on a mass-point directory name",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="local output directory for file lists (default: current directory)",
    )
    return parser.parse_args(argv)


def build_filelists(args):
    start_dir = normalize_store_path(args.start_dir).rstrip("/")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Scanning campaigns under: {start_dir}")
    campaigns = sorted(
        name for name in eos_ls(start_dir) if args.campaign_pattern in name
    )
    if not campaigns:
        print(
            f"[WARN] No campaigns containing {args.campaign_pattern!r} "
            f"were found under {start_dir}"
        )
        return 0

    lists_written = 0
    files_written = 0

    for campaign in campaigns:
        campaign_dir = f"{start_dir}/{campaign}"
        print(f"[INFO] Scanning mass points under: {campaign_dir}")

        mass_points = sorted(
            name
            for name in eos_ls(campaign_dir)
            if name.startswith(args.mass_point_prefix)
        )
        if not mass_points:
            print(f"[WARN]   No {args.mass_point_prefix}* directories found")
            continue

        for mass_point in mass_points:
            mass_point_dir = f"{campaign_dir}/{mass_point}"
            mini_dir = f"{mass_point_dir}/MINI"

            try:
                tiers = eos_ls(mass_point_dir)
            except RuntimeError as error:
                print(f"[WARN]   {error}", file=sys.stderr)
                continue

            if "MINI" not in tiers:
                print(f"[WARN]   Missing MINI directory: {mass_point_dir}")
                continue

            print(f"[INFO] Finding .root files under: {mini_dir}")
            try:
                root_files = eos_find_rootfiles(mini_dir)
            except RuntimeError as error:
                print(f"[WARN]   {error}", file=sys.stderr)
                continue

            if not root_files:
                print(f"[WARN]   No .root files found under {mini_dir}")
                continue

            output_name = (
                f"filelist_{safe_filename(campaign)}__"
                f"{safe_filename(mass_point)}_MINI.txt"
            )
            output_path = output_dir / output_name
            output_path.write_text(
                "".join(f"{xrootd_url(path)}\n" for path in root_files),
                encoding="utf-8",
            )

            print(f"[INFO]   Wrote {len(root_files)} files to {output_path}")
            lists_written += 1
            files_written += len(root_files)

    print(
        f"[INFO] Done: wrote {lists_written} list(s) containing "
        f"{files_written} ROOT file(s)"
    )
    return lists_written


def main(argv=None):
    args = parse_args(argv)
    try:
        build_filelists(args)
    except RuntimeError as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
