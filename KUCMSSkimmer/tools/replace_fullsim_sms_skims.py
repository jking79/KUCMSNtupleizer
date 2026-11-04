#!/usr/bin/env python3
"""Copy local v52 full-simulation SMS skims to their EOS v52 names.

Run from KUCMSSkimmer on LPC after all rows in the v34 SMS list have finished:

    python3 tools/replace_fullsim_sms_skims.py --source-dir .

The default is a preview. Add --execute to upload; existing target files are
overwritten only after an exact confirmation prompt. Local files are retained.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


EOS_HOST = "root://cmseos.fnal.gov"
DEFAULT_EOS_DIR = "/store/group/lpcsusylep/malazaro/KUCMSSkims/skims_v53"
DEFAULT_INPUT_LIST = "backups/KUCMS_Ntuple_Master_v34_SMS_Sig_Files_List.txt"
LOCAL_SUFFIX = "_rjrskim_v52.root"
EOS_SUFFIX = "_rjrskim_SVHPM100_Full_v52.root"


def target_basename(sample: str) -> str:
    """Translate an input-list sample name to the legacy EOS rjrskim name."""
    if sample.startswith("GlGl_"):
        # This matches the existing EOS v52 names, e.g.
        # GlGl_mGl-2300_..._GZ_N2ctau-0p1 ->
        # gogoGZ_mGl-2300_..._GZ_ct0p1_rjrskim_SVHPM100_Full_v52.root
        eos_sample = sample.removeprefix("GlGl_").replace(
            "_GZ_N2ctau-", "_GZ_ct", 1
        )
        if eos_sample == sample.removeprefix("GlGl_"):
            raise ValueError(f"GlGl sample lacks expected _GZ_N2ctau- token: {sample}")
        return f"gogoGZ_{eos_sample}{EOS_SUFFIX}"

    if sample.startswith("T5Zg_"):
        # No T5Zg v52 rjrskim example was present in the supplied EOS listing;
        # use the confirmed legacy rjrskim suffix while preserving the sample ID.
        return f"{sample}{EOS_SUFFIX}"

    raise ValueError(f"Unsupported full-simulation SMS sample name: {sample}")


def read_mapping(input_list: Path) -> list[tuple[str, str]]:
    mapping: list[tuple[str, str]] = []
    seen_sources: set[str] = set()
    seen_targets: set[str] = set()
    for line_number, raw_line in enumerate(input_list.read_text().splitlines(), 1):
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) < 2:
            raise ValueError(f"Malformed row {line_number} in {input_list}: {raw_line}")
        sample = fields[1]
        source = f"{sample}{LOCAL_SUFFIX}"
        target = target_basename(sample)
        if source in seen_sources:
            raise ValueError(f"Duplicate local source in input list: {source}")
        if target in seen_targets:
            raise ValueError(f"Multiple samples map to EOS target: {target}")
        seen_sources.add(source)
        seen_targets.add(target)
        mapping.append((source, target))
    if not mapping:
        raise ValueError(f"No active rows found in {input_list}")
    return mapping


def run(command: list[str], *, capture: bool = False) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-dir", type=Path, default=Path("."),
        help="directory containing the local *_rjrskim_v52.root outputs (default: .)",
    )
    parser.add_argument(
        "--input-list", type=Path, default=Path(DEFAULT_INPUT_LIST),
        help=f"SMS master list relative to the working directory (default: {DEFAULT_INPUT_LIST})",
    )
    parser.add_argument(
        "--eos-dir", default=DEFAULT_EOS_DIR,
        help=f"EOS destination directory (default: {DEFAULT_EOS_DIR})",
    )
    parser.add_argument(
        "--execute", action="store_true",
        help="upload files and overwrite matching EOS targets after confirmation",
    )
    args = parser.parse_args()

    try:
        mapping = read_mapping(args.input_list)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    print(f"Input list: {args.input_list}")
    print(f"Local source directory: {args.source_dir}")
    print(f"EOS destination: {args.eos_dir}")
    print(f"Mapped full-sim SMS files: {len(mapping)}")

    missing: list[str] = []
    for source_name, target_name in mapping:
        source = args.source_dir / source_name
        if not source.is_file() or source.stat().st_size == 0:
            missing.append(source_name)
            state = "MISSING/EMPTY"
        else:
            state = f"{source.stat().st_size} bytes"
        print(f"{state}: {source_name}  ->  {target_name}")

    if missing:
        print(f"ERROR: {len(missing)} expected local output(s) are missing or empty.", file=sys.stderr)
        if args.execute:
            print("No files uploaded.", file=sys.stderr)
            return 2

    if not args.execute:
        print("PREVIEW ONLY: no EOS files were changed. Re-run with --execute after reviewing the map.")
        return 0 if not missing else 2

    for executable in ("xrdfs", "xrdcp"):
        if not shutil.which(executable):
            print(f"ERROR: required command not found: {executable}", file=sys.stderr)
            return 2

    listing = run(["xrdfs", EOS_HOST, "ls", args.eos_dir], capture=True)
    if listing.returncode != 0:
        print("ERROR: could not list the EOS destination; no files uploaded.", file=sys.stderr)
        if listing.stderr:
            print(listing.stderr.rstrip(), file=sys.stderr)
        return 2
    existing = {Path(line.rstrip("/")).name for line in listing.stdout.splitlines()}
    overwrite_count = sum(target in existing for _, target in mapping)
    print(f"EOS targets already present and to be overwritten: {overwrite_count}")

    confirmation = f"UPLOAD {len(mapping)} FULLSIM SMS FILES"
    print(f"This will upload {len(mapping)} files; {overwrite_count} existing target(s) will be replaced.")
    if input(f'Type exactly "{confirmation}" to continue: ').strip() != confirmation:
        print("Cancelled; no files uploaded.")
        return 1

    failures: list[str] = []
    for source_name, target_name in mapping:
        source = args.source_dir / source_name
        destination = f"{EOS_HOST}//{args.eos_dir.strip('/')}/{target_name}"
        print(f"Uploading {source} -> {destination}", flush=True)
        result = run(["xrdcp", "-f", str(source), destination])
        if result.returncode != 0:
            failures.append(source_name)
            print(f"ERROR: upload failed for {source_name}", file=sys.stderr)
    if failures:
        print(f"Completed with {len(failures)} upload failure(s):", file=sys.stderr)
        for name in failures:
            print(f"  {name}", file=sys.stderr)
        return 1

    print(f"Successfully uploaded {len(mapping)} full-sim SMS files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
