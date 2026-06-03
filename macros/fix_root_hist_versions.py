#!/usr/bin/env python3

import argparse
import glob
import os
import re
import sys

import ROOT


VERSION_RE = re.compile(r"_v\d+")

def get_file_version_tag(path):
    """
    Extract the final _v#### tag from the ROOT filename.
    """
    basename = os.path.basename(path)
    stem = basename
    if stem.endswith(".root"):
        stem = stem[:-5]

    matches = list(VERSION_RE.finditer(stem))
    if not matches:
        raise ValueError("Could not find version tag like _v0501 in filename: {}".format(basename))

    return matches[-1].group(0)


def replace_last_version_tag(name, new_version):
    """
    Replace the last _v#### tag in a histogram name.

    Returns:
      new_name, old_version

    If no version tag is found:
      original_name, None
    """
    matches = list(VERSION_RE.finditer(name))
    if not matches:
        return name, None

    last = matches[-1]
    old_version = last.group(0)

    new_name = name[:last.start()] + new_version + name[last.end():]
    return new_name, old_version


def is_histogram(obj):
    """
    True for TH1/TH2/TH3-derived objects.
    """
    return obj.InheritsFrom("TH1")

def process_directory(directory, file_version: str, dry_run: bool, overwrite_existing: bool) -> tuple[int, int]:
    """
    Recursively process histograms in a ROOT directory.

    Returns:
      n_changed, n_skipped
    """
    n_changed = 0
    n_skipped = 0

    # Important: make a static list of keys before writing new objects.
    # Otherwise newly written keys can appear during iteration.
    keys = [key for key in directory.GetListOfKeys()]

    for key in keys:
        obj = key.ReadObj()
        obj_name = obj.GetName()

        if obj.InheritsFrom("TDirectory"):
            subdir = directory.Get(obj_name)
            changed, skipped = process_directory(
                subdir,
                file_version,
                dry_run=dry_run,
                overwrite_existing=overwrite_existing,
            )
            n_changed += changed
            n_skipped += skipped
            continue

        if not is_histogram(obj):
            continue

        new_name, old_version = replace_last_version_tag(obj_name, file_version)

        if old_version is None:
            print(f"[skip:no version] {directory.GetPath()} :: {obj_name}")
            n_skipped += 1
            continue

        if old_version == file_version:
            continue

        existing_key = directory.GetKey(new_name)
        if existing_key and not overwrite_existing:
            print(
                f"[skip:exists] {directory.GetPath()} :: "
                f"{obj_name} -> {new_name}"
            )
            n_skipped += 1
            continue

        print(
            f"[rename clone] {directory.GetPath()} :: "
            f"{obj_name} -> {new_name}"
        )

        if not dry_run:
            directory.cd()

            cloned = obj.Clone(new_name)
            cloned.SetName(new_name)
            cloned.SetTitle(obj.GetTitle())

            write_option = ROOT.TObject.kOverwrite if overwrite_existing else 0
            cloned.Write(new_name, write_option)

            # Avoid PyROOT ownership / directory weirdness.
            cloned.SetDirectory(0)

        n_changed += 1

    return n_changed, n_skipped


def process_file(path: str, dry_run: bool, overwrite_existing: bool) -> tuple[int, int]:
    file_version = get_file_version_tag(path)

    print(f"\n=== {path}")
    print(f"File version tag: {file_version}")

    if dry_run:
        root_file = ROOT.TFile.Open(path, "READ")
    else:
        root_file = ROOT.TFile.Open(path, "UPDATE")

    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {path}")

    changed, skipped = process_directory(
        root_file,
        file_version,
        dry_run=dry_run,
        overwrite_existing=overwrite_existing,
    )

    if not dry_run:
        root_file.Write("", ROOT.TObject.kOverwrite)

    root_file.Close()

    return changed, skipped


def expand_inputs(patterns: list[str]) -> list[str]:
    files = []

    for pattern in patterns:
        matched = glob.glob(pattern)
        if matched:
            files.extend(matched)
        elif os.path.exists(pattern):
            files.append(pattern)
        else:
            print(f"[warning] No files matched: {pattern}")

    files = sorted(set(files))
    return files


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Clone histograms in ROOT files whose internal _v#### tag does not "
            "match the _v#### tag in the ROOT filename."
        )
    )

    parser.add_argument(
        "inputs",
        nargs="+",
        help="ROOT files or glob patterns, e.g. 'res2droot_*/res2dPlots*.root'",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be changed, but do not modify files.",
    )

    parser.add_argument(
        "--overwrite-existing",
        action="store_true",
        help="Overwrite the corrected histogram name if it already exists.",
    )

    args = parser.parse_args()

    files = expand_inputs(args.inputs)

    if not files:
        print("No input ROOT files found.")
        return 1

    total_changed = 0
    total_skipped = 0

    for path in files:
        changed, skipped = process_file(
            path,
            dry_run=args.dry_run,
            overwrite_existing=args.overwrite_existing,
        )
        total_changed += changed
        total_skipped += skipped

    print("\n=== Summary")
    print(f"Files processed: {len(files)}")
    print(f"Histograms cloned/renamed: {total_changed}")
    print(f"Skipped: {total_skipped}")

    if args.dry_run:
        print("Dry run only. No files were modified.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
