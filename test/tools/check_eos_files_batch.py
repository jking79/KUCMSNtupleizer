#!/usr/bin/env python3

import ROOT
import sys
import os
import argparse
from pathlib import Path

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

EOS_BASE = "root://cmseos.fnal.gov//store/user/lpcsusylep/jaking"
EXPECTED_TTREE = "tree/llpgtree"
REQUIRED_BRANCHES = ["Vertex_px", "Evt_run", "Met_sumEt", "Photon_nPhoton"]

# -----------------------------------------------------------------------------

def check_root_file(fullpath):
    result = {
        "path": fullpath,
        "ok": False,
        "issues": []
    }

    f = ROOT.TFile.Open(fullpath, "READ")

    if not f or f.IsZombie():
        result["issues"].append("File unreadable or zombie")
        return result

    if f.TestBit(ROOT.TFile.kRecovered):
        result["issues"].append("File was recovered/incomplete write")

    t = f.Get(EXPECTED_TTREE)
    if not t or not isinstance(t, ROOT.TTree):
        result["issues"].append(f"Missing or invalid TTree '{EXPECTED_TTREE}'")
        f.Close()
        return result

    if t.GetEntries() <= 0:
        result["issues"].append("TTree has zero entries")

    branches = [b.GetName() for b in t.GetListOfBranches()]
    missing = [b for b in REQUIRED_BRANCHES if b not in branches]
    if missing:
        result["issues"].append(f"Missing branches: {missing}")

    try:
        t.GetEntry(0)
    except Exception as e:
        result["issues"].append(f"Error reading entry 0: {e}")

    f.Close()
    result["ok"] = len(result["issues"]) == 0
    return result

# -----------------------------------------------------------------------------

def parse_master_list(master_file):
    samples = []

    with open(master_file) as f:
        for lineno, line in enumerate(f, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) < 2:
                print(f"[WARN] malformed master line {lineno}: {line}")
                continue

            eos_subdir = parts[0].strip("/")
            sample_list = parts[1]

            samples.append({
                "lineno": lineno,
                "eos_subdir": eos_subdir,
                "sample_list": sample_list,
                "raw": line,
            })

    return samples

# -----------------------------------------------------------------------------

def read_sample_file_list(sample_list_path):
    if not os.path.exists(sample_list_path):
        return None, [f"Sample list file not found: {sample_list_path}"]

    files = []

    with open(sample_list_path) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            files.append(line)

    return files, []

# -----------------------------------------------------------------------------

def build_eos_path(eos_subdir, root_relpath):
    return f"{EOS_BASE}/{eos_subdir.strip('/')}/{root_relpath.lstrip('/')}"

# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Validate ROOT files referenced by a KUCMS master list."
    )

    parser.add_argument("master_list", help="Master list file, e.g. KUCMS_Ntuple_Master_SMS_Sig_Files_List_v34p1.txt")
    parser.add_argument(
        "--list-dir",
        default=".",
        help="Directory containing per-sample list files. Default: current directory."
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=-1,
        help="Maximum ROOT files to check per sample list. Default: all."
    )
    parser.add_argument(
        "--only",
        default="",
        help="Only process master-list rows whose sample-list filename contains this string."
    )
    parser.add_argument(
        "--quiet-ok",
        action="store_true",
        help="Do not print every OK file."
    )
    parser.add_argument(
        "--bad-out",
        default="bad_root_files.txt",
        help="Output text file listing bad ROOT files."
    )

    args = parser.parse_args()

    samples = parse_master_list(args.master_list)

    if args.only:
        samples = [s for s in samples if args.only in s["sample_list"]]

    print("------------------------------------------------------")
    print(f"Master list: {args.master_list}")
    print(f"Samples found: {len(samples)}")
    print(f"EOS base: {EOS_BASE}")
    print(f"Sample list dir: {args.list_dir}")
    print("------------------------------------------------------")

    total_files = 0
    total_bad = 0
    missing_lists = []

    bad_records = []

    for isample, sample in enumerate(samples, start=1):
        sample_list_path = Path(args.list_dir) / sample["sample_list"]

        print()
        print("======================================================")
        print(f"[Sample {isample}/{len(samples)}]")
        print(f"Master line: {sample['lineno']}")
        print(f"EOS subdir:  {sample['eos_subdir']}")
        print(f"List file:   {sample_list_path}")
        print("======================================================")

        root_files, list_issues = read_sample_file_list(sample_list_path)

        if list_issues:
            for issue in list_issues:
                print(f"  BAD: {issue}")
            missing_lists.append(sample_list_path)
            total_bad += 1
            bad_records.append({
                "sample": str(sample_list_path),
                "path": "MISSING_LIST",
                "issues": list_issues,
            })
            continue

        if args.max_files > 0:
            root_files = root_files[:args.max_files]

        print(f"Files in list: {len(root_files)}")

        sample_bad = 0

        for i, relpath in enumerate(root_files, start=1):
            total_files += 1

            full_eos_path = build_eos_path(sample["eos_subdir"], relpath)
            res = check_root_file(full_eos_path)

            if res["ok"]:
                if not args.quiet_ok:
                    print(f"[{i:04d}] OK   {relpath}")
            else:
                sample_bad += 1
                total_bad += 1
                print(f"[{i:04d}] BAD  {relpath}")
                for issue in res["issues"]:
                    print(f"       -> {issue}")

                bad_records.append({
                    "sample": str(sample_list_path),
                    "path": full_eos_path,
                    "issues": res["issues"],
                })

        print(f"Sample summary: checked={len(root_files)} bad={sample_bad}")

    with open(args.bad_out, "w") as f:
        for rec in bad_records:
            f.write(f"SAMPLE_LIST: {rec['sample']}\n")
            f.write(f"PATH: {rec['path']}\n")
            for issue in rec["issues"]:
                f.write(f"ISSUE: {issue}\n")
            f.write("\n")

    print()
    print("------------------------------------------------------")
    print("Final summary")
    print("------------------------------------------------------")
    print(f"Samples processed: {len(samples)}")
    print(f"ROOT files checked: {total_files}")
    print(f"Problematic files/lists: {total_bad}")
    print(f"Missing sample lists: {len(missing_lists)}")
    print(f"Bad-file report: {args.bad_out}")
    print("------------------------------------------------------")

# -----------------------------------------------------------------------------

if __name__ == "__main__":
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    main()
