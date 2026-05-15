#!/usr/bin/env python3

import ROOT
import argparse
import subprocess
import fnmatch

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

EOS_MGM = "root://cmseos.fnal.gov"
XRD_PREFIX = "root://cmseos.fnal.gov/"

EXPECTED_TTREE = "tree/llpgtree"
REQUIRED_BRANCHES = ["Vertex_px", "Evt_run", "Met_sumEt", "Photon_nPhoton"]

# -----------------------------------------------------------------------------


def eos_cmd(cmd):
    try:
        return subprocess.check_output(
            cmd,
            shell=True,
            text=True,
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as e:
        print("[ERROR] EOS command failed")
        print(f"[ERROR] Command: {cmd}")
        print(f"[ERROR] Exit code: {e.returncode}")
        print("[ERROR] Output:")
        print(e.output)
        raise


def list_matching_subdirs(eos_parent_dir, match_pattern):
    cmd = f"eos {EOS_MGM} ls {eos_parent_dir}"
    print(f"[INFO] Running: {cmd}")

    out = eos_cmd(cmd)

    subdirs = []

    for line in out.splitlines():
        name = line.strip()

        if not name:
            continue

        if fnmatch.fnmatch(name, match_pattern):
            subdirs.append(f"{eos_parent_dir.rstrip('/')}/{name}")

    return sorted(subdirs)


def find_root_files_recursive(eos_dir, match_pattern=""):

    if match_pattern:
        search_dirs = list_matching_subdirs(eos_dir, match_pattern)
    else:
        search_dirs = [eos_dir]

    print(f"[INFO] Matched search directories: {len(search_dirs)}")
    for d in search_dirs:
        print(f"  {d}")

    files = []

    for search_dir in search_dirs:
        cmd = f"eos {EOS_MGM} find -f {search_dir}"
        print(f"[INFO] Running: {cmd}")

        out = eos_cmd(cmd)

        for line in out.splitlines():
            line = line.strip()

            if not line.endswith(".root"):
                continue

            if line.startswith("/store/"):
                files.append(XRD_PREFIX + line)
            else:
                files.append(line)

    return sorted(files)


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


def main():
    parser = argparse.ArgumentParser(
        description="Recursively validate ntuple ROOT files under an EOS directory."
    )

    parser.add_argument(
        "eos_dir",
        help="EOS parent directory to scan, e.g. /store/user/lpcsusylep/jaking/KUCMSNtuple"
    )

    parser.add_argument(
        "--match",
        default="",
        help=(
            "Optional wildcard pattern for immediate subdirectories under eos_dir, "
            "e.g. 'kucmsntuple_GJets_R18_SVHPM100_MiniAOD*'"
        )
    )

    parser.add_argument(
        "--max-files",
        type=int,
        default=-1,
        help="Maximum ROOT files to check. Default: all."
    )

    parser.add_argument(
        "--only",
        default="",
        help="Only process ROOT files whose full path contains this string."
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

    root_files = find_root_files_recursive(args.eos_dir, args.match)

    if args.only:
        root_files = [f for f in root_files if args.only in f]

    if args.max_files > 0:
        root_files = root_files[:args.max_files]

    print("------------------------------------------------------")
    print(f"EOS directory: {args.eos_dir}")
    print(f"Subdir match:   {args.match if args.match else '[none]'}")
    print(f"ROOT files found: {len(root_files)}")
    print(f"Expected TTree: {EXPECTED_TTREE}")
    print(f"Required branches: {REQUIRED_BRANCHES}")
    print("------------------------------------------------------")

    total_files = 0
    total_bad = 0
    bad_records = []

    for i, full_eos_path in enumerate(root_files, start=1):
        total_files += 1

        res = check_root_file(full_eos_path)

        if res["ok"]:
            if not args.quiet_ok:
                print(f"[{i:05d}/{len(root_files):05d}] OK   {full_eos_path}")
        else:
            total_bad += 1

            print(f"[{i:05d}/{len(root_files):05d}] BAD  {full_eos_path}")
            for issue in res["issues"]:
                print(f"       -> {issue}")

            bad_records.append({
                "path": full_eos_path,
                "issues": res["issues"],
            })

    with open(args.bad_out, "w") as f:
        for rec in bad_records:
            f.write(f"PATH: {rec['path']}\n")
            for issue in rec["issues"]:
                f.write(f"ISSUE: {issue}\n")
            f.write("\n")

    print()
    print("------------------------------------------------------")
    print("Final summary")
    print("------------------------------------------------------")
    print(f"ROOT files checked: {total_files}")
    print(f"Problematic files: {total_bad}")
    print(f"Bad-file report: {args.bad_out}")
    print("------------------------------------------------------")


if __name__ == "__main__":
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    main()
