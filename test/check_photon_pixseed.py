#!/usr/bin/env python3

from ROOT import TFile
import argparse
import subprocess
import os


def get_eos_files(eos_dir, xrootd_prefix="root://cmseos.fnal.gov/"):

    cmd = f'eos root://cmseos.fnal.gov find -f {eos_dir}'
    print(f"[INFO] Running: {cmd}")

    out = subprocess.check_output(
        cmd,
        shell=True,
        text=True,
        stderr=subprocess.STDOUT,
    )

    files = []

    for line in out.splitlines():
        line = line.strip()

        if not line.endswith(".root"):
            continue

        if line.startswith("/store/"):
            files.append(xrootd_prefix + line)
        else:
            files.append(line)

    return sorted(files)


def check_file(
    filename,
    tree_name="kuSkimTree",
    branch_name="photon_PixSeed",
):

    f = TFile.Open(filename, "READ")

    if not f or f.IsZombie():
        return {
            "status": "bad_file",
            "file": filename,
            "n_true": 0,
            "n_false": 0,
            "n_photons": 0,
        }

    t = f.Get(tree_name)

    if not t:
        f.Close()
        return {
            "status": "missing_tree",
            "file": filename,
            "n_true": 0,
            "n_false": 0,
            "n_photons": 0,
        }

    if not t.GetBranch(branch_name):
        f.Close()
        return {
            "status": "missing_branch",
            "file": filename,
            "n_true": 0,
            "n_false": 0,
            "n_photons": 0,
        }

    t.SetBranchStatus("*", 0)
    t.SetBranchStatus(branch_name, 1)

    n_true = 0
    n_false = 0

    nentries = t.GetEntries()

    for ientry in range(nentries):

        t.GetEntry(ientry)

        pixseed_vec = getattr(t, branch_name)

        # vector<bool>
        for i in range(pixseed_vec.size()):
            if bool(pixseed_vec.at(i)):
                n_true += 1
            else:
                n_false += 1

    f.Close()

    n_photons = n_true + n_false

    if n_photons == 0:
        status = "no_photons"
    elif n_true > 0 and n_false > 0:
        status = "ok"
    elif n_true > 0 and n_false == 0:
        status = "all_true"
    elif n_true == 0 and n_false > 0:
        status = "all_false"
    else:
        status = "unknown"

    return {
        "status": status,
        "file": filename,
        "n_true": n_true,
        "n_false": n_false,
        "n_photons": n_photons,
    }


def main():

    parser = argparse.ArgumentParser(
        description="Check photon_PixSeed content in skim ROOT files"
    )

    parser.add_argument(
        "-d",
        "--eosDir",
        default="/store/user/lpcsusylep/malazaro/KUCMSSkims/skims_v50",
        help="EOS directory containing skim ROOT files",
    )

    parser.add_argument(
        "-t",
        "--tree",
        default="kuSkimTree",
        help="Tree name",
    )

    parser.add_argument(
        "-b",
        "--branch",
        default="photon_PixSeed",
        help="Branch name",
    )

    parser.add_argument(
        "-n",
        "--maxFiles",
        type=int,
        default=-1,
        help="Maximum number of files to check, useful for testing",
    )

    args = parser.parse_args()

    files = get_eos_files(args.eosDir)

    if args.maxFiles > 0:
        files = files[:args.maxFiles]

    print(f"[INFO] Found {len(files)} ROOT files")
    print()

    bad_results = []
    no_photon_files = 0
    ok_files = 0

    for idx, filename in enumerate(files, start=1):

        print(f"[CHECK] {idx}/{len(files)} {filename}")

        result = check_file(
            filename=filename,
            tree_name=args.tree,
            branch_name=args.branch,
        )

        status = result["status"]

        if status == "ok":
            ok_files += 1
        elif status == "no_photons":
            no_photon_files += 1
        else:
            bad_results.append(result)

            print(
                f"[BAD] {status}: "
                f"nTrue={result['n_true']} "
                f"nFalse={result['n_false']} "
                f"nPhotons={result['n_photons']}"
            )

    print()
    print("============================================================")
    print("Summary")
    print("============================================================")
    print(f"Files checked       : {len(files)}")
    print(f"OK files            : {ok_files}")
    print(f"No-photon files     : {no_photon_files}")
    print(f"Problem files       : {len(bad_results)}")
    print()

    if bad_results:
        print("Problem file list:")
        for r in bad_results:
            print(
                f"{r['status']:15s} "
                f"nTrue={r['n_true']:10d} "
                f"nFalse={r['n_false']:10d} "
                f"nPhotons={r['n_photons']:10d} "
                f"{r['file']}"
            )
    else:
        print("[PASS] No all-true/all-false photon_PixSeed files found.")


if __name__ == "__main__":
    main()
