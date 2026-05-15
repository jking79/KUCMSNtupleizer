#!/usr/bin/env python3

from ROOT import TFile
import argparse
import math


def open_tree(filename, treename="tree/llpgtree"):

    f = TFile.Open(filename, "READ")

    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open file: {filename}")

    t = f.Get(treename)

    if not t:
        raise RuntimeError(f"Could not find tree '{treename}' in file: {filename}")

    return f, t


def is_float_like(x):
    return isinstance(x, float)

def values_equal(a, b, rel_tol=1e-6, abs_tol=1e-9):

    # ----------------------------------------
    # Handle NaN explicitly
    # ----------------------------------------
    try:
        if math.isnan(a) and math.isnan(b):
            return True
    except TypeError:
        pass

    # ----------------------------------------
    # Handle exact equality
    # ----------------------------------------
    if a == b:
        return True

    # ----------------------------------------
    # Float comparison
    # ----------------------------------------
    if isinstance(a, float) or isinstance(b, float):
        try:
            return math.isclose(float(a), float(b), rel_tol=rel_tol, abs_tol=abs_tol)
        except Exception:
            return False

    return False

def object_to_list(obj):

    # Handles std::vector-like objects
    if hasattr(obj, "size") and hasattr(obj, "at"):
        return [obj.at(i) for i in range(obj.size())]

    # Handles scalar ROOT proxy values
    return obj

def branch_values_equal(v1, v2, rel_tol=1e-6, abs_tol=1e-9):

    v1 = object_to_list(v1)
    v2 = object_to_list(v2)

    if isinstance(v1, list) or isinstance(v2, list):

        if not isinstance(v1, list) or not isinstance(v2, list):
            return False

        if len(v1) != len(v2):
            return False

        for a, b in zip(v1, v2):
            if not values_equal(a, b, rel_tol, abs_tol):
                return False

        return True

    return values_equal(v1, v2, rel_tol, abs_tol)


def format_value(v, max_items=8):

    v = object_to_list(v)

    if isinstance(v, list):
        if len(v) > max_items:
            return f"{v[:max_items]} ... size={len(v)}"
        return str(v)

    return str(v)


def looks_like_gen_branch(branch_name):

    gen_tokens = [
        "gen",
        "Gen",
        "GEN",
        "susy",
        "SUSY",
        "xsec",
        "XSec",
        "lhe",
        "LHE",
    ]

    return any(tok in branch_name for tok in gen_tokens)


def compare_trees(
    file1,
    file2,
    treename="tree/llpgtree",
    max_diffs_per_branch=5,
    max_events=-1,
    rel_tol=1e-6,
    abs_tol=1e-9,
):

    f1, t1 = open_tree(file1, treename)
    f2, t2 = open_tree(file2, treename)

    branches1 = sorted([b.GetName() for b in t1.GetListOfBranches()])
    branches2 = sorted([b.GetName() for b in t2.GetListOfBranches()])

    set1 = set(branches1)
    set2 = set(branches2)

    only1 = sorted(set1 - set2)
    only2 = sorted(set2 - set1)
    common = sorted(set1 & set2)

    print("============================================================")
    print("Tree comparison")
    print("============================================================")
    print(f"File 1 : {file1}")
    print(f"File 2 : {file2}")
    print(f"Tree   : {treename}")
    print()

    print(f"[INFO] Entries file 1 : {t1.GetEntries()}")
    print(f"[INFO] Entries file 2 : {t2.GetEntries()}")
    print(f"[INFO] Branches file 1: {len(branches1)}")
    print(f"[INFO] Branches file 2: {len(branches2)}")
    print(f"[INFO] Common branches : {len(common)}")
    print()

    if only1:
        print("[WARN] Branches only in file 1:")
        for b in only1:
            print(f"  {b}")
        print()

    if only2:
        print("[WARN] Branches only in file 2:")
        for b in only2:
            print(f"  {b}")
        print()

    nentries = min(t1.GetEntries(), t2.GetEntries())

    if max_events > 0:
        nentries = min(nentries, max_events)

    if t1.GetEntries() != t2.GetEntries():
        print("[WARN] Trees have different numbers of entries")
        print()

    differing_branches = {}
    checked_events = 0

    # Disable all branches first, then enable one branch at a time
    t1.SetBranchStatus("*", 0)
    t2.SetBranchStatus("*", 0)

    for bname in common:

        t1.SetBranchStatus(bname, 1)
        t2.SetBranchStatus(bname, 1)

        diffs = []

        for ientry in range(nentries):

            t1.GetEntry(ientry)
            t2.GetEntry(ientry)

            v1 = getattr(t1, bname)
            v2 = getattr(t2, bname)

            if not branch_values_equal(v1, v2, rel_tol, abs_tol):

                if len(diffs) < max_diffs_per_branch:
                    diffs.append(
                        {
                            "entry": ientry,
                            "file1": format_value(v1),
                            "file2": format_value(v2),
                        }
                    )

        if diffs:
            differing_branches[bname] = diffs

        t1.SetBranchStatus(bname, 0)
        t2.SetBranchStatus(bname, 0)

    print("============================================================")
    print("Differing branches")
    print("============================================================")

    if not differing_branches:
        print("[PASS] No differing common branches found.")
    else:
        for bname, diffs in differing_branches.items():

            gen_flag = "GEN-LIKE" if looks_like_gen_branch(bname) else "NON-GEN"

            print(f"[DIFF] {bname}  ({gen_flag})")

            for d in diffs:
                print(f"  entry {d['entry']}")
                print(f"    file1: {d['file1']}")
                print(f"    file2: {d['file2']}")

            print()

    non_gen_diffs = [
        b for b in differing_branches
        if not looks_like_gen_branch(b)
    ]

    print("============================================================")
    print("Summary")
    print("============================================================")
    print(f"Compared entries          : {nentries}")
    print(f"Differing common branches : {len(differing_branches)}")
    print(f"Non-Gen-like differences  : {len(non_gen_diffs)}")

    if non_gen_diffs:
        print()
        print("[WARN] Non-Gen-like branches differ:")
        for b in non_gen_diffs:
            print(f"  {b}")
    else:
        print()
        print("[PASS] All differing common branches look Gen-related.")

    f1.Close()
    f2.Close()


def parse_args():

    parser = argparse.ArgumentParser(
        description="Compare all branches in tree/llpgtree between two ROOT files"
    )

    parser.add_argument("file1", help="First ROOT file")
    parser.add_argument("file2", help="Second ROOT file")

    parser.add_argument(
        "-t",
        "--tree",
        default="tree/llpgtree",
        help="Tree path inside ROOT file"
    )

    parser.add_argument(
        "-n",
        "--maxEvents",
        type=int,
        default=-1,
        help="Maximum number of events to compare"
    )

    parser.add_argument(
        "--maxDiffs",
        type=int,
        default=5,
        help="Maximum number of example differences to print per branch"
    )

    parser.add_argument(
        "--relTol",
        type=float,
        default=1e-6,
        help="Relative tolerance for float comparisons"
    )

    parser.add_argument(
        "--absTol",
        type=float,
        default=1e-9,
        help="Absolute tolerance for float comparisons"
    )

    return parser.parse_args()


def main():

    args = parse_args()

    compare_trees(
        file1=args.file1,
        file2=args.file2,
        treename=args.tree,
        max_diffs_per_branch=args.maxDiffs,
        max_events=args.maxEvents,
        rel_tol=args.relTol,
        abs_tol=args.absTol,
    )


if __name__ == "__main__":
    main()
