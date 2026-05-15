#!/usr/bin/env python3

import argparse
import subprocess
from collections import defaultdict
import ROOT

EOS_HOST = "cmseos.fnal.gov"
EOS_BASE = "/store/user/lpcsusylep/jaking/KUCMSNtuple"
TREE_NAME = "tree/llpgtree"
BRANCH_NAME = "Photon_pixelSeed"


def make_empty_summary():
    return {
        "files_found": 0,
        "ok": 0,
        "bad": 0,
        "unreadable": 0,
        "missing_tree": 0,
        "missing_branch": 0,
        "no_photons": 0,
        "all_false": 0,
        "entries": 0,
        "photons": 0,
        "pixel_seed_true": 0,
    }


def add_result_to_summary(summary, res):
    summary["files_found"] += 1
    summary["entries"] += res["nentries"]
    summary["photons"] += res["nphotons"]
    summary["pixel_seed_true"] += res["n_pixel_seed_true"]

    if res["ok"]:
        summary["ok"] += 1
    else:
        summary["bad"] += 1

        if res["issue"] == "file unreadable or zombie":
            summary["unreadable"] += 1
        elif res["issue"].startswith("missing tree"):
            summary["missing_tree"] += 1
        elif res["issue"].startswith("missing branch"):
            summary["missing_branch"] += 1
        elif res["issue"] == "no photons found":
            summary["no_photons"] += 1
        elif res["issue"] == "all Photon_pixelSeed values are false":
            summary["all_false"] += 1


def eos_find_root_files(eos_subdir):
    eos_dir = f"{EOS_BASE.rstrip('/')}/{eos_subdir.strip('/')}"
    cmd = ["xrdfs", EOS_HOST, "ls", "-R", eos_dir]
    out = subprocess.check_output(cmd, text=True)

    root_files = [
        f"root://{EOS_HOST}/{line.strip()}"
        for line in out.splitlines()
        if line.strip().endswith(".root")
    ]

    return eos_dir, root_files


def get_first_subdir_under_input(root_path, eos_dir):
    """
    Example:
      eos_dir:
        /store/.../KUCMSNtuple/kucmsntuple_GJets_R18_SVHPM100_MiniAOD_v34

      root_path:
        root://cmseos.fnal.gov//store/.../KUCMSNtuple/kucmsntuple_GJets_R18_SVHPM100_MiniAOD_v34/GJets_HT-100To200.../.../kucmsntuple_512.root

      returns:
        GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8
    """

    prefix = f"root://{EOS_HOST}/"
    path = root_path

    if path.startswith(prefix):
        path = path[len(prefix):]

    # normalize leading double slash cases
    path = "/" + path.lstrip("/")
    eos_dir = "/" + eos_dir.strip("/")

    rel = path[len(eos_dir):].strip("/") if path.startswith(eos_dir) else ""

    if not rel:
        return "__NO_SUBDIR__"

    return rel.split("/", 1)[0]


def check_photon_pixel_seed(root_path, max_entries=-1):
    result = {
        "path": root_path,
        "ok": False,
        "nentries": 0,
        "nphotons": 0,
        "n_pixel_seed_true": 0,
        "issue": "",
    }

    f = ROOT.TFile.Open(root_path, "READ")
    if not f or f.IsZombie():
        result["issue"] = "file unreadable or zombie"
        return result

    tree = f.Get(TREE_NAME)
    if not tree:
        result["issue"] = f"missing tree {TREE_NAME}"
        f.Close()
        return result

    if not tree.GetBranch(BRANCH_NAME):
        result["issue"] = f"missing branch {BRANCH_NAME}"
        f.Close()
        return result

    nentries = tree.GetEntries()
    result["nentries"] = nentries

    n_to_check = nentries if max_entries < 0 else min(nentries, max_entries)

    for i in range(n_to_check):
        tree.GetEntry(i)
        pixel_seed_vec = getattr(tree, BRANCH_NAME)

        for val in pixel_seed_vec:
            result["nphotons"] += 1
            if bool(val):
                result["n_pixel_seed_true"] += 1

    f.Close()

    if result["nphotons"] == 0:
        result["issue"] = "no photons found"
        return result

    if result["n_pixel_seed_true"] == 0:
        result["issue"] = "all Photon_pixelSeed values are false"
        return result

    result["ok"] = True
    return result


def print_summary_block(title, summary):
    print(title)
    print("-" * len(title))
    print(f"  files found:                    {summary['files_found']}")
    print(f"  files OK:                       {summary['ok']}")
    print(f"  files BAD:                      {summary['bad']}")
    print()
    print("  BAD breakdown:")
    print(f"    unreadable/zombie:            {summary['unreadable']}")
    print(f"    missing tree:                 {summary['missing_tree']}")
    print(f"    missing Photon_pixelSeed:     {summary['missing_branch']}")
    print(f"    no photons:                   {summary['no_photons']}")
    print(f"    all Photon_pixelSeed false:   {summary['all_false']}")
    print()
    print("  counts scanned:")
    print(f"    entries:                      {summary['entries']}")
    print(f"    photons:                      {summary['photons']}")
    print(f"    Photon_pixelSeed true:        {summary['pixel_seed_true']}")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Scan EOS ntuple ROOT files and check Photon_pixelSeed is not always false."
    )

    parser.add_argument(
        "eos_subdir",
        help="Subdirectory under /store/user/lpcsusylep/jaking/KUCMSNtuple/"
    )
    parser.add_argument(
        "--max-entries",
        type=int,
        default=-1,
        help="Maximum entries to scan per file. Default: all."
    )
    parser.add_argument(
        "--bad-out",
        default="bad_pixelSeed_files.txt",
        help="Output file listing bad ROOT files."
    )
    parser.add_argument(
        "--section-out",
        default="pixelSeed_section_summary.txt",
        help="Output file containing compact per-section summary."
    )

    args = parser.parse_args()

    eos_dir, root_files = eos_find_root_files(args.eos_subdir)

    total = make_empty_summary()
    sections = defaultdict(make_empty_summary)
    bad = []

    for root_path in root_files:
        section = get_first_subdir_under_input(root_path, eos_dir)
        res = check_photon_pixel_seed(root_path, args.max_entries)

        add_result_to_summary(total, res)
        add_result_to_summary(sections[section], res)

        if not res["ok"]:
            bad.append({
                "section": section,
                **res
            })

    with open(args.bad_out, "w") as f:
        for res in bad:
            f.write(f"SECTION: {res['section']}\n")
            f.write(f"PATH: {res['path']}\n")
            f.write(f"ISSUE: {res['issue']}\n")
            f.write(f"ENTRIES: {res['nentries']}\n")
            f.write(f"PHOTONS: {res['nphotons']}\n")
            f.write(f"PIXELSEED_TRUE: {res['n_pixel_seed_true']}\n")
            f.write("\n")

    with open(args.section_out, "w") as f:
        f.write("section files_found ok bad unreadable missing_tree missing_branch no_photons all_false entries photons pixelSeedTrue\n")
        for section in sorted(sections):
            s = sections[section]
            f.write(
                f"{section} "
                f"{s['files_found']} {s['ok']} {s['bad']} "
                f"{s['unreadable']} {s['missing_tree']} {s['missing_branch']} "
                f"{s['no_photons']} {s['all_false']} "
                f"{s['entries']} {s['photons']} {s['pixel_seed_true']}\n"
            )

    print()
    print("======================================================")
    print("Photon_pixelSeed EOS scan summary")
    print("======================================================")
    print(f"EOS base:       {EOS_BASE}")
    print(f"EOS subdir:     {args.eos_subdir}")
    print(f"Full EOS dir:   {eos_dir}")
    print(f"Sections found: {len(sections)}")
    print()

    print_summary_block("TOTAL", total)

    print("Section summaries")
    print("=================")
    for section in sorted(sections):
        s = sections[section]
        print(
            f"{section}: "
            f"files={s['files_found']} "
            f"OK={s['ok']} "
            f"BAD={s['bad']} "
            f"allFalse={s['all_false']} "
            f"missingBranch={s['missing_branch']} "
            f"photons={s['photons']} "
            f"pixelSeedTrue={s['pixel_seed_true']}"
        )

    print()
    print(f"Bad-file report:       {args.bad_out}")
    print(f"Section summary table: {args.section_out}")
    print("======================================================")
    print()


if __name__ == "__main__":
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    main()
