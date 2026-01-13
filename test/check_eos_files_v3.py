#!/usr/bin/env python3
import ROOT
import sys
import os

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
EOS_BASE = "/store/user/lpcsusylep/jaking/KUCMSNtuple"
EXPECTED_TTREE = "tree/llpgtree"
REQUIRED_BRANCHES = ["Vertex_px", "Evt_run", "Met_sumEt", "Photon_nPhoton"]  # customize as needed

# -----------------------------------------------------------------------------
def check_root_file(fullpath):
    """Perform detailed validity checks on a ROOT file via xrootd."""

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
        result["issues"].append("File was recovered (incomplete write)")

    # --- Tree presence ---
    t = f.Get(EXPECTED_TTREE)
    if not t or not isinstance(t, ROOT.TTree):
        result["issues"].append(f"Missing or invalid TTree '{EXPECTED_TTREE}'")
        f.Close()
        return result

    # --- Entry count ---
    nentries = t.GetEntries()
    if nentries <= 0:
        result["issues"].append("TTree has zero entries")

    # --- Required branches ---
    branches = [b.GetName() for b in t.GetListOfBranches()]
    missing = [b for b in REQUIRED_BRANCHES if b not in branches]
    if missing:
        result["issues"].append(f"Missing branches: {missing}")

    # --- Try reading first entry ---
    try:
        t.GetEntry(0)
    except Exception as e:
        result["issues"].append(f"Error reading entry 0: {e}")

    f.Close()

    result["ok"] = (len(result["issues"]) == 0)
    return result


# -----------------------------------------------------------------------------
def main():
    if len(sys.argv) not in (3, 4):
        print("Usage: python3 check_eos_files_v2 <subpath> <filelist.txt> [start_line]")
        print("Example:")
        print("  python3 check_eos_files_v2 "
              "kucmsntuple_MET_R18_SVIPM100_v31/ "
              "ntuple_master_lists/MET_R18_SVIPM100_v31_MET_AOD_Run2018A.txt 501")
        sys.exit(1)

    subpath = sys.argv[1].strip("/")
    listfile = sys.argv[2]

    # Optional: starting line number (1-based)
    start_line = 1
    if len(sys.argv) == 4:
        try:
            start_line = int(sys.argv[3])
        except ValueError:
            print(f"ERROR: start_line must be an integer, got '{sys.argv[3]}'")
            sys.exit(1)

    if start_line < 1:
        print("ERROR: start_line must be >= 1")
        sys.exit(1)

    # Load list of ROOT files
    if not os.path.exists(listfile):
        print(f"Error: file list '{listfile}' not found")
        sys.exit(1)

    with open(listfile) as f:
        all_files = [line.strip() for line in f if line.strip()]

    if not all_files:
        print("No files found in list.")
        sys.exit(1)

    # Slice starting at requested line
    # start_line is 1-based; Python list is 0-based
    start_idx = start_line - 1
    if start_idx >= len(all_files):
        print(f"ERROR: start_line={start_line} exceeds file length ({len(all_files)} lines)")
        sys.exit(1)

    files = all_files[start_idx:]
    skipped = start_idx

    print(f"\n>>> Checking ROOT files under EOS:")
    print(f"    Base:   {EOS_BASE}")
    print(f"    Subdir: {subpath}")
    print(f"    List:   {listfile}")
    print(f"    Start:  line {start_line} (skipped {skipped} entries)")
    print(f"    Will check: {len(files)} files\n")

    bad_files = []

    # Option 1: Keep original line numbering in the list
    for j, filename in enumerate(files, start=start_line):
        full_eos_path = f"root://cmseos.fnal.gov/{EOS_BASE}/{subpath}/{filename}"
        res = check_root_file(full_eos_path)
        status = "OK" if res["ok"] else "BAD"
        print(f"[{j:05d}] {status:3s}  {filename}")
        if not res["ok"]:
            for issue in res["issues"]:
                print(f"         - {issue}")
            bad_files.append(res)

    # --- Summary ---
    print("\n------------------------------------------------------")
    print(f"Checked {len(files)} files (from line {start_line} onward).")
    if bad_files:
        print(f"WARNING: {len(bad_files)} problematic files found:")
        for idx, b in enumerate(bad_files, 1):
            print(f"   {idx:03d}: {b['path']}")
            for issue in b["issues"]:
                print(f"         - {issue}")
    else:
        print("All files readable and valid.")
    print("------------------------------------------------------\n")


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    ROOT.gErrorIgnoreLevel = ROOT.kWarning  # suppress noisy ROOT messages
    main()

