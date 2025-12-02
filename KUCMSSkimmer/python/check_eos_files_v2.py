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

    # --- Optional: checksum check ---
    #if not f.TestBit(ROOT.TFile.kWriteChecksum):
    #    result["issues"].append("Checksum not written (cannot verify integrity)")

    f.Close()

    result["ok"] = (len(result["issues"]) == 0)
    return result

# -----------------------------------------------------------------------------
def main():
    if len(sys.argv) != 3:
        print("Usage: python3 check_eos_files_v2 <subpath> <filelist.txt>")
        print("Example:")
        print("  python3 check_eos_files_v2 "
              "kucmsntuple_MET_R18_SVIPM100_v31/ "
              "ntuple_master_lists/MET_R18_SVIPM100_v31_MET_AOD_Run2018A.txt")
        sys.exit(1)

    subpath = sys.argv[1].strip("/")
    listfile = sys.argv[2]

    # Load list of ROOT files
    if not os.path.exists(listfile):
        print(f"Error: file list '{listfile}' not found")
        sys.exit(1)

    with open(listfile) as f:
        files = [line.strip() for line in f if line.strip()]
    if not files:
        print("No files found in list.")
        sys.exit(1)

    print(f"\n🔍 Checking {len(files)} ROOT files under EOS:")
    print(f"   Base: {EOS_BASE}")
    print(f"   Subdir: {subpath}\n")

    bad_files = []
    for i, filename in enumerate(files, start=1):
        full_eos_path = f"root://cmseos.fnal.gov/{EOS_BASE}/{subpath}/{filename}"
        res = check_root_file(full_eos_path)
        status = "✅ OK" if res["ok"] else "❌ BAD"
        print(f"[{i:03d}] {status}  {filename}")
        if not res["ok"]:
            for issue in res["issues"]:
                print(f"        ⚠️  {issue}")
            bad_files.append(res)

    # --- Summary ---
    print("\n------------------------------------------------------")
    print(f"Checked {len(files)} files total.")
    if bad_files:
        print(f"⚠️  {len(bad_files)} problematic files found:")
        for idx, b in enumerate(bad_files, 1):
            print(f"   {idx:03d}: {b['path']}")
            for issue in b["issues"]:
                print(f"         ↳ {issue}")
    else:
        print("✅ All files readable and valid.")
    print("------------------------------------------------------\n")

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    ROOT.gErrorIgnoreLevel = ROOT.kWarning  # suppress noisy ROOT messages
    main()

