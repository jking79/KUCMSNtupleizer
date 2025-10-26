#!/usr/bin/env python3
import ROOT
import sys
import os

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
EOS_BASE = "/store/user/lpcsusylep/jaking/KUCMSNtuple"

# -----------------------------------------------------------------------------
def check_root_file(fullpath):

    """Try to open and read a ROOT file via xrootd."""
    f = ROOT.TFile.Open(fullpath, "READ")
    if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
        return False, "Unreadable or corrupted"
    try:
        keys = f.GetListOfKeys()
        _ = [key.GetName() for key in keys]
    except Exception as e:
        f.Close()
        return False, f"Read error: {e}"
    f.Close()
    return True, "OK"

# -----------------------------------------------------------------------------
def main():
    if len(sys.argv) != 3:
        print("Usage: python3 check_root_files_lpc.py <subpath> <filelist.txt>")
        print("Example:")
        print("  python3 check_root_files_lpc.py "
              "kucmsntuple_MET_R18_SVIPM100_v31/MET/kucmsntuple_MET_R18_SVIPM100_v31_MET_AOD_Run2018A "
              "filelist.txt")
        sys.exit(1)

    subpath = sys.argv[1].strip("/")
    listfile = sys.argv[2]

    # Load file names from the provided text list
    if not os.path.exists(listfile):
        print(f"Error: file list {listfile} not found")
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
        # Construct the full EOS path
        full_eos_path = f"root://cmseos.fnal.gov/{EOS_BASE}/{subpath}/{filename}"
        ok, msg = check_root_file(full_eos_path)
        status = "✅ OK" if ok else "❌ BAD"
        print(f"[{i:03d}] {status}  {filename}  -->  {msg}")
        if not ok:
            bad_files.append((i, filename, msg))

    # Summary
    print("\n------------------------------------------------------")
    print(f"Checked {len(files)} files total.")
    if bad_files:
        print(f"⚠️  {len(bad_files)} problematic files found:")
        for idx, name, reason in bad_files:
            print(f"   ID {idx:03d}: {name} ({reason})")
    else:
        print("✅ All files readable and valid.")
    print("------------------------------------------------------\n")

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    ROOT.gErrorIgnoreLevel = ROOT.kWarning  # suppress ROOT noise
    main()

