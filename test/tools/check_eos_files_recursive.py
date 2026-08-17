#!/usr/bin/env python3
import ROOT
import sys
import subprocess

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
EOS_HOST = "cmseos.fnal.gov"
EOS_BASE = "/store/user/lpcsusylep/jaking"
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
def find_root_files(eos_dir):
    """Return all .root files in eos_dir and all directories below it."""

    cmd = ["xrdfs", EOS_HOST, "ls", "-R", eos_dir]

    try:
        proc = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
    except FileNotFoundError:
        print("Error: 'xrdfs' was not found in PATH.")
        print("Run this script in an environment with the XRootD client tools available.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Error scanning EOS directory: {eos_dir}")
        if e.stderr:
            print(e.stderr.strip())
        sys.exit(1)

    return sorted(
        line.strip()
        for line in proc.stdout.splitlines()
        if line.strip().lower().endswith(".root")
    )

# -----------------------------------------------------------------------------
def main():
    if len(sys.argv) != 2:
        print("Usage: python3 check_eos_files_recursive.py <directory>")
        print("Examples:")
        print("  python3 check_eos_files_recursive.py KUCMSNtuple/gogoGZ_2022_Fast_PUlib1")
        print("  python3 check_eos_files_recursive.py /store/user/lpcsusylep/jaking/KUCMSNtuple/gogoGZ_2022_Fast_PUlib1")
        sys.exit(1)

    requested_dir = sys.argv[1].rstrip("/")

    # Accept either a full EOS /store/... path or a path relative to EOS_BASE.
    if requested_dir.startswith("/store/"):
        eos_dir = requested_dir
    else:
        eos_dir = f"{EOS_BASE}/{requested_dir.lstrip('/')}"

    print(f"\n🔎 Scanning recursively for ROOT files under:")
    print(f"   {eos_dir}\n")

    files = find_root_files(eos_dir)

    if not files:
        print("No .root files found in that directory or below it.")
        sys.exit(1)

    print(f"Found {len(files)} ROOT files. Beginning checks...\n")

    bad_files = []
    for i, eos_path in enumerate(files, start=1):
        full_eos_path = f"root://{EOS_HOST}/{eos_path}"
        res = check_root_file(full_eos_path)
        status = "✅ OK" if res["ok"] else "❌ BAD"

        print(f"[{i:04d}/{len(files):04d}] {status}  {eos_path}")

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
            print(f"   {idx:04d}: {b['path']}")
            for issue in b["issues"]:
                print(f"         ↳ {issue}")
    else:
        print("✅ All files readable and valid.")
    print("------------------------------------------------------\n")

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    ROOT.gErrorIgnoreLevel = ROOT.kWarning  # suppress noisy ROOT messages
    main()
