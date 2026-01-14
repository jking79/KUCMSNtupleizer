#!/usr/bin/env python3
import ROOT
import sys
import os
from datetime import datetime

def get_runs_from_file(filename):
    """Extract unique run numbers from a ROOT file using PyROOT."""
    runs = set()
    try:
        f = ROOT.TFile.Open(filename)
        if not f or f.IsZombie():
            print(f"[ERROR] Cannot open file: {filename}")
            return runs

        tree = f.Get("Events")
        if not tree:
            print(f"[WARNING] No 'Events' tree found in: {filename}")
            f.Close()
            return runs

        # Prepare variable for reading branch
        run_num = ROOT.std.vector('unsigned int')()
        tree.SetBranchStatus("*", 0)
        tree.SetBranchStatus("run", 1)
        tree.SetBranchAddress("run", run_num)

        # Loop over entries
        for entry in range(tree.GetEntries()):
            tree.GetEntry(entry)
            runs.add(run_num[0])  # store run number
            if entry % 100000 == 0 and entry != 0:
                print(f"  Processed {entry} entries...")

        f.Close()
    except Exception as e:
        print(f"[ERROR] Failed reading {filename}: {e}")
    return runs


def main():
    if len(sys.argv) < 2:
        print("Usage: python check_runs_in_files_root.py <file_list.txt>")
        sys.exit(1)

    filelist_path = sys.argv[1]

    # Read list of files
    with open(filelist_path) as f:
        files = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    if not files:
        print("[ERROR] No valid file paths found.")
        sys.exit(1)

    all_runs = set()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_name = f"run_check_report_{timestamp}.txt"

    with open(report_name, "w") as out:
        out.write(f"Run check report generated: {timestamp}\n")
        out.write(f"Input list: {filelist_path}\n\n")

        for filename in files:
            print(f"Processing: {filename}")
            runs = get_runs_from_file(filename)
            if runs:
                out.write(f"{filename}:\n  Runs: {sorted(runs)}\n\n")
            else:
                out.write(f"{filename}:\n  [No runs or unreadable file]\n\n")
            all_runs.update(runs)

        out.write(f"\nAll unique runs across all files:\n{sorted(all_runs)}\n")

    print(f"\n✅ Report saved to {report_name}")


if __name__ == "__main__":
    main()

