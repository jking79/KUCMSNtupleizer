#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path

EOS_MGM = "root://cmseos.fnal.gov"
EOS_BASE = "/store/user/lpcsusylep/jaking/KUCMSNtuple"

def run_cmd(cmd, dry_run: bool):
    if dry_run:
        print("[DRY-RUN]", " ".join(cmd))
        return 0
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        print("[ERROR]", " ".join(cmd))
        if r.stderr.strip():
            print("   ", r.stderr.strip())
    return r.returncode

def main():
    ap = argparse.ArgumentParser(
        description="Remove EOS files listed in a text file from LPC EOS (cmseos.fnal.gov)."
    )
    ap.add_argument("listfile", help="Text file with one relative path per line.")
    ap.add_argument("start_dir", help='Start directory under EOS_BASE, e.g. "kucmsntuple_MET_R17_SVIPM100_v31/MET/".')
    ap.add_argument("--execute", action="store_true", help="Actually delete files (default is dry-run).")
    ap.add_argument("--recursive", action="store_true",
                    help="Use recursive remove (-r). Only needed if your list includes directories.")
    ap.add_argument("--stop-on-error", action="store_true", help="Stop at first failure.")
    args = ap.parse_args()

    list_path = Path(args.listfile)
    if not list_path.exists():
        raise SystemExit(f"ERROR: listfile not found: {list_path}")

    # Normalize start_dir and build base LFN prefix
    start_dir = args.start_dir.strip().strip("/")
    base_prefix = f"{EOS_BASE}/{start_dir}".rstrip("/")

    # Read list entries
    entries = []
    with list_path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            entries.append(s)

    if not entries:
        print("No entries found in list file.")
        return

    dry_run = not args.execute
    print("EOS MGM:", EOS_MGM)
    print("EOS_BASE:", EOS_BASE)
    print("start_dir:", start_dir)
    print("entries:", len(entries))
    print("mode:", "EXECUTE (deleting!)" if not dry_run else "DRY-RUN (no deletes)")
    print("recursive:", args.recursive)
    print("")

    # Build and run eosrm commands (expects LFN /store/... paths)
    n_ok = 0
    n_fail = 0

    for i, rel in enumerate(entries, 1):
        rel = rel.lstrip("/")  # ensure it doesn't wipe base_prefix
        lfn = f"{base_prefix}/{rel}"

        # eosrm is an LPC alias/wrapper; call through a login shell so aliases resolve.
        # If you prefer not to rely on aliases, see note below.
        cmd = ["bash", "-lc", f"eosrm {'-r ' if args.recursive else ''}{lfn}"]

        rc = run_cmd(cmd, dry_run=dry_run)
        if rc == 0:
            n_ok += 1
        else:
            n_fail += 1
            if args.stop_on_error:
                break

        if i % 200 == 0:
            print(f"...processed {i}/{len(entries)}")

    print("\nSummary:")
    print("  attempted:", n_ok + n_fail)
    print("  ok:", n_ok)
    print("  failed:", n_fail)

if __name__ == "__main__":
    main()

