#!/usr/bin/env python3
"""Read-only KUCMS EOS diagnostic scan; explicitly approved deletion is separate."""
import argparse
import csv
import math
import re
import subprocess
from collections import Counter
from pathlib import Path

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

# ===================== USER-EDITABLE DEFAULTS =====================
HOST = "cmseos.fnal.gov"
DIRECTORY = (
    "/store/user/lpcsusylep/jaking/KUCMSNtuple/"
    "kucmsntuple_MET_R18D_SVHPM100_MiniAODv2_v34/MET/"
    "kucmsntuple_MET_R18D_SVHPM100_MiniAODv2_v34_"
    "MET_Run2018D-UL2018_MiniAODv2_GT36-v1/260831_230957/0001"
)
TREE_NAME = "tree/llpgtree"
REQUIRED_BRANCHES = ("Met_pt", "Jet_pt", "Photon_nPhoton", "ECALRecHit_energy")
OPTIONAL_BRANCHES = ("Met_CPt", "Evt_run", "Evt_luminosityBlock", "Evt_event")
LARGE_MB = 500.0
LARGE_KB_PER_EVENT = 100.0  # tune after seeing normal production-file values
MET_GEV = 1000.0
N_JETS = 50
N_PHOTONS = 20
N_RECHITS = 10000
SAMPLE_EVENTS = 25
REPORT = "eos_root_scan.csv"
EVENT_REPORT = "eos_suspicious_events.csv"
APPROVALS = "eos_delete_approvals.csv"
# ==================================================================

FIELDS = [
    "path", "filename", "size_bytes", "size_mb", "entries", "kb_per_event",
    "status", "reasons", "branches_missing", "sampled_events", "full_scan",
    "events_checked", "problem_events", "max_met", "max_corrected_met",
    "max_jets", "max_photons", "max_rechits", "error"
]
EVENT_FIELDS = ["path", "entry", "run", "lumi", "event",
                "met", "corrected_met", "jets", "photons", "rechits", "reasons"]
APPROVAL_FIELDS = ["path", "size_bytes", "delete"]


def xrd(*args):
    p = subprocess.run(["xrdfs", HOST, *args], text=True, capture_output=True)
    if p.returncode:
        raise RuntimeError((p.stderr or p.stdout).strip())
    return p.stdout


def remote_size(path):
    output = xrd("stat", path)
    m = re.search(r"(?m)^\s*Size:\s*(\d+)\b", output)
    if not m:
        raise RuntimeError("Cannot parse xrdfs stat Size: " + output[:300])
    return int(m.group(1))


def paths():
    return sorted(p.strip() for p in xrd("ls", DIRECTORY).splitlines()
                  if p.strip().endswith(".root"))


def indices(n, count):
    if n == 0:
        return []
    if n <= count:
        return range(n)
    return sorted(set(round(k * (n - 1) / (count - 1))
                      for k in range(count)))


def valid_number(v):
    v = float(v)
    return v if math.isfinite(v) else None


def inspect(path, event_writer, sample_count):
    row = dict.fromkeys(FIELDS, "")
    row.update(path=path, filename=Path(path).name, status="OK",
               sampled_events=0, full_scan="NO", events_checked=0,
               problem_events=0)
    reasons = set()
    file = None
    try:
        size = remote_size(path)
        row["size_bytes"] = size
        row["size_mb"] = round(size / 1048576, 3)
        if row["size_mb"] > LARGE_MB:
            reasons.add("LARGE_FILE")
        file = ROOT.TFile.Open(f"root://{HOST}/{path}", "READ")
        if not file or file.IsZombie():
            row["status"] = "INVALID"
            reasons.add("OPEN_FAILED")
            return row
        if file.TestBit(ROOT.TFile.kRecovered):
            reasons.add("RECOVERED_FILE")
            row["status"] = "STRUCTURE_WARNING"
        tree = file.Get(TREE_NAME)
        if not tree or not tree.InheritsFrom("TTree"):
            row["status"] = "INVALID"
            reasons.add("MISSING_EVENT_TREE")
            return row
        branches = {b.GetName() for b in tree.GetListOfBranches()}
        missing = sorted(set(REQUIRED_BRANCHES) - branches)
        row["branches_missing"] = ";".join(missing)
        if missing:
            row["status"] = "INVALID"
            reasons.add("MISSING_REQUIRED_BRANCHES")
            return row
        n = int(tree.GetEntries())
        row["entries"] = n
        row["kb_per_event"] = round(size / (1024 * n), 3) if n else ""
        if n == 0:
            reasons.add("EMPTY_TREE")
        if n and float(row["kb_per_event"]) > LARGE_KB_PER_EVENT:
            reasons.add("LARGE_BYTES_PER_EVENT")
        if n == 0:
            return row

        selected = [b for b in REQUIRED_BRANCHES + OPTIONAL_BRANCHES if b in branches]
        tree.SetBranchStatus("*", 0)
        for b in selected:
            tree.SetBranchStatus(b, 1)

        maxima = {"max_met": None, "max_corrected_met": None,
                  "max_jets": None, "max_photons": None, "max_rechits": None}

        def event_at(i, write=False):
            if tree.GetEntry(i) <= 0:
                raise RuntimeError(f"Failed reading selected branches at entry {i}")
            met = valid_number(getattr(tree, "Met_pt"))
            corrected = (valid_number(getattr(tree, "Met_CPt"))
                         if "Met_CPt" in branches else None)
            jets = len(getattr(tree, "Jet_pt"))
            photons = int(getattr(tree, "Photon_nPhoton"))
            rechits = len(getattr(tree, "ECALRecHit_energy"))
            vals = (met, corrected, jets, photons, rechits)
            for key, val in zip(maxima, vals):
                if val is not None:
                    maxima[key] = max(val, maxima[key]) if maxima[key] is not None else val
            problems = []
            if met is None or (corrected is None and "Met_CPt" in branches):
                problems.append("NONFINITE_MET")
            if met is not None and met > MET_GEV:
                problems.append("HIGH_MET")
            if jets > N_JETS:
                problems.append("MANY_JETS")
            if photons > N_PHOTONS:
                problems.append("MANY_PHOTONS")
            if rechits > N_RECHITS:
                problems.append("MANY_RECHITS")
            if write and problems:
                event_writer.writerow({
                    "path": path, "entry": i,
                    "run": int(getattr(tree, "Evt_run")) if "Evt_run" in branches else "",
                    "lumi": int(getattr(tree, "Evt_luminosityBlock"))
                    if "Evt_luminosityBlock" in branches else "",
                    "event": int(getattr(tree, "Evt_event"))
                    if "Evt_event" in branches else "",
                    "met": met, "corrected_met": corrected,
                    "jets": jets, "photons": photons, "rechits": rechits,
                    "reasons": ";".join(problems)
                })
            return problems

        sampled = list(indices(n, sample_count))
        sample_problems = 0
        for i in sampled:
            sample_problems += bool(event_at(i))
        row["sampled_events"] = len(sampled)
        candidate = bool(reasons) or sample_problems > 0
        if candidate:
            row["full_scan"] = "YES"
            row["problem_events"] = 0
            for i in range(n):
                row["problem_events"] += bool(event_at(i, write=True))
            row["events_checked"] = n
        else:
            row["events_checked"] = len(sampled)
        row.update({k: round(v, 3) if isinstance(v, float) else v
                    for k, v in maxima.items() if v is not None})
        if row["problem_events"] or sample_problems:
            reasons.add("SUSPICIOUS_EVENT_CONTENT")
        if row["status"] == "OK" and reasons:
            row["status"] = "REVIEW"
    except Exception as exc:
        row["error"] = str(exc)
        row["status"] = "ERROR"
        reasons.add("SCAN_ERROR")
    finally:
        row["reasons"] = ";".join(sorted(reasons))
        if file:
            file.Close()
    return row


def scan(args):
    file_paths = paths()
    print(f"Found {len(file_paths)} ROOT files in {DIRECTORY}", flush=True)
    counts = Counter()
    with open(REPORT, "w", newline="") as fout, \
         open(EVENT_REPORT, "w", newline="") as fevents, \
         open(APPROVALS, "w", newline="") as fapprove:
        report = csv.DictWriter(fout, FIELDS)
        events = csv.DictWriter(fevents, EVENT_FIELDS)
        approvals = csv.DictWriter(fapprove, APPROVAL_FIELDS)
        report.writeheader()
        events.writeheader()
        approvals.writeheader()
        for i, path in enumerate(file_paths, 1):
            row = inspect(path, events, args.sample_events)
            report.writerow(row)
            approvals.writerow({"path": path, "size_bytes": row["size_bytes"],
                                "delete": "NO"})
            fout.flush()
            fevents.flush()
            fapprove.flush()
            counts[row["status"]] += 1
            print(f"[{i}/{len(file_paths)}] {row['status']:8} "
                  f"{str(row['size_mb']):>10} MB  "
                  f"{str(row['entries']):>8} entries  {row['filename']}",
                  flush=True)
    print(f"Results: {REPORT}, {EVENT_REPORT}, {APPROVALS}")
    print("Status counts:", dict(counts))
    print("No files deleted.")


def delete(args):
    with open(APPROVALS, newline="") as f:
        rows = list(csv.DictReader(f))
    prefix = DIRECTORY.rstrip("/") + "/"
    approved = [r for r in rows if r["delete"].strip().upper() == "YES"]
    if any(r["delete"].strip().upper() not in ("YES", "NO") for r in rows):
        raise RuntimeError("Approvals CSV must contain only YES or NO in delete column")
    if not approved:
        print("No files approved.")
        return
    if len({r["path"] for r in approved}) != len(approved):
        raise RuntimeError("Duplicate approved paths")
    # Check all approved files before allowing any deletion.
    for row in approved:
        path = row["path"]
        if (not path.startswith(prefix) or "/" in path[len(prefix):]
                or not path.endswith(".root") or not row["size_bytes"].isdigit()):
            raise RuntimeError(f"Unsafe or unverified approval: {path}")
        current = remote_size(path)
        if current != int(row["size_bytes"]):
            raise RuntimeError(f"File changed since scan: {path} "
                               f"({row['size_bytes']} -> {current} bytes)")
    for row in approved:
        print(f"  {row['size_bytes']:>14} bytes  {row['path']}")
    if not args.execute:
        print(f"\nDRY RUN: {len(approved)} files approved. Nothing deleted.")
        return
    if input(f"\nType DELETE {len(approved)} FILES to confirm: ") != \
            f"DELETE {len(approved)} FILES":
        print("Cancelled.")
        return
    audit = "eos_deletion_audit.csv"
    with open(audit, "a", newline="") as f:
        writer = csv.writer(f)
        if f.tell() == 0:
            writer.writerow(["path", "scanned_size", "result", "error"])
        for row in approved:
            path = row["path"]
            try:
                if remote_size(path) != int(row["size_bytes"]):
                    raise RuntimeError("Size changed immediately before deletion")
                xrd("rm", path)
                writer.writerow([path, row["size_bytes"], "DELETED", ""])
                print("DELETED", path)
            except Exception as exc:
                writer.writerow([path, row["size_bytes"], "FAILED", str(exc)])
                print("FAILED", path, exc)
            f.flush()
    print("Audit:", audit)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    s = sub.add_parser("scan")
    s.add_argument("--sample-events", type=int, default=SAMPLE_EVENTS)
    d = sub.add_parser("delete")
    d.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.command == "scan":
        if args.sample_events < 1:
            parser.error("--sample-events must be >= 1")
        scan(args)
    else:
        delete(args)


if __name__ == "__main__":
    main()

