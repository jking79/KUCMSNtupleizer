#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
import time
from collections import defaultdict
import ROOT

EOS_HOST = "cmseos.fnal.gov"
EOS_BASE = "/store/user/lpcsusylep/jaking/KUCMSNtuple"
TREE_NAME = "tree/llpgtree"
BRANCH_NAME = "Photon_pixelSeed"
DEFAULT_VERSION_SUFFIX = "_v34"
CHECKPOINT_SCHEMA_VERSION = 2


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


def summary_from_maybe_dict(obj):
    s = make_empty_summary()
    if isinstance(obj, dict):
        for key in s:
            try:
                s[key] = int(obj.get(key, 0))
            except (TypeError, ValueError):
                s[key] = 0
    return s


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


def merge_summary(dst, src):
    for key in make_empty_summary():
        dst[key] += int(src.get(key, 0))


def run_xrdfs_ls(eos_dir, recursive=False):
    cmd = ["xrdfs", EOS_HOST, "ls"]
    if recursive:
        cmd.append("-R")
    cmd.append(eos_dir)
    return subprocess.check_output(cmd, text=True)


def eos_find_matching_version_dirs(version_suffix=DEFAULT_VERSION_SUFFIX):
    """Return immediate children of EOS_BASE whose directory names end with version_suffix."""
    out = run_xrdfs_ls(EOS_BASE, recursive=False)

    eos_dirs = []
    for line in out.splitlines():
        path = line.strip().rstrip("/")
        if not path:
            continue
        if os.path.basename(path).endswith(version_suffix):
            eos_dirs.append(path)

    return sorted(eos_dirs)


def eos_find_root_files_in_dir(eos_dir):
    """Return all ROOT files recursively under one fully qualified EOS directory."""
    out = run_xrdfs_ls(eos_dir, recursive=True)

    root_files = [
        f"root://{EOS_HOST}/{line.strip()}"
        for line in out.splitlines()
        if line.strip().endswith(".root")
    ]

    return root_files


def eos_dir_from_subdir(eos_subdir):
    return f"{EOS_BASE.rstrip('/')}/{eos_subdir.strip('/')}"


def get_first_subdir_under_input(root_path, eos_dir):
    prefix = f"root://{EOS_HOST}/"
    path = root_path

    if path.startswith(prefix):
        path = path[len(prefix):]

    path = "/" + path.lstrip("/")
    eos_dir = "/" + eos_dir.strip("/")

    rel = path[len(eos_dir):].strip("/") if path.startswith(eos_dir) else ""

    if not rel:
        return "__NO_SUBDIR__"

    return rel.split("/", 1)[0]


def get_scan_section(root_path, eos_dir):
    """Label results as top-level v34 directory plus the first child section/dataset."""
    top_dir = os.path.basename(eos_dir.rstrip("/"))
    child = get_first_subdir_under_input(root_path, eos_dir)

    if child == "__NO_SUBDIR__":
        return top_dir

    return f"{top_dir}/{child}"


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


def format_elapsed(seconds):
    seconds = int(max(0, seconds))
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)

    if hours:
        return f"{hours:d}h {minutes:02d}m {secs:02d}s"
    if minutes:
        return f"{minutes:d}m {secs:02d}s"
    return f"{secs:d}s"


def print_status(message):
    print(message, flush=True)


def should_print_progress(done, total_files, last_progress_time, args):
    if total_files <= 0:
        return False
    if done >= total_files:
        return True

    by_files = (
        args.progress_every_files > 0
        and done > 0
        and done % args.progress_every_files == 0
    )
    by_time = (
        args.progress_every_seconds > 0
        and time.monotonic() - last_progress_time >= args.progress_every_seconds
    )

    return by_files or by_time


def should_checkpoint(done, total_files, last_checkpoint_time, args):
    if total_files <= 0:
        return False
    if done >= total_files:
        return True

    by_files = (
        args.checkpoint_every_files > 0
        and done > 0
        and done % args.checkpoint_every_files == 0
    )
    by_time = (
        args.checkpoint_every_seconds > 0
        and time.monotonic() - last_checkpoint_time >= args.checkpoint_every_seconds
    )

    return by_files or by_time


def safe_file_component(name):
    return name.replace("/", "__").replace(" ", "_")


def report_paths(report_dir, top_dir_name):
    base = safe_file_component(top_dir_name)
    return {
        "live": os.path.join(report_dir, f"{base}.live.txt"),
        "done": os.path.join(report_dir, f"{base}.done.txt"),
    }


def read_json_checkpoint(path):
    try:
        with open(path, "r") as f:
            first = f.read(1)
            f.seek(0)
            if first != "{":
                return None
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return None


def report_is_completed(done_path):
    """Return True if done report explicitly says the directory completed."""
    data = read_json_checkpoint(done_path)
    if data is not None:
        return data.get("status") in ("DONE", "DONE_NO_ROOT_FILES")

    # Backward-compatible check for older text-style .done.txt files.
    if not os.path.exists(done_path):
        return False
    try:
        with open(done_path, "r") as f:
            for line in f:
                line = line.strip()
                if line in ("STATUS: DONE", "STATUS: DONE_NO_ROOT_FILES"):
                    return True
                if line.startswith("STATUS:"):
                    return False
    except OSError:
        return False
    return False


def load_live_checkpoint(live_path, top_dir_name, eos_dir):
    """Load a JSON .live.txt checkpoint. Returns None if unavailable or incompatible."""
    data = read_json_checkpoint(live_path)
    if data is None:
        return None

    if data.get("schema_version") != CHECKPOINT_SCHEMA_VERSION:
        return None
    if data.get("top_dir_name") != top_dir_name:
        return None
    if data.get("eos_dir") != eos_dir:
        return None
    if data.get("tree") != TREE_NAME or data.get("branch") != BRANCH_NAME:
        return None
    if not isinstance(data.get("root_files"), list):
        return None
    if not isinstance(data.get("processed_paths"), list):
        return None

    return data


TEXT_SUMMARY_KEY_MAP = {
    "files_found": "files_found",
    "ok": "ok",
    "bad": "bad",
    "unreadable": "unreadable",
    "missing_tree": "missing_tree",
    "missing_branch": "missing_branch",
    "no_photons": "no_photons",
    "all_false": "all_false",
    "entries": "entries",
    "photons": "photons",
    "pixelSeedTrue": "pixel_seed_true",
    "pixel_seed_true": "pixel_seed_true",
}


def parse_key_value_line(line):
    """Parse text report lines like 'FILES_SCANNED: 11000'."""
    if ":" not in line:
        return None, None
    key, value = line.split(":", 1)
    return key.strip(), value.strip()


def int_or_zero(value):
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return 0


def load_text_live_report(live_path, top_dir_name, eos_dir):
    """
    Load an older text-style .live.txt report.

    This format does not contain processed_paths, so recovery is by position:
    FILES_SCANNED gives the number of leading ROOT files to skip after the fresh
    xrdfs discovery list is rebuilt.
    """
    if not os.path.exists(live_path):
        return None

    try:
        with open(live_path, "r") as f:
            lines = [line.rstrip("\n") for line in f]
    except OSError:
        return None

    header = {}
    summary = make_empty_summary()
    bad_results = []
    current_bad = None
    in_summary = False
    in_bad_files = False

    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue

        if line == "SUMMARY":
            in_summary = True
            in_bad_files = False
            continue
        if line == "BAD FILES":
            in_summary = False
            in_bad_files = True
            continue
        if set(line) == {"-"}:
            continue

        key, value = parse_key_value_line(line)
        if key is None:
            continue

        if in_bad_files:
            if key == "SECTION":
                if current_bad:
                    bad_results.append(current_bad)
                current_bad = {"section": value}
            elif current_bad is not None:
                if key == "PATH":
                    current_bad["path"] = value
                elif key == "ISSUE":
                    current_bad["issue"] = value
                elif key == "ENTRIES":
                    current_bad["nentries"] = int_or_zero(value)
                elif key == "PHOTONS":
                    current_bad["nphotons"] = int_or_zero(value)
                elif key == "PIXELSEED_TRUE":
                    current_bad["n_pixel_seed_true"] = int_or_zero(value)
            continue

        if in_summary:
            mapped = TEXT_SUMMARY_KEY_MAP.get(key)
            if mapped is not None:
                summary[mapped] = int_or_zero(value)
            continue

        header[key] = value

    if current_bad:
        bad_results.append(current_bad)

    status = header.get("STATUS", "")
    if status not in ("RUNNING", "RUNNING_DISCOVERY_DONE"):
        return None
    if header.get("DIRECTORY") != top_dir_name:
        return None
    if header.get("EOS_PATH") != eos_dir:
        return None
    if header.get("TREE") != TREE_NAME or header.get("BRANCH") != BRANCH_NAME:
        return None

    files_scanned = int_or_zero(header.get("FILES_SCANNED", 0))
    files_total_discovered = int_or_zero(header.get("FILES_TOTAL_DISCOVERED", 0))
    if files_scanned <= 0:
        return None

    # If the SUMMARY block is absent/incomplete, at least preserve the count of
    # files skipped so the per-top-directory progress is not reset to zero.
    if summary["files_found"] == 0:
        summary["files_found"] = files_scanned

    normalized_bad_results = []
    for res in bad_results:
        normalized = {
            "section": res.get("section", top_dir_name),
            "path": res.get("path", ""),
            "ok": False,
            "nentries": int_or_zero(res.get("nentries", 0)),
            "nphotons": int_or_zero(res.get("nphotons", 0)),
            "n_pixel_seed_true": int_or_zero(res.get("n_pixel_seed_true", 0)),
            "issue": res.get("issue", "recovered from text live report"),
        }
        normalized_bad_results.append(normalized)

    return {
        "files_scanned": files_scanned,
        "files_total_discovered": files_total_discovered,
        "summary": summary,
        "bad_results": normalized_bad_results,
    }


def make_checkpoint_data(
    status,
    top_dir_name,
    eos_dir,
    root_files,
    processed_paths,
    summary,
    sections,
    bad_results,
    elapsed_seconds,
):
    files_scanned = len(processed_paths)
    files_total = len(root_files)
    pct = 0.0 if files_total <= 0 else 100.0 * files_scanned / files_total

    return {
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "status": status,
        "top_dir_name": top_dir_name,
        "eos_dir": eos_dir,
        "tree": TREE_NAME,
        "branch": BRANCH_NAME,
        "files_scanned": files_scanned,
        "files_total_discovered": files_total,
        "percent_scanned": pct,
        "elapsed_seconds": int(max(0, elapsed_seconds)),
        "elapsed": format_elapsed(elapsed_seconds),
        "root_files": root_files,
        "processed_paths": processed_paths,
        "summary": dict(summary),
        "sections": {k: dict(v) for k, v in sections.items()},
        "bad_results": bad_results,
    }


def write_checkpoint(path, data):
    """Write JSON checkpoint atomically."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp_path = f"{path}.tmp"
    with open(tmp_path, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write("\n")
    os.replace(tmp_path, path)


def write_topdir_checkpoint(
    path,
    status,
    top_dir_name,
    eos_dir,
    root_files,
    processed_paths,
    summary,
    sections,
    bad_results,
    elapsed_seconds,
):
    data = make_checkpoint_data(
        status,
        top_dir_name,
        eos_dir,
        root_files,
        processed_paths,
        summary,
        sections,
        bad_results,
        elapsed_seconds,
    )
    write_checkpoint(path, data)


def write_global_reports(bad_out, section_out, topdir_out, bad, sections, top_dirs):
    with open(bad_out, "w") as f:
        for res in bad:
            f.write(f"SECTION: {res['section']}\n")
            f.write(f"PATH: {res['path']}\n")
            f.write(f"ISSUE: {res['issue']}\n")
            f.write(f"ENTRIES: {res['nentries']}\n")
            f.write(f"PHOTONS: {res['nphotons']}\n")
            f.write(f"PIXELSEED_TRUE: {res['n_pixel_seed_true']}\n")
            f.write("\n")

    with open(section_out, "w") as f:
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

    with open(topdir_out, "w") as f:
        f.write("top_dir files_found ok bad unreadable missing_tree missing_branch no_photons all_false entries photons pixelSeedTrue\n")
        for top_dir in sorted(top_dirs):
            s = top_dirs[top_dir]
            f.write(
                f"{top_dir} "
                f"{s['files_found']} {s['ok']} {s['bad']} "
                f"{s['unreadable']} {s['missing_tree']} {s['missing_branch']} "
                f"{s['no_photons']} {s['all_false']} "
                f"{s['entries']} {s['photons']} {s['pixel_seed_true']}\n"
            )


def recover_state_from_live_checkpoint(data):
    root_files = list(data["root_files"])
    processed_paths = list(data["processed_paths"])
    summary = summary_from_maybe_dict(data.get("summary", {}))

    sections = defaultdict(make_empty_summary)
    for section, s in data.get("sections", {}).items():
        sections[section] = summary_from_maybe_dict(s)

    bad_results = list(data.get("bad_results", []))
    return root_files, processed_paths, summary, sections, bad_results


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Scan EOS ntuple ROOT files and check Photon_pixelSeed is not always false. "
            "By default, scans only immediate EOS_BASE child directories ending in _v34."
        )
    )

    parser.add_argument(
        "eos_subdir",
        nargs="?",
        default=None,
        help=(
            "Optional single subdirectory under /store/user/lpcsusylep/jaking/KUCMSNtuple/. "
            "If omitted, all immediate children ending in --version-suffix are scanned."
        ),
    )
    parser.add_argument(
        "--version-suffix",
        default=DEFAULT_VERSION_SUFFIX,
        help="Only scan top-level EOS directories whose names end with this suffix. Default: _v34.",
    )
    parser.add_argument(
        "--max-entries",
        type=int,
        default=-1,
        help="Maximum entries to scan per file. Default: all.",
    )
    parser.add_argument(
        "--bad-out",
        default="bad_pixelSeed_files.txt",
        help="Output file listing bad ROOT files accumulated during this run, including resumed checkpoint state.",
    )
    parser.add_argument(
        "--section-out",
        default="pixelSeed_section_summary.txt",
        help="Output file containing compact per-section summary accumulated during this run, including resumed checkpoint state.",
    )
    parser.add_argument(
        "--topdir-out",
        default="pixelSeed_topdir_summary.txt",
        help="Output file containing compact per-top-directory summary accumulated during this run, including resumed checkpoint state.",
    )
    parser.add_argument(
        "--report-dir",
        default="pixelSeed_scan_reports",
        help="Directory for per-top-directory JSON .live.txt and .done.txt reports. Default: pixelSeed_scan_reports.",
    )
    parser.add_argument(
        "--skip-completed",
        action="store_true",
        help=(
            "Skip .done.txt directories and resume compatible JSON .live.txt checkpoints. "
            "If an older text-style .live.txt report is found, use FILES_SCANNED "
            "to skip that many leading ROOT files and continue the scan."
        ),
    )
    parser.add_argument(
        "--checkpoint-every-files",
        type=int,
        default=500,
        help=(
            "Rewrite the current directory JSON .live.txt every N newly scanned ROOT files. "
            "Use 0 or a negative value to disable file-count checkpointing. Default: 25."
        ),
    )
    parser.add_argument(
        "--checkpoint-every-seconds",
        type=int,
        default=900,
        help=(
            "Rewrite the current directory JSON .live.txt at least every N seconds. "
            "Use 0 or a negative value to disable time-based checkpointing. Default: 300."
        ),
    )
    parser.add_argument(
        "--progress-every-files",
        type=int,
        default=500,
        help=(
            "Print an in-directory progress update every N total completed ROOT files. "
            "Use 0 or a negative value to disable file-count progress. Default: 25."
        ),
    )
    parser.add_argument(
        "--progress-every-seconds",
        type=int,
        default=900,
        help=(
            "Print an in-directory progress update at least every N seconds, even if "
            "the file-count interval has not been reached. Use 0 or a negative value "
            "to disable time-based progress. Default: 60."
        ),
    )

    args = parser.parse_args()
    os.makedirs(args.report_dir, exist_ok=True)

    if args.eos_subdir:
        eos_dir = eos_dir_from_subdir(args.eos_subdir)
        if not os.path.basename(eos_dir.rstrip("/")).endswith(args.version_suffix):
            raise SystemExit(
                f"Refusing to scan non-matching directory: {eos_dir}\n"
                f"Directory basename must end with {args.version_suffix!r}."
            )
        scan_dirs = [eos_dir]
    else:
        scan_dirs = eos_find_matching_version_dirs(args.version_suffix)

    if not scan_dirs:
        raise SystemExit(
            f"No directories ending in {args.version_suffix!r} found under {EOS_BASE}"
        )

    print_status("======================================================")
    print_status("Starting Photon_pixelSeed EOS scan")
    print_status("======================================================")
    print_status(f"EOS base:          {EOS_BASE}")
    print_status(f"Version suffix:    {args.version_suffix}")
    print_status(f"Report directory:  {args.report_dir}")
    print_status(f"Skip completed:    {args.skip_completed}")
    print_status(f"Directories found: {len(scan_dirs)}")
    for idx, d in enumerate(scan_dirs, start=1):
        print_status(f"  [{idx}/{len(scan_dirs)}] {d}")
    print_status("")

    total = make_empty_summary()
    sections = defaultdict(make_empty_summary)
    top_dirs = defaultdict(make_empty_summary)
    bad = []
    scan_start_time = time.monotonic()
    skipped = []
    resumed = []

    for idir, eos_dir in enumerate(scan_dirs, start=1):
        top_dir_name = os.path.basename(eos_dir.rstrip("/"))
        paths = report_paths(args.report_dir, top_dir_name)

        if args.skip_completed and report_is_completed(paths["done"]):
            skipped.append(top_dir_name)
            print_status(f"[{idir}/{len(scan_dirs)}] Skipping completed directory: {top_dir_name}")
            continue

        dir_start_time = time.monotonic()
        dir_sections = defaultdict(make_empty_summary)
        dir_bad = []
        processed_paths = []
        root_files = None

        live_data = None
        text_live_data = None
        if args.skip_completed:
            live_data = load_live_checkpoint(paths["live"], top_dir_name, eos_dir)
            if live_data is None:
                text_live_data = load_text_live_report(paths["live"], top_dir_name, eos_dir)

        if live_data is not None:
            root_files, processed_paths, dir_summary, dir_sections, dir_bad = recover_state_from_live_checkpoint(live_data)
            n_root_files = len(root_files)
            processed_set = set(processed_paths)
            resumed.append(top_dir_name)

            top_dirs[top_dir_name] = dir_summary
            merge_summary(total, dir_summary)
            for section, s in dir_sections.items():
                merge_summary(sections[section], s)
            bad.extend(dir_bad)

            print_status("------------------------------------------------------")
            print_status(f"[{idir}/{len(scan_dirs)}] Resuming directory from JSON checkpoint: {top_dir_name}")
            print_status(f"EOS path: {eos_dir}")
            print_status(f"Checkpoint: {paths['live']}")
            print_status(
                f"Already processed {len(processed_set)}/{n_root_files} files "
                f"({100.0 * len(processed_set) / n_root_files if n_root_files else 0.0:.1f}%)"
            )
        elif text_live_data is not None:
            print_status("------------------------------------------------------")
            print_status(f"[{idir}/{len(scan_dirs)}] Resuming directory from text live report: {top_dir_name}")
            print_status(f"EOS path: {eos_dir}")
            print_status(f"Text live report: {paths['live']}")
            print_status("Finding ROOT files recursively with xrdfs ls -R ...")

            root_files = eos_find_root_files_in_dir(eos_dir)
            n_root_files = len(root_files)
            skip_count = min(text_live_data["files_scanned"], n_root_files)

            processed_paths = root_files[:skip_count]
            processed_set = set(processed_paths)
            dir_summary = summary_from_maybe_dict(text_live_data.get("summary", {}))
            dir_sections = defaultdict(make_empty_summary)
            dir_bad = list(text_live_data.get("bad_results", []))
            resumed.append(top_dir_name)

            top_dirs[top_dir_name] = dir_summary
            merge_summary(total, dir_summary)
            bad.extend(dir_bad)

            # Convert the old text report into the newer JSON checkpoint format
            # immediately. Future resumes will use exact processed_paths instead
            # of relying on positional skipping.
            write_topdir_checkpoint(
                paths["live"],
                "RUNNING",
                top_dir_name,
                eos_dir,
                root_files,
                processed_paths,
                dir_summary,
                dir_sections,
                dir_bad,
                time.monotonic() - dir_start_time,
            )

            pct = 100.0 * skip_count / n_root_files if n_root_files else 0.0
            print_status(
                f"Recovered FILES_SCANNED={text_live_data['files_scanned']}; "
                f"skipping first {skip_count}/{n_root_files} files ({pct:.1f}%)."
            )
            if text_live_data.get("files_total_discovered") and text_live_data["files_total_discovered"] != n_root_files:
                print_status(
                    f"NOTE: old report saw FILES_TOTAL_DISCOVERED={text_live_data['files_total_discovered']}, "
                    f"fresh xrdfs discovery found {n_root_files}. Positional resume assumes ordering is stable."
                )
        else:
            processed_set = set()
            dir_summary = top_dirs[top_dir_name]

            print_status("------------------------------------------------------")
            print_status(f"[{idir}/{len(scan_dirs)}] Starting directory: {top_dir_name}")
            print_status(f"EOS path: {eos_dir}")

            if args.skip_completed and os.path.exists(paths["live"]):
                print_status(
                    f"Found existing .live.txt but it is neither a compatible JSON checkpoint nor a usable text live report; restarting this directory: {paths['live']}"
                )

            print_status("Finding ROOT files recursively with xrdfs ls -R ...")
            root_files = eos_find_root_files_in_dir(eos_dir)
            n_root_files = len(root_files)

            write_topdir_checkpoint(
                paths["live"],
                "RUNNING_DISCOVERY_DONE",
                top_dir_name,
                eos_dir,
                root_files,
                processed_paths,
                dir_summary,
                dir_sections,
                dir_bad,
                time.monotonic() - dir_start_time,
            )

        n_root_files = len(root_files)
        print_status(f"Found {n_root_files} ROOT files in {top_dir_name}")

        if n_root_files == 0:
            elapsed = time.monotonic() - dir_start_time
            write_topdir_checkpoint(
                paths["done"],
                "DONE_NO_ROOT_FILES",
                top_dir_name,
                eos_dir,
                root_files,
                processed_paths,
                dir_summary,
                dir_sections,
                dir_bad,
                elapsed,
            )
            try:
                os.remove(paths["live"])
            except FileNotFoundError:
                pass
            write_global_reports(args.bad_out, args.section_out, args.topdir_out, bad, sections, top_dirs)
            print_status(f"[{idir}/{len(scan_dirs)}] Finished directory: {top_dir_name} -- no ROOT files found")
            continue

        last_progress_time = time.monotonic()
        last_checkpoint_time = time.monotonic()

        for root_path in root_files:
            if root_path in processed_set:
                continue

            section = get_scan_section(root_path, eos_dir)
            res = check_photon_pixel_seed(root_path, args.max_entries)

            add_result_to_summary(total, res)
            add_result_to_summary(sections[section], res)
            add_result_to_summary(dir_sections[section], res)
            add_result_to_summary(top_dirs[top_dir_name], res)

            if not res["ok"]:
                bad_res = {
                    "section": section,
                    **res,
                }
                bad.append(bad_res)
                dir_bad.append(bad_res)

            processed_paths.append(root_path)
            processed_set.add(root_path)
            files_done = len(processed_paths)

            if should_checkpoint(files_done, n_root_files, last_checkpoint_time, args):
                write_topdir_checkpoint(
                    paths["live"],
                    "RUNNING",
                    top_dir_name,
                    eos_dir,
                    root_files,
                    processed_paths,
                    top_dirs[top_dir_name],
                    dir_sections,
                    dir_bad,
                    time.monotonic() - dir_start_time,
                )
                last_checkpoint_time = time.monotonic()

            if should_print_progress(files_done, n_root_files, last_progress_time, args):
                s = top_dirs[top_dir_name]
                elapsed = time.monotonic() - dir_start_time
                rate = (files_done - len(live_data.get("processed_paths", [])) if live_data else files_done) / elapsed if elapsed > 0 else 0.0
                print_status(
                    f"[{idir}/{len(scan_dirs)}] {top_dir_name}: "
                    f"scanned {files_done}/{n_root_files} files "
                    f"({100.0 * files_done / n_root_files:.1f}%), "
                    f"OK={s['ok']}, BAD={s['bad']}, "
                    f"elapsed={format_elapsed(elapsed)}, rate={rate:.2f} new files/s"
                )
                last_progress_time = time.monotonic()

        s = top_dirs[top_dir_name]
        dir_elapsed = time.monotonic() - dir_start_time
        write_topdir_checkpoint(
            paths["done"],
            "DONE",
            top_dir_name,
            eos_dir,
            root_files,
            processed_paths,
            s,
            dir_sections,
            dir_bad,
            dir_elapsed,
        )
        try:
            os.remove(paths["live"])
        except FileNotFoundError:
            pass

        write_global_reports(args.bad_out, args.section_out, args.topdir_out, bad, sections, top_dirs)

        print_status(
            f"[{idir}/{len(scan_dirs)}] Finished directory: {top_dir_name} -- "
            f"files={s['files_found']}, OK={s['ok']}, BAD={s['bad']}, "
            f"elapsed={format_elapsed(dir_elapsed)}"
        )
        print_status(f"Wrote: {paths['done']}")
        print_status("")

    scan_elapsed = time.monotonic() - scan_start_time
    write_global_reports(args.bad_out, args.section_out, args.topdir_out, bad, sections, top_dirs)
    print_status(f"Finished this run in {format_elapsed(scan_elapsed)}")
    if skipped:
        print_status(f"Skipped completed directories: {len(skipped)}")
    if resumed:
        print_status(f"Resumed partial directories: {len(resumed)}")

    print()
    print("======================================================")
    print("Photon_pixelSeed EOS scan summary for this run")
    print("======================================================")
    print(f"EOS base:          {EOS_BASE}")
    print(f"Version suffix:    {args.version_suffix}")
    print(f"Directories found: {len(scan_dirs)}")
    print(f"Directories skipped as completed: {len(skipped)}")
    print(f"Directories resumed from live checkpoint: {len(resumed)}")
    print(f"Report directory:  {args.report_dir}")
    print()

    print_summary_block("TOTAL SCANNED/RECOVERED IN THIS RUN", total)

    print("Top-level directory summaries scanned/recovered in this run")
    print("===========================================================")
    for top_dir in sorted(top_dirs):
        s = top_dirs[top_dir]
        print(
            f"{top_dir}: "
            f"files={s['files_found']} "
            f"OK={s['ok']} "
            f"BAD={s['bad']} "
            f"allFalse={s['all_false']} "
            f"missingBranch={s['missing_branch']} "
            f"photons={s['photons']} "
            f"pixelSeedTrue={s['pixel_seed_true']}"
        )

    print()
    print("Section summaries scanned/recovered in this run")
    print("===============================================")
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
    print(f"Per-directory reports: {args.report_dir}/*.done.txt")
    print(f"Bad-file report for this run:       {args.bad_out}")
    print(f"Section summary table for this run: {args.section_out}")
    print(f"Top-dir summary table for this run: {args.topdir_out}")
    print("======================================================")
    print()


if __name__ == "__main__":
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    main()

