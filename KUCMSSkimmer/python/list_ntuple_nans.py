#!/usr/bin/env python3
"""
List ntuple branches containing NaN values.

This is a one-campaign companion to diagnose_ntuple_diffs.py.  It uses the
same EOS discovery, branch selection, chunked uproot iteration, and optional
golden-lumi masking, but reports absolute NaN rates instead of old/new deltas.
"""

import argparse
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

from diagnose_ntuple_diffs import (
    DEFAULT_EXCLUDE,
    EVENT_TREE,
    LUMI_BRANCH_CANDIDATES,
    RUN_BRANCH_CANDIDATES,
    WEIGHT_BRANCH_CANDIDATES,
    _UPROOT_LOCK,
    branch_allowed,
    branch_manifest,
    chunk_lumi_mask,
    chunks,
    count_lumi_sections,
    discover,
    is_vector_type,
    load_golden_lumi,
    normalize_eos_path,
    pick_first_branch,
    remap_groups,
    tree_specs,
    type_family,
)


def clean_float_values(values):
    return np.asarray(values).astype(np.float64, copy=False)


class NanStats:
    def __init__(self, name, typename):
        self.name = name
        self.typename = typename
        self.family = type_family(typename)
        self.is_vector = is_vector_type(typename)
        self.chunks = 0
        self.events = 0
        self.elements = 0
        self.nan = 0
        self.events_with_nan = 0
        self.empty_events = 0
        self.size_sum = 0.0
        self.size_min = None
        self.size_max = None
        self.errors = []

    def update(self, arr):
        self.chunks += 1
        self.events += len(arr)

        collection_like = self.is_vector
        try:
            sizes = ak.to_numpy(ak.num(arr, axis=1)).astype(np.float64, copy=False)
            collection_like = True
            self.is_vector = True
            if sizes.size:
                self.size_sum += float(np.sum(sizes))
                self.empty_events += int(np.count_nonzero(sizes == 0))
                cur_min = float(np.min(sizes))
                cur_max = float(np.max(sizes))
                self.size_min = cur_min if self.size_min is None else min(self.size_min, cur_min)
                self.size_max = cur_max if self.size_max is None else max(self.size_max, cur_max)
        except Exception:
            collection_like = False

        if self.family not in {"float", "int"}:
            return

        try:
            flat_values = clean_float_values(ak.to_numpy(ak.flatten(arr, axis=None)))
        except Exception as exc:
            self.errors.append(f"flatten: {exc}")
            return

        self.elements += int(flat_values.size)
        if flat_values.size == 0:
            return

        self.nan += int(np.count_nonzero(np.isnan(flat_values)))

        try:
            if collection_like:
                per_event_nan = ak.any(np.isnan(arr), axis=1)
                self.events_with_nan += int(np.count_nonzero(ak.to_numpy(per_event_nan)))
            else:
                event_values = clean_float_values(ak.to_numpy(arr))
                self.events_with_nan += int(np.count_nonzero(np.isnan(event_values)))
        except Exception as exc:
            self.errors.append(f"event nan count: {exc}")

    def as_dict(self):
        nan_fraction = self.nan / self.elements if self.elements else 0.0
        event_fraction = self.events_with_nan / self.events if self.events else 0.0
        size_mean = self.size_sum / self.events if self.events else None
        nan_per_event = self.nan / self.events if self.events else 0.0
        return {
            "type": self.typename,
            "family": self.family,
            "is_vector": self.is_vector,
            "chunks": self.chunks,
            "events": self.events,
            "elements": self.elements,
            "nan": self.nan,
            "nan_fraction": nan_fraction,
            "nan_per_event": nan_per_event,
            "events_with_nan": self.events_with_nan,
            "event_fraction": event_fraction,
            "empty_events": self.empty_events if self.is_vector else None,
            "_size_sum": self.size_sum if self.is_vector else None,
            "size_mean": size_mean if self.is_vector else None,
            "size_min": self.size_min if self.is_vector else None,
            "size_max": self.size_max if self.is_vector else None,
            "errors": self.errors[:5],
        }


def summarize_nan_branches(
    files,
    manifest,
    branches,
    step_size,
    branch_batch,
    jobs=1,
    golden_lumi=None,
    run_branch=None,
    lumi_branch=None,
    quiet=False,
):
    branch_groups = list(chunks(branches, branch_batch))
    if jobs > 1 and len(files) > 1:
        file_groups = [group for group in chunks(files, max(1, math.ceil(len(files) / jobs))) if group]
        try:
            worker_stats = []
            errors = []
            total_events = 0
            with ProcessPoolExecutor(max_workers=jobs) as executor:
                futures = {
                    executor.submit(
                        summarize_nan_file_group,
                        file_group,
                        manifest,
                        branches,
                        branch_batch,
                        step_size,
                        golden_lumi,
                        run_branch,
                        lumi_branch,
                    ): (idx, file_group)
                    for idx, file_group in enumerate(file_groups, start=1)
                }
                for future in as_completed(futures):
                    idx, file_group = futures[future]
                    if not quiet:
                        print(
                            f"  file shard {idx}: {len(file_group)} file(s) done",
                            file=sys.stderr,
                        )
                    try:
                        group_stats, group_events, group_errors = future.result()
                        worker_stats.append(group_stats)
                        total_events += group_events
                        errors.extend(group_errors)
                    except Exception as exc:
                        errors.append({"files": file_group[:5], "error": str(exc)})
            return combine_stats_by_branch(worker_stats), total_events, errors
        except OSError as exc:
            if not quiet:
                print(
                    f"WARNING: --jobs {jobs} unavailable ({exc}); falling back to serial scan",
                    file=sys.stderr,
                )

    stats = {name: NanStats(name, manifest[name]) for name in branches}
    errors = []
    total_events = 0
    specs = tree_specs(files, EVENT_TREE)
    count_events_group = branch_groups[0] if branch_groups else []
    for idx, group in enumerate(branch_groups, start=1):
        count_events = group == count_events_group
        if not quiet:
            print(f"  branch batch {idx}: {len(group)} branch(es)", file=sys.stderr)
        try:
            expressions = list(group)
            if golden_lumi:
                for name in (run_branch, lumi_branch):
                    if name and name not in expressions:
                        expressions.append(name)
            with _UPROOT_LOCK:
                for arrays in uproot.iterate(
                    specs,
                    expressions=expressions,
                    step_size=step_size,
                    library="ak",
                ):
                    keep = None
                    if golden_lumi:
                        keep = chunk_lumi_mask(arrays, run_branch, lumi_branch, golden_lumi)
                    if count_events:
                        total_events += int(np.count_nonzero(keep)) if keep is not None else len(arrays[group[0]])
                    for name in group:
                        if name not in arrays.fields:
                            continue
                        arr = arrays[name]
                        if keep is not None:
                            arr = arr[keep]
                        stats[name].update(arr)
        except Exception as exc:
            errors.append({"branches": group, "error": str(exc)})

    return strip_internal_stats({name: obj.as_dict() for name, obj in stats.items()}), total_events, errors


def summarize_nan_file_group(
    files,
    manifest,
    branches,
    branch_batch,
    step_size,
    golden_lumi=None,
    run_branch=None,
    lumi_branch=None,
):
    stats = {name: NanStats(name, manifest[name]) for name in branches}
    errors = []
    total_events = 0
    specs = tree_specs(files, EVENT_TREE)
    branch_groups = list(chunks(branches, branch_batch))
    count_events_group = branch_groups[0] if branch_groups else []
    for group in branch_groups:
        count_events = group == count_events_group
        try:
            expressions = list(group)
            if golden_lumi:
                for name in (run_branch, lumi_branch):
                    if name and name not in expressions:
                        expressions.append(name)
            with _UPROOT_LOCK:
                for arrays in uproot.iterate(
                    specs,
                    expressions=expressions,
                    step_size=step_size,
                    library="ak",
                ):
                    keep = None
                    if golden_lumi:
                        keep = chunk_lumi_mask(arrays, run_branch, lumi_branch, golden_lumi)
                    if count_events:
                        total_events += int(np.count_nonzero(keep)) if keep is not None else len(arrays[group[0]])
                    for name in group:
                        if name not in arrays.fields:
                            continue
                        arr = arrays[name]
                        if keep is not None:
                            arr = arr[keep]
                        stats[name].update(arr)
        except Exception as exc:
            errors.append({"branches": group, "error": str(exc)})
    return {name: obj.as_dict() for name, obj in stats.items()}, total_events, errors


def combine_stats_by_branch(worker_stats):
    combined = {}
    for stats in worker_stats:
        for name, payload in stats.items():
            if name not in combined:
                combined[name] = dict(payload)
                continue
            target = combined[name]
            for key in ["chunks", "events", "elements", "nan", "events_with_nan", "empty_events"]:
                target[key] = (target.get(key) or 0) + (payload.get(key) or 0)
            if payload.get("size_min") is not None:
                target["size_min"] = (
                    payload["size_min"]
                    if target.get("size_min") is None
                    else min(target["size_min"], payload["size_min"])
                )
            if payload.get("size_max") is not None:
                target["size_max"] = (
                    payload["size_max"]
                    if target.get("size_max") is None
                    else max(target["size_max"], payload["size_max"])
                )
            total_events = target.get("events") or 0
            target["nan_fraction"] = target["nan"] / target["elements"] if target.get("elements") else 0.0
            target["nan_per_event"] = target["nan"] / total_events if total_events else 0.0
            target["event_fraction"] = (
                target["events_with_nan"] / total_events if total_events else 0.0
            )
            target["_size_sum"] = (target.get("_size_sum") or 0.0) + (
                payload.get("_size_sum") or 0.0
            )
            target["size_mean"] = None
            if target.get("is_vector") and total_events:
                target["size_mean"] = target["_size_sum"] / total_events
            target["errors"] = (target.get("errors") or []) + (payload.get("errors") or [])
            target["errors"] = target["errors"][:5]
    for payload in combined.values():
        payload.pop("_size_sum", None)
    return combined


def strip_internal_stats(stats):
    for payload in stats.values():
        payload.pop("_size_sum", None)
    return stats


def sort_nan_rows(stats):
    rows = [
        {"name": name, **payload}
        for name, payload in stats.items()
        if payload.get("nan", 0) > 0
    ]
    rows.sort(
        key=lambda item: (
            item.get("events_with_nan", 0),
            item.get("nan", 0),
            item.get("nan_fraction", 0.0),
        ),
        reverse=True,
    )
    return rows


def format_pct(value):
    if value is None:
        return "n/a"
    return f"{100.0 * value:.4g}%"


def format_float(value):
    if value is None:
        return "n/a"
    if isinstance(value, float) and math.isfinite(value):
        return f"{value:.6g}"
    return str(value)


def write_text_report(path, report):
    rows = report["nan_branches"]
    lines = [
        "KUCMS ntuple NaN scan",
        "=" * 72,
        f"Label: {report['label']}",
        f"Input path: {report['input_path']}",
        f"Campaign files found: {report['file_count']}",
        f"Branch-scan files: {report['branch_file_count']}",
        f"Branch-scan events: {report['event_entries']}",
        f"Branches scanned: {report['branches_scanned']}",
        f"Branches with NaNs: {len(rows)}",
        "",
    ]
    if report["branch_file_count"] != report["file_count"]:
        lines.extend(
            [
                "Sampling",
                "-" * 72,
                (
                    "NaN counts and event fractions are for the branch-scan file subset, "
                    "not the full campaign."
                ),
                (
                    "Run without --branch-max-files-per-group and with "
                    "--branch-file-stride 1 for full-campaign NaN rates."
                ),
                "",
            ]
        )
    if report.get("golden_lumi", {}).get("enabled"):
        lumi = report["golden_lumi"]
        lines.extend(
            [
                "Golden lumi",
                "-" * 72,
                f"mask: {lumi['path']}",
                f"runs: {lumi['runs']} lumi sections: {lumi['sections']}",
                f"run branch: {lumi['run_branch']} lumi branch: {lumi['lumi_branch']}",
                "",
            ]
        )

    lines.extend(
        [
            "Branches with NaNs",
            "-" * 72,
            (
                f"{'branch':38s} {'nan':>12s} {'nan/elem':>11s} "
                f"{'evts_nan':>12s} {'evts_nan/all':>13s} {'nan/evt':>10s} "
                f"{'size_mean':>10s}"
            ),
        ]
    )
    for row in rows:
        lines.append(
            f"{row['name'][:38]:38s} "
            f"{row['nan']:12d} "
            f"{format_pct(row['nan_fraction']):>11s} "
            f"{row['events_with_nan']:12d} "
            f"{format_pct(row['event_fraction']):>13s} "
            f"{format_float(row['nan_per_event']):>10s} "
            f"{format_float(row.get('size_mean')):>10s}"
        )

    if report.get("errors", {}).get("branches"):
        lines.extend(["", "Branch read errors", "-" * 72])
        for error in report["errors"]["branches"][:20]:
            lines.append(str(error))

    lines.extend(
        [
            "",
            "Notes",
            "-" * 72,
            "nan/elem is NaN values divided by flattened branch elements.",
            "evts_nan/all is events with at least one NaN in that branch divided by scanned events.",
            "size_mean is the mean collection length per event for vector branches.",
        ]
    )
    Path(path).write_text("\n".join(lines) + "\n")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Scan one KUCMS ntuple campaign and list branches containing NaNs."
    )
    parser.add_argument("path", help="EOS directory, local directory, or ROOT file to scan.")
    parser.add_argument("--label", default="sample", help="Label used in reports.")
    parser.add_argument("--json", default="ntuple_nan_scan.json", help="JSON report path.")
    parser.add_argument("--text", default="ntuple_nan_scan.txt", help="Text report path.")
    parser.add_argument(
        "--golden-lumi",
        default=None,
        help="Golden lumi JSON mask applied to event counts and branch summaries.",
    )
    parser.add_argument("--step-size", default="25 MB", help="uproot step_size, e.g. '10 MB'.")
    parser.add_argument("--branch-batch", type=int, default=6, help="Branches read per pass.")
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Parallel file shards to scan. Values above 1 use separate processes.",
    )
    parser.add_argument("--max-files", type=int, default=-1, help="Limit files per discovered group.")
    parser.add_argument("--max-branches", type=int, default=-1, help="Debug limit on summarized branches.")
    parser.add_argument(
        "--branch-max-files-per-group",
        type=int,
        default=-1,
        help="Limit files per group for event-tree branch summaries.",
    )
    parser.add_argument(
        "--branch-file-stride",
        type=int,
        default=1,
        help="Use every Nth file per group for event-tree branch summaries.",
    )
    parser.add_argument("--include", action="append", default=[], help="Shell pattern for branches to include.")
    parser.add_argument(
        "--exclude",
        action="append",
        default=list(DEFAULT_EXCLUDE),
        help="Shell pattern for branches to exclude. Can be repeated.",
    )
    parser.add_argument("--verbose", action="store_true", help="Print xrdfs commands.")
    parser.add_argument("--quiet", action="store_true", help="Reduce progress messages.")
    args = parser.parse_args()
    if args.branch_batch < 1:
        parser.error("--branch-batch must be >= 1")
    if args.branch_file_stride < 1:
        parser.error("--branch-file-stride must be >= 1")
    if args.jobs < 1:
        parser.error("--jobs must be >= 1")
    args.golden_lumi_mask = None
    if args.golden_lumi:
        args.golden_lumi_mask = load_golden_lumi(args.golden_lumi)
    return args


def main():
    args = parse_args()
    input_path = normalize_eos_path(args.path)
    print(f"Discovering files under {input_path}", file=sys.stderr)
    datasets = discover(input_path, max_files=args.max_files, verbose=args.verbose)
    if not datasets:
        sys.exit("ERROR: no ROOT files found")
    datasets = remap_groups(datasets, "auto")
    all_files = sorted(path for files in datasets.values() for path in files)

    branch_datasets = {}
    for group_name, files in sorted(datasets.items()):
        group_files = sorted(files)
        if args.branch_file_stride > 1:
            group_files = group_files[:: args.branch_file_stride]
        if args.branch_max_files_per_group > 0:
            group_files = group_files[: args.branch_max_files_per_group]
        branch_datasets[group_name] = group_files
    branch_files = sorted(path for files in branch_datasets.values() for path in files)
    if not branch_files:
        branch_files = all_files

    pre_manifest = branch_manifest(branch_files, label=args.label)
    run_branch = None
    lumi_branch = None
    weight_branch = None
    if args.golden_lumi_mask:
        run_branch = pick_first_branch(pre_manifest, RUN_BRANCH_CANDIDATES)
        lumi_branch = pick_first_branch(pre_manifest, LUMI_BRANCH_CANDIDATES)
        weight_branch = pick_first_branch(pre_manifest, WEIGHT_BRANCH_CANDIDATES)
        if not run_branch or not lumi_branch:
            raise RuntimeError(
                "--golden-lumi requires run/lumi branches in "
                f"{EVENT_TREE}; tried run={RUN_BRANCH_CANDIDATES}, "
                f"lumi={LUMI_BRANCH_CANDIDATES}"
            )

    print(
        f"[{args.label}] files: {len(all_files)} branch-summary files: {len(branch_files)}",
        file=sys.stderr,
    )
    selected_branches = sorted(
        name
        for name in pre_manifest
        if branch_allowed(name, args.include, args.exclude)
        and type_family(pre_manifest[name]) == "float"
    )
    if args.max_branches > 0:
        selected_branches = selected_branches[: args.max_branches]

    print(
        f"[{args.label}] scanning {len(selected_branches)} numeric branch(es) "
        f"with {args.jobs} job(s)",
        file=sys.stderr,
    )
    stats, event_entries, branch_errors = summarize_nan_branches(
        branch_files,
        pre_manifest,
        selected_branches,
        args.step_size,
        args.branch_batch,
        jobs=args.jobs,
        golden_lumi=args.golden_lumi_mask,
        run_branch=run_branch,
        lumi_branch=lumi_branch,
        quiet=args.quiet,
    )
    rows = sort_nan_rows(stats)

    report = {
        "label": args.label,
        "input_path": input_path,
        "file_count": len(all_files),
        "branch_file_count": len(branch_files),
        "event_entries": event_entries,
        "groups": {name: len(files) for name, files in sorted(datasets.items())},
        "branch_groups": {name: len(files) for name, files in sorted(branch_datasets.items())},
        "branches_scanned": len(selected_branches),
        "nan_branch_count": len(rows),
        "nan_branches": rows,
        "branch_stats": stats,
        "golden_lumi": {
            "enabled": bool(args.golden_lumi_mask),
            "path": args.golden_lumi,
            "runs": len(args.golden_lumi_mask or {}),
            "sections": count_lumi_sections(args.golden_lumi_mask or {}),
            "run_branch": run_branch,
            "lumi_branch": lumi_branch,
            "weight_branch": weight_branch,
        },
        "settings": {
            "step_size": args.step_size,
            "branch_batch": args.branch_batch,
            "jobs": args.jobs,
            "branch_max_files_per_group": args.branch_max_files_per_group,
            "branch_file_stride": args.branch_file_stride,
            "max_files": args.max_files,
            "max_branches": args.max_branches,
            "include": args.include,
            "exclude": args.exclude,
        },
        "errors": {
            "events": [],
            "branches": branch_errors[:20],
        },
    }
    Path(args.json).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    write_text_report(args.text, report)
    print(f"Wrote JSON report: {args.json}", file=sys.stderr)
    print(f"Wrote text report: {args.text}", file=sys.stderr)


if __name__ == "__main__":
    main()
