#!/usr/bin/env python3
"""Collect and summarize one user's live and historical LPC HTCondor jobs.

Run on a cmslpc interactive node with the LPC HTCondor environment available.
The report is written as Markdown to --output (default: a timestamped file in
 the current directory) and also printed to stdout.
"""

from __future__ import annotations

import argparse
import datetime as dt
import errno
import getpass
import json
import os
import shlex
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

STATUS_NAMES = {
    1: "Idle",
    2: "Running",
    3: "Removed",
    4: "Completed",
    5: "Held",
    6: "Transferring output",
    7: "Suspended",
}


def run_condor(args: list[str], *, timeout: int) -> tuple[list[dict[str, Any]] | None, str | None]:
    """Run a Condor command expected to return JSON ClassAds."""
    try:
        result = subprocess.run(
            args,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )
    except OSError as exc:
        # Some LPC Condor command wrappers are shell-readable scripts that
        # interactive shells can run but execve cannot launch directly.
        if exc.errno != errno.ENOEXEC:
            return None, f"Could not run {args[0]}: {exc}"
        try:
            result = subprocess.run(
                ["bash", "-c", shlex.join(args)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=timeout,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as fallback_exc:
            return None, f"Could not run {args[0]} directly or through bash: {fallback_exc}"
    except subprocess.TimeoutExpired:
        return None, f"Timed out after {timeout}s: {shlex.join(args)}"

    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip() or f"exit code {result.returncode}"
        return None, f"{shlex.join(args)} failed: {detail}"

    raw = result.stdout.strip()
    if not raw:
        return [], None
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError as exc:
        preview = raw[:500].replace("\n", " ")
        return None, f"Could not parse JSON from {args[0]}: {exc}; output begins: {preview!r}"
    if isinstance(parsed, dict):
        parsed = [parsed]
    if not isinstance(parsed, list):
        return None, f"Unexpected JSON format from {args[0]}"
    return [ad for ad in parsed if isinstance(ad, dict)], None


def get_attr(ad: dict[str, Any], *names: str, default: Any = None) -> Any:
    """HTCondor JSON attribute names are usually capitalized; tolerate variants."""
    for name in names:
        for key, value in ad.items():
            if key.lower() == name.lower() and value is not None:
                return value
    return default


def discover_schedds(timeout: int) -> tuple[list[str], str | None]:
    try:
        result = subprocess.run(
            ["condor_status", "-schedd", "-af", "Name"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return [], f"Could not discover schedds: {exc}"
    if result.returncode != 0:
        return [], result.stderr.strip() or "condor_status failed"
    # LPC's pool advertises other schedds and CE endpoints too. Jobs submitted
    # through the LPC farm are hosted by lpcscheddN; querying unrelated
    # services can generate authorization/connection failures and noise.
    names = sorted({
        line.strip() for line in result.stdout.splitlines()
        if line.strip().startswith("lpcschedd") and line.strip().endswith(".fnal.gov")
    })
    return names, None


def timestamp(value: Any) -> str:
    try:
        number = int(value)
        if number <= 0:
            return "—"
        return dt.datetime.fromtimestamp(number).astimezone().strftime("%Y-%m-%d %H:%M %Z")
    except (TypeError, ValueError, OSError, OverflowError):
        return "—"


def elapsed(seconds: Any) -> str:
    try:
        total = max(0, int(float(seconds)))
    except (TypeError, ValueError):
        return "—"
    days, rem = divmod(total, 86400)
    hours, rem = divmod(rem, 3600)
    minutes, secs = divmod(rem, 60)
    if days:
        return f"{days}d {hours:02}:{minutes:02}:{secs:02}"
    return f"{hours:02}:{minutes:02}:{secs:02}"


def job_id(ad: dict[str, Any]) -> str:
    cluster = get_attr(ad, "ClusterId", default="?")
    proc = get_attr(ad, "ProcId", default="?")
    return f"{cluster}.{proc}"


def job_label(ad: dict[str, Any]) -> str:
    return str(get_attr(ad, "JobBatchName", "Cmd", default="(command unavailable)"))


def history_outcome(ad: dict[str, Any]) -> str:
    status = int(get_attr(ad, "JobStatus", default=0) or 0)
    if status == 3:
        return "Removed"
    if status != 4:
        return STATUS_NAMES.get(status, "Unknown")
    signaled = bool(get_attr(ad, "ExitBySignal", default=False))
    signal = get_attr(ad, "ExitSignal")
    code = get_attr(ad, "ExitCode")
    if signaled:
        return f"Failed (signal {signal if signal is not None else '?'})"
    if code is None:
        return "Completed (exit unknown)"
    if int(code) != 0:
        return f"Failed (exit {code})"
    return "Succeeded"


def md(value: Any) -> str:
    return str(value if value is not None else "—").replace("|", "\\|").replace("\n", " ").replace("\r", " ")


def summarize(user: str, schedds: list[str], history_limit: int | None, details: int, timeout: int) -> tuple[str, list[str]]:
    issues: list[str] = []
    live_by_schedd: dict[str, list[dict[str, Any]]] = {}
    history_by_schedd: dict[str, list[dict[str, Any]]] = {}

    for schedd in schedds:
        base = ["-name", schedd]
        q_args = ["condor_q", *base, "-json", user]
        q_ads, q_error = run_condor(q_args, timeout=timeout)
        if q_error:
            issues.append(f"Queue query failed for `{schedd}`: {q_error}")
            live_by_schedd[schedd] = []
        else:
            live_by_schedd[schedd] = q_ads or []

        h_args = ["condor_history", *base, "-json"]
        if history_limit is not None:
            h_args.extend(["-limit", str(history_limit)])
        h_args.append(user)
        h_ads, h_error = run_condor(h_args, timeout=timeout)
        if h_error:
            issues.append(f"History query failed for `{schedd}`: {h_error}")
            history_by_schedd[schedd] = []
        else:
            history_by_schedd[schedd] = h_ads or []

    current_ads = [(schedd, ad) for schedd, ads in live_by_schedd.items() for ad in ads]
    history_ads = [(schedd, ad) for schedd, ads in history_by_schedd.items() for ad in ads]

    # A schedd's history can occasionally include duplicate records. Keep one
    # record per schedd and job ID for counts and report tables.
    unique_history: dict[tuple[str, str], tuple[str, dict[str, Any]]] = {}
    for schedd, ad in history_ads:
        unique_history[(schedd, job_id(ad))] = (schedd, ad)
    history_ads = list(unique_history.values())

    current_counts = Counter(int(get_attr(ad, "JobStatus", default=0) or 0) for _, ad in current_ads)
    history_counts = Counter(history_outcome(ad) for _, ad in history_ads)

    now = dt.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    lines = [
        f"# LPC HTCondor report for `{user}`",
        "",
        f"Generated: {now}  ",
        f"Schedds queried: {len(schedds)}  ",
        f"Current jobs found: {len(current_ads)}  ",
        f"History records found: {len(history_ads)}",
        "",
        "## Summary",
        "",
        "| Current queue status | Jobs |",
        "|---|---:|",
    ]
    for code in (1, 2, 5, 6, 7, 3, 4, 0):
        count = current_counts.get(code, 0)
        if count or code in (1, 2, 5):
            lines.append(f"| {STATUS_NAMES.get(code, 'Unknown')} | {count} |")
    lines.extend([
        "",
        "| History outcome | Jobs |",
        "|---|---:|",
    ])
    preferred_history = ["Succeeded", "Failed", "Completed (exit unknown)", "Removed"]
    for outcome in preferred_history:
        lines.append(f"| {outcome} | {history_counts.get(outcome, 0)} |")
    for outcome, count in sorted(history_counts.items()):
        if outcome not in preferred_history:
            lines.append(f"| {outcome} | {count} |")

    lines.extend(["", "## Per-schedd totals", "", "| Schedd | Idle | Running | Held | Other current | Succeeded | Failed | Removed |", "|---|---:|---:|---:|---:|---:|---:|---:|"])
    for schedd in schedds:
        live = live_by_schedd.get(schedd, [])
        hist = history_by_schedd.get(schedd, [])
        lc = Counter(int(get_attr(ad, "JobStatus", default=0) or 0) for ad in live)
        hc = Counter(history_outcome(ad) for ad in hist)
        other_live = sum(n for status, n in lc.items() if status not in (1, 2, 5))
        failed = sum(n for outcome, n in hc.items() if outcome.startswith("Failed"))
        lines.append(f"| `{schedd}` | {lc[1]} | {lc[2]} | {lc[5]} | {other_live} | {hc['Succeeded']} | {failed} | {hc['Removed']} |")

    attention_live = [(s, a) for s, a in current_ads if int(get_attr(a, "JobStatus", default=0) or 0) == 5]
    attention_hist = [(s, a) for s, a in history_ads if history_outcome(a).startswith("Failed")]
    attention_hist.sort(key=lambda pair: int(get_attr(pair[1], "CompletionDate", "EnteredCurrentStatus", default=0) or 0), reverse=True)

    lines.extend(["", "## Jobs needing attention", ""])
    if not attention_live and not attention_hist:
        lines.append("No held current jobs or failed historical jobs were found in the queried records.")
    if attention_live:
        lines.extend(["### Held jobs", "", "| Job ID | Schedd | Submitted | Hold reason | Command |", "|---|---|---|---|---|"])
        for schedd, ad in attention_live:
            lines.append("| {} | `{}` | {} | {} | {} |".format(
                md(job_id(ad)), md(schedd), md(timestamp(get_attr(ad, "QDate"))),
                md(get_attr(ad, "HoldReason", default="(no reason provided)")), md(job_label(ad))))
        lines.append("")
    if attention_hist:
        lines.extend([f"### Failed historical jobs (latest {min(details, len(attention_hist))})", "", "| Job ID | Schedd | Finished | Outcome | Runtime | Command |", "|---|---|---|---|---:|---|"])
        for schedd, ad in attention_hist[:details]:
            lines.append("| {} | `{}` | {} | {} | {} | {} |".format(
                md(job_id(ad)), md(schedd), md(timestamp(get_attr(ad, "CompletionDate", "EnteredCurrentStatus"))),
                md(history_outcome(ad)), md(elapsed(get_attr(ad, "RemoteWallClockTime"))), md(job_label(ad))))
        lines.append("")

    recent = sorted(history_ads, key=lambda pair: int(get_attr(pair[1], "CompletionDate", "EnteredCurrentStatus", default=0) or 0), reverse=True)
    lines.extend([f"## Latest history records (up to {min(details, len(recent))})", "", "| Job ID | Schedd | Finished | Outcome | Runtime | Command |", "|---|---|---|---|---:|---|"])
    if recent:
        for schedd, ad in recent[:details]:
            lines.append("| {} | `{}` | {} | {} | {} | {} |".format(
                md(job_id(ad)), md(schedd), md(timestamp(get_attr(ad, "CompletionDate", "EnteredCurrentStatus"))),
                md(history_outcome(ad)), md(elapsed(get_attr(ad, "RemoteWallClockTime"))), md(job_label(ad))))
    else:
        lines.append("| — | — | — | No history records returned | — | — |")

    lines.extend(["", "## Coverage notes", ""])
    if history_limit is None:
        lines.append("History was queried without a record limit; the report includes records retained by each schedd's history service.")
    else:
        lines.append(f"History was limited to the most recent {history_limit} records per schedd. Pass `--history-limit 0` to request all retained records.")
    lines.append("Condor history reflects records retained by the schedulers; it may not include jobs whose history has expired or was otherwise purged.")
    if issues:
        lines.extend(["", "## Query issues", ""])
        lines.extend(f"- {issue}" for issue in issues)
    lines.append("")
    return "\n".join(lines), issues


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize your current and historical LPC Condor jobs.")
    parser.add_argument("-u", "--user", default=getpass.getuser(), help="HTCondor owner/login name (default: current user)")
    parser.add_argument("-o", "--output", help="Markdown report path (default: timestamped file in current directory)")
    parser.add_argument("--history-limit", type=int, default=0,
                        help="Maximum history records per schedd; 0 means all retained history (default: 0)")
    parser.add_argument("--details", type=int, default=25, help="Maximum detailed recent/failure rows (default: 25)")
    parser.add_argument("--timeout", type=int, default=120, help="Per-command timeout in seconds (default: 120)")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.history_limit < 0 or args.details < 1 or args.timeout < 1:
        print("ERROR: history-limit must be nonnegative; details and timeout must be positive.", file=sys.stderr)
        return 2
    for command in ("condor_status", "condor_q", "condor_history"):
        if not shutil_which(command):
            print(f"ERROR: {command} is unavailable. Run this from a cmslpc interactive node with Condor configured.", file=sys.stderr)
            return 1

    schedds, error = discover_schedds(args.timeout)
    if error or not schedds:
        print(f"ERROR: Could not discover LPC schedds: {error or 'none found'}", file=sys.stderr)
        return 1

    limit = args.history_limit or None
    report, issues = summarize(args.user, schedds, limit, args.details, args.timeout)
    output = Path(args.output) if args.output else Path.cwd() / f"lpc_condor_report_{dt.datetime.now().astimezone().strftime('%Y%m%d_%H%M%S')}.md"
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(report, encoding="utf-8")
    except OSError as exc:
        print(f"ERROR: Could not write report to {output}: {exc}", file=sys.stderr)
        return 1

    print(report)
    print(f"\nSaved report: {output}")
    return 1 if issues else 0


def shutil_which(command: str) -> str | None:
    """Avoid importing shutil solely for a simple PATH lookup."""
    for directory in os.environ.get("PATH", "").split(os.pathsep):
        candidate = Path(directory) / command
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return None


if __name__ == "__main__":
    raise SystemExit(main())

