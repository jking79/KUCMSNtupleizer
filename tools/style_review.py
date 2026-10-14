#!/usr/bin/env python3
"""Audit and mechanically normalize formatting hygiene in tracked source files.

The built-in cleanup is deliberately narrow: it normalizes line endings,
removes ordinary trailing horizontal whitespace, and adds a final newline. It
never changes indentation. Review the diff because whitespace can be meaningful
inside multiline strings or heredocs. Trailing whitespace after a backslash is
reported but not fixed because removing it can activate line continuation.

External language formatters are opt-in with --formatters. They run through
stdin and are never allowed to edit the working tree directly.
"""

from __future__ import annotations

import argparse
import collections
import difflib
import fnmatch
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import tempfile
from typing import Iterable, NamedTuple, Optional, Sequence


CONFIG_NAME = ".style-review.json"
DEFAULT_MAX_PATHS = 30

CPP_SUFFIXES = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".chh"}
PYTHON_SUFFIXES = {".py"}
SHELL_SUFFIXES = {".sh"}


class ReviewError(RuntimeError):
    """A user-facing configuration or repository error."""


class FileReview(NamedTuple):
    path: str
    original: bytes
    cleaned: bytes
    issues: tuple[str, ...]

    @property
    def changed(self) -> bool:
        return self.original != self.cleaned


def run_git(root: Path, args: Sequence[str], *, check: bool = True) -> bytes:
    proc = subprocess.run(
        ["git", "-C", str(root), *args],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if check and proc.returncode:
        detail = proc.stderr.decode("utf-8", errors="replace").strip()
        raise ReviewError(detail or f"git {' '.join(args)} failed")
    return proc.stdout


def find_repo_root(start: Path) -> Path:
    proc = subprocess.run(
        ["git", "-C", str(start), "rev-parse", "--show-toplevel"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        text=True,
    )
    if proc.returncode:
        raise ReviewError(f"not inside a Git repository: {start}")
    return Path(proc.stdout.strip()).resolve()


def load_config(root: Path, config_path: Optional[Path]) -> dict:
    path = config_path or root / CONFIG_NAME
    try:
        config = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ReviewError(f"configuration not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ReviewError(f"invalid JSON in {path}: {exc}") from exc

    if config.get("version") != 1:
        raise ReviewError(f"{path}: expected version 1")
    for key in ("include", "exclude"):
        if not isinstance(config.get(key), list) or not all(
            isinstance(item, str) for item in config[key]
        ):
            raise ReviewError(f"{path}: {key!r} must be a list of strings")
    return config


def matches_any(path: str, patterns: Iterable[str]) -> bool:
    return any(fnmatch.fnmatchcase(path, pattern) for pattern in patterns)


def tracked_paths(root: Path) -> list[str]:
    raw = run_git(root, ["ls-files", "-z"])
    return sorted(
        item.decode("utf-8", errors="surrogateescape")
        for item in raw.split(b"\0")
        if item
    )


def select_paths(
    paths: Iterable[str], config: dict, requested: Sequence[str]
) -> tuple[list[str], list[str]]:
    candidates: list[str] = []
    excluded: list[str] = []
    for path in paths:
        selected = matches_any(path, config["include"])
        if requested:
            selected = selected and any(
                path == item or path.startswith(item.rstrip("/") + "/")
                for item in requested
            )
        if selected and not matches_any(path, config["exclude"]):
            candidates.append(path)
        else:
            excluded.append(path)
    return candidates, excluded


def is_probably_binary(data: bytes) -> bool:
    return b"\0" in data[:8192]


def normalize_bytes(data: bytes) -> tuple[bytes, tuple[str, ...]]:
    issues: list[str] = []
    if is_probably_binary(data):
        return data, ("binary",)

    if b"\r\n" in data:
        issues.append("crlf")
    residual_cr = data.replace(b"\r\n", b"")
    if b"\r" in residual_cr:
        issues.append("bare-cr")

    normalized = data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    had_final_newline = not normalized or normalized.endswith(b"\n")
    lines = normalized.split(b"\n")
    cleaned_lines: list[bytes] = []
    trailing = False
    guarded = False

    for index, line in enumerate(lines):
        is_terminal_empty = (
            index == len(lines) - 1 and line == b"" and normalized.endswith(b"\n")
        )
        if is_terminal_empty:
            cleaned_lines.append(line)
            continue
        stripped = line.rstrip(b" \t")
        if stripped != line:
            if stripped.endswith(b"\\"):
                guarded = True
                cleaned_lines.append(line)
            else:
                trailing = True
                cleaned_lines.append(stripped)
        else:
            cleaned_lines.append(line)

    cleaned = b"\n".join(cleaned_lines)
    if cleaned and not cleaned.endswith(b"\n"):
        cleaned += b"\n"

    if trailing:
        issues.append("trailing-whitespace")
    if guarded:
        issues.append("guarded-backslash-whitespace")
    if not had_final_newline and normalized:
        issues.append("missing-final-newline")
    return cleaned, tuple(issues)


def formatter_kind(path: str) -> Optional[str]:
    source = Path(path)
    suffix = source.suffix.lower()
    if suffix in CPP_SUFFIXES:
        return "clang-format"
    if suffix in PYTHON_SUFFIXES:
        return "black"
    if suffix in SHELL_SUFFIXES:
        return "shfmt"
    return None


def formatter_status(paths: Iterable[str]) -> dict[str, dict[str, object]]:
    counts = collections.Counter(
        kind for path in paths if (kind := formatter_kind(path))
    )
    return {
        name: {"files": counts[name], "available": shutil.which(name) is not None}
        for name in sorted(counts)
    }


def run_formatter(path: str, data: bytes) -> bytes:
    kind = formatter_kind(path)
    if kind is None:
        return data
    executable = shutil.which(kind)
    if executable is None:
        raise ReviewError(f"{kind} is required for --formatters but was not found")

    if kind == "clang-format":
        command = [executable, f"--assume-filename={path}"]
    elif kind == "black":
        command = [executable, "--quiet", "--stdin-filename", path, "-"]
    elif kind == "shfmt":
        command = [executable, "-filename", path]
    else:  # pragma: no cover - formatter_kind controls this branch
        return data

    proc = subprocess.run(
        command, input=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False
    )
    if proc.returncode:
        detail = proc.stderr.decode("utf-8", errors="replace").strip()
        raise ReviewError(f"{kind} failed for {path}: {detail}")
    return proc.stdout


def review_file(root: Path, path: str, with_formatters: bool) -> FileReview:
    absolute = root / path
    if absolute.is_symlink():
        data = absolute.read_bytes()
        return FileReview(path, data, data, ("symlink",))
    original = absolute.read_bytes()
    cleaned, issues = normalize_bytes(original)
    if "binary" not in issues:
        try:
            cleaned.decode("utf-8")
        except UnicodeDecodeError:
            return FileReview(path, original, original, ("non-utf8",))
    if with_formatters and "binary" not in issues:
        formatted = run_formatter(path, cleaned)
        if formatted != cleaned:
            issues = (*issues, f"{formatter_kind(path)}-difference")
        cleaned = formatted
    return FileReview(path, original, cleaned, issues)


def review_files(
    root: Path, paths: Iterable[str], with_formatters: bool
) -> list[FileReview]:
    return [review_file(root, path, with_formatters) for path in paths]


def indentation_profile(root: Path, paths: Iterable[str]) -> dict[str, int]:
    result: collections.Counter[str] = collections.Counter()
    for path in paths:
        if formatter_kind(path) not in {"clang-format", "black", "shfmt"}:
            continue
        data = (root / path).read_bytes()
        if is_probably_binary(data):
            continue
        for line in data.splitlines():
            prefix = line[: len(line) - len(line.lstrip(b" \t"))]
            if not prefix:
                continue
            if b" " in prefix and b"\t" in prefix:
                result["mixed-leading-whitespace-lines"] += 1
            elif b"\t" in prefix:
                result["tab-indented-lines"] += 1
            else:
                result["space-indented-lines"] += 1
    return dict(sorted(result.items()))


def unified_diff(review: FileReview) -> str:
    before = review.original.decode("utf-8", errors="replace").splitlines(keepends=True)
    after = review.cleaned.decode("utf-8", errors="replace").splitlines(keepends=True)
    return "".join(
        difflib.unified_diff(
            before, after, fromfile=f"a/{review.path}", tofile=f"b/{review.path}"
        )
    )


def write_atomic(path: Path, data: bytes) -> None:
    original_mode = stat.S_IMODE(path.stat().st_mode)
    with tempfile.NamedTemporaryFile(
        dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        temporary = Path(handle.name)
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    try:
        temporary.chmod(original_mode)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def print_group(title: str, items: Iterable[str], limit: int) -> None:
    values = list(items)
    print(f"{title}: {len(values)}")
    for value in values[:limit]:
        print(f"  {value}")
    if len(values) > limit:
        print(f"  ... {len(values) - limit} more")


def command_survey(
    root: Path,
    all_paths: list[str],
    candidates: list[str],
    excluded: list[str],
    args: argparse.Namespace,
) -> int:
    reviews = review_files(root, candidates, with_formatters=False)
    issue_counts = collections.Counter(
        issue for review in reviews for issue in review.issues
    )
    changed = [review.path for review in reviews if review.changed]
    issue_paths = [
        f"{review.path}: {', '.join(review.issues)}"
        for review in reviews
        if review.issues
    ]
    payload = {
        "repository": str(root),
        "tracked_files": len(all_paths),
        "cleanup_candidates": len(candidates),
        "excluded_files": len(excluded),
        "proposed_cleanup_changes": len(changed),
        "issue_counts": dict(sorted(issue_counts.items())),
        "indentation": indentation_profile(root, candidates),
        "formatters": formatter_status(candidates),
        "changed_paths": changed,
        "issue_paths": issue_paths,
    }
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    print(f"Repository: {root}")
    print(f"Tracked files: {len(all_paths)}")
    print(f"Maintained-source candidates: {len(candidates)}")
    print(f"Excluded from cleanup: {len(excluded)}")
    print(f"Files with proposed cleanup changes: {len(changed)}")
    print("\nIssue counts:")
    if issue_counts:
        for issue, count in sorted(issue_counts.items()):
            print(f"  {issue}: {count}")
    else:
        print("  none")
    print("\nIndentation inventory (reported only; never auto-fixed):")
    profile = payload["indentation"]
    if profile:
        for label, count in profile.items():
            print(f"  {label}: {count}")
    else:
        print("  no indented source lines")
    print("\nExternal formatter readiness:")
    for name, status in payload["formatters"].items():
        state = "available" if status["available"] else "missing"
        print(f"  {name}: {state} ({status['files']} files)")
    if not payload["formatters"]:
        print("  no formatter-managed files")
    print()
    print_group("Files with issues", issue_paths, args.max_paths)
    return 0


def command_diff(reviews: list[FileReview]) -> int:
    changes = [review for review in reviews if review.changed]
    for review in changes:
        sys.stdout.write(unified_diff(review))
    if not changes:
        print("No cleanup differences.")
    return 1 if changes else 0


def command_check(reviews: list[FileReview], max_paths: int) -> int:
    changes = [review.path for review in reviews if review.changed]
    guarded = [
        review.path
        for review in reviews
        if "guarded-backslash-whitespace" in review.issues
    ]
    if changes:
        print_group("Files requiring cleanup", changes, max_paths)
    if guarded:
        print_group(
            "Manual review required for whitespace after backslash", guarded, max_paths
        )
    if changes or guarded:
        return 1
    print("Style hygiene check passed.")
    return 0


def require_clean_worktree(root: Path) -> None:
    dirty = run_git(root, ["status", "--porcelain", "--untracked-files=no"])
    if dirty:
        raise ReviewError(
            "refusing to apply changes to a dirty tracked worktree; use --allow-dirty to override"
        )


def command_apply(
    root: Path, reviews: list[FileReview], allow_dirty: bool, max_paths: int
) -> int:
    if not allow_dirty:
        require_clean_worktree(root)
    changes = [review for review in reviews if review.changed]
    guarded = [
        review.path
        for review in reviews
        if "guarded-backslash-whitespace" in review.issues
    ]
    for review in changes:
        write_atomic(root / review.path, review.cleaned)
    print_group("Cleaned files", (review.path for review in changes), max_paths)
    if guarded:
        print_group("Left unchanged for manual backslash review", guarded, max_paths)
    diff_check = subprocess.run(
        ["git", "-C", str(root), "diff", "--check"], check=False
    )
    if diff_check.returncode:
        print("git diff --check reported remaining whitespace errors.", file=sys.stderr)
        return diff_check.returncode
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo",
        type=Path,
        default=Path.cwd(),
        help="repository path (default: current directory)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        help=f"configuration path (default: repository/{CONFIG_NAME})",
    )
    parser.add_argument(
        "--path",
        dest="paths",
        action="append",
        default=[],
        metavar="PATH",
        help="restrict to a repository-relative file or directory; repeat as needed",
    )
    parser.add_argument(
        "--max-paths",
        type=int,
        default=DEFAULT_MAX_PATHS,
        help="maximum paths printed per report section",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    survey = subparsers.add_parser(
        "survey", help="report current style hygiene without changing files"
    )
    survey.add_argument(
        "--json", action="store_true", help="emit machine-readable JSON"
    )

    for name, help_text in (
        ("diff", "show the cleanup patch without changing files"),
        ("check", "fail if cleanup or manual review is required"),
        ("apply", "apply cleanup to selected files"),
    ):
        command = subparsers.add_parser(name, help=help_text)
        command.add_argument(
            "--formatters",
            action="store_true",
            help="also run available language formatters",
        )
        if name == "apply":
            command.add_argument(
                "--allow-dirty",
                action="store_true",
                help="permit changes when tracked files are already modified",
            )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        root = find_repo_root(args.repo.resolve())
        config_path = args.config.resolve() if args.config else None
        config = load_config(root, config_path)
        all_paths = tracked_paths(root)
        candidates, excluded = select_paths(all_paths, config, args.paths)
        if not candidates:
            raise ReviewError("no tracked files matched the configured scope")

        if args.command == "survey":
            return command_survey(root, all_paths, candidates, excluded, args)

        reviews = review_files(root, candidates, with_formatters=args.formatters)
        if args.command == "diff":
            return command_diff(reviews)
        if args.command == "check":
            return command_check(reviews, args.max_paths)
        if args.command == "apply":
            return command_apply(root, reviews, args.allow_dirty, args.max_paths)
        raise ReviewError(f"unsupported command: {args.command}")
    except ReviewError as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
