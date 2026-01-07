#!/usr/bin/env python3
"""
Merge 8 kucmsntuple filelists for Run2017F into one, while removing duplicates
defined in duplicates_detailed_2017F.txt.

Rule:
- In duplicates_detailed_2017F.txt, each /store/...root block lists one-or-more
  occurrences like:  Met2017F_files_pX.json  (key: NNNN)
- The *first* (pX,key) listed for that block is the original and is kept.
- All subsequent (pX,key) entries for that same block are duplicates and must be removed
  from the merged kucmsntuple lists.

Removal target in kucmsntuple lists:
- Remove the path that contains:
    /kucmsntuple_MET_R17_SVIPM100_pX_v31_MET_AOD_Run2017F/
  and ends with:
    kucmsntuple_KEY.root
  (KEY is the integer from "(key: KEY)" in duplicates_detailed_2017F.txt)
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple


DUP_FILE_RE = re.compile(r"^\s*Met2017F_files_p(?P<p>\d+)\.json\s*\(key:\s*(?P<key>\d+)\)\s*$")
STORE_ROOT_RE = re.compile(r"^\s*/store/.*\.root\s*$")


def extract_root_path_token(line: str) -> str | None:
    """
    Extract the first token in the line that contains '.root'.
    Keeps any protocol prefix (root://, etc.) if present in the token.
    """
    parts = line.strip().split()
    for tok in parts:
        if ".root" in tok:
            # Strip trailing punctuation that sometimes appears
            return tok.strip(",;")
    return None


def parse_duplicates(dup_path: Path) -> Set[Tuple[int, int]]:
    """
    Returns a set of (p, key) pairs to DROP (i.e., all but the first entry per /store/... block).
    """
    to_drop: Set[Tuple[int, int]] = set()

    current_block: List[Tuple[int, int]] = []
    in_block = False

    def flush_block():
        nonlocal current_block
        if len(current_block) > 1:
            # keep first, drop the rest
            for p, key in current_block[1:]:
                to_drop.add((p, key))
        current_block = []

    with dup_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue

            if STORE_ROOT_RE.match(line):
                # New /store/... root begins -> flush previous block, start new
                flush_block()
                in_block = True
                continue

            m = DUP_FILE_RE.match(line)
            if m:
                in_block = True
                p = int(m.group("p"))
                key = int(m.group("key"))
                current_block.append((p, key))
                continue

            # Any other line: ignore but do not necessarily end a block.
            # Blocks are separated by the next /store/... line (handled above).

    flush_block()
    return to_drop


def build_drop_matchers(sample_tag: str, pairs_to_drop: Set[Tuple[int, int]]) -> List[re.Pattern]:
    """
    Build regex patterns that match the exact kucmsntuple file to drop, per (p,key).

    We match any path containing the dataset tag segment, and ending with kucmsntuple_<key>.root.
    """
    pats: List[re.Pattern] = []
    for p, key in sorted(pairs_to_drop):
        # Example segment:
        #   /kucmsntuple_MET_R17_SVIPM100_p7_v31_MET_AOD_Run2017F/
        seg = re.escape(f"/kucmsntuple_{sample_tag}_p{p}_")
        # Anchor on file name
        fname = re.escape(f"kucmsntuple_{key}.root")
        # Allow anything between segment and filename, but require filename at end of token
        pats.append(re.compile(seg + r".*" + fname + r"$"))
    return pats


def should_drop(path_token: str, drop_patterns: List[re.Pattern]) -> bool:
    # Strip quotes if any
    tok = path_token.strip().strip('"').strip("'")
    return any(p.search(tok) for p in drop_patterns)


def iter_input_lists(paths: Iterable[Path]) -> Iterable[Tuple[Path, int, str]]:
    """
    Yields (file_path, line_no, root_token) for each token that looks like a .root path.
    """
    for p in paths:
        with p.open("r", encoding="utf-8", errors="replace") as f:
            for i, line in enumerate(f, start=1):
                tok = extract_root_path_token(line)
                if tok:
                    yield (p, i, tok)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Merge kucmsntuple filelists and remove duplicates defined by duplicates_detailed_2017F.txt."
    )
    ap.add_argument(
        "--duplicates",
        default="duplicates_detailed_2017F.txt",
        help="Duplicates detail file (default: duplicates_detailed_2017F.txt)",
    )
    ap.add_argument(
        "--in-prefix",
        default="kucmsntuple_MET_R17_SVIPM100_p",
        help="Input list file prefix (default: kucmsntuple_MET_R17_SVIPM100_p)",
    )
    ap.add_argument(
        "--in-suffix",
        default="_v31_MET_AOD_Run2017F.txt",
        help="Input list file suffix (default: _v31_MET_AOD_Run2017F.txt)",
    )
    ap.add_argument(
        "--parts",
        type=int,
        default=8,
        help="Number of parts p1..pN (default: 8)",
    )
    ap.add_argument(
        "--sample-tag",
        default="MET_R17_SVIPM100",
        help="Tag used inside path segment '/kucmsntuple_<sample-tag>_pX_...'(default: MET_R17_SVIPM100)",
    )
    ap.add_argument(
        "--out",
        default="kucmsntuple_MET_R17_SVIPM100_v31_MET_AOD_Run2017F_merged_nodup.txt",
        help="Output merged file name",
    )
    ap.add_argument(
        "--report",
        default="kucmsntuple_MET_R17_SVIPM100_v31_MET_AOD_Run2017F_merged_nodup.report.txt",
        help="Write a small report here (default: ...report.txt)",
    )
    args = ap.parse_args()

    dup_path = Path(args.duplicates)
    if not dup_path.exists():
        print(f"ERROR: duplicates file not found: {dup_path}", file=sys.stderr)
        return 2

    input_lists: List[Path] = []
    for pidx in range(1, args.parts + 1):
        fp = Path(f"{args.in_prefix}{pidx}{args.in_suffix}")
        if not fp.exists():
            print(f"ERROR: input list not found: {fp}", file=sys.stderr)
            return 2
        input_lists.append(fp)

    pairs_to_drop = parse_duplicates(dup_path)
    drop_patterns = build_drop_matchers(args.sample_tag, pairs_to_drop)

    seen: Set[str] = set()
    kept: List[str] = []
    dropped: List[str] = []
    dup_across_inputs: List[str] = []

    for src_file, line_no, tok in iter_input_lists(input_lists):
        clean_tok = tok.strip().strip('"').strip("'")

        if should_drop(clean_tok, drop_patterns):
            dropped.append(f"{src_file}:{line_no}  {clean_tok}")
            continue

        if clean_tok in seen:
            dup_across_inputs.append(f"{src_file}:{line_no}  {clean_tok}")
            continue

        seen.add(clean_tok)
        kept.append(clean_tok)

    out_path = Path(args.out)
    out_path.write_text("\n".join(kept) + ("\n" if kept else ""), encoding="utf-8")

    report_lines = []
    report_lines.append(f"duplicates file: {dup_path}")
    report_lines.append(f"input lists: {', '.join(str(p) for p in input_lists)}")
    report_lines.append(f"output: {out_path}")
    report_lines.append("")
    report_lines.append(f"(p,key) pairs to DROP from duplicates_detailed file: {len(pairs_to_drop)}")
    report_lines.append(f"lines dropped due to duplicates_detailed rules: {len(dropped)}")
    report_lines.append(f"lines dropped due to duplicate paths across merged inputs: {len(dup_across_inputs)}")
    report_lines.append(f"final kept paths: {len(kept)}")
    report_lines.append("")
    if dropped:
        report_lines.append("=== Dropped due to duplicates_detailed rules ===")
        report_lines.extend(dropped[:2000])  # keep report bounded
        if len(dropped) > 2000:
            report_lines.append(f"... truncated ({len(dropped)-2000} more)")
        report_lines.append("")
    if dup_across_inputs:
        report_lines.append("=== Dropped because already seen in merged list (exact path match) ===")
        report_lines.extend(dup_across_inputs[:2000])
        if len(dup_across_inputs) > 2000:
            report_lines.append(f"... truncated ({len(dup_across_inputs)-2000} more)")

    Path(args.report).write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(f"Wrote merged list: {out_path}  (kept {len(kept)} paths)")
    print(f"Wrote report: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

