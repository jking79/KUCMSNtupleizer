#!/usr/bin/env python3
import json
import re
import sys
from pathlib import Path
from collections import Counter, defaultdict
from typing import Dict, List, Tuple, Optional


EOS_PREFIX = ""  # not used, placeholder


def read_nonblank_lines(path: Path) -> List[str]:
    out = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s:
                out.append(s)
    return out


def basename(p: str) -> str:
    return p.rsplit("/", 1)[-1]


# Expect: .../pX/.../.../kucmsntuple_N.root  (your example matches this)
RE_P_AND_N = re.compile(r"/(p[1-9])/.*/kucmsntuple_(\d+)\.root$")


def parse_processed_line(line: str) -> Optional[Tuple[str, str]]:
    m = RE_P_AND_N.search(line)
    if not m:
        return None
    return (m.group(1), m.group(2))


def load_mapping_json(path: Path) -> Dict[str, str]:
    """
    JSON looks like: {"19": ["/store/data/.../AAA.root"], ...}
    Return dict: N(string) -> source_path(string)
    """
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    out: Dict[str, str] = {}
    for k, v in data.items():
        ks = str(k)
        if isinstance(v, list) and len(v) >= 1 and isinstance(v[0], str):
            out[ks] = v[0]
        elif isinstance(v, str):
            out[ks] = v
        # else: ignore malformed entries
    return out


def main():
    if len(sys.argv) != 5:
        print(
            "Usage:\n"
            f"  {sys.argv[0]} NEEDED.txt PROCESSED.txt JSON_DIR METTAG\n\n"
            "Example:\n"
            f"  {sys.argv[0]} Run2017F_Fullcrabfilelist.txt processed.txt ./ Met2017F_files\n\n"
            "Expects JSONs named: METTAG_p1.json ... METTAG_p9.json\n"
        )
        sys.exit(1)

    needed_txt = Path(sys.argv[1])
    processed_txt = Path(sys.argv[2])
    json_dir = Path(sys.argv[3])
    mettag = sys.argv[4]

    if not needed_txt.exists():
        print(f"ERROR: needed file list not found: {needed_txt}")
        sys.exit(2)
    if not processed_txt.exists():
        print(f"ERROR: processed file list not found: {processed_txt}")
        sys.exit(2)
    if not json_dir.exists():
        print(f"ERROR: json_dir not found: {json_dir}")
        sys.exit(2)

    # Load needed AOD basenames (UUID.root)
    needed_lines = read_nonblank_lines(needed_txt)
    needed_basenames = [basename(x) for x in needed_lines]
    needed_set = set(needed_basenames)
    needed_counts = Counter(needed_basenames)

    # Load processed lines
    processed_lines = read_nonblank_lines(processed_txt)

    # Parse processed lines -> (pX, N)
    parsed: List[Tuple[str, str, str]] = []  # (orig, pX, N)
    unparseable_processed: List[str] = []
    for line in processed_lines:
        pn = parse_processed_line(line)
        if pn is None:
            unparseable_processed.append(line)
        else:
            parsed.append((line, pn[0], pn[1]))

    # Load mappings
    mappings: Dict[str, Dict[str, str]] = {}
    missing_jsons: List[str] = []
    for i in range(1, 10):
        p = f"p{i}"
        jpath = json_dir / f"{mettag}_{p}.json"
        if not jpath.exists():
            missing_jsons.append(str(jpath))
            continue
        mappings[p] = load_mapping_json(jpath)

    # Map processed -> needed AOD basename
    # Keep the FIRST processed line for each mapped AOD basename.
    keep_processed_lines: List[str] = []
    dup_processed_lines: List[str] = []
    unmapped_processed_lines: List[str] = []  # includes parseable-but-unmapped and missing-pX-json
    mapped_aod_for_processed: Dict[str, str] = {}  # processed line -> aod basename
    first_seen_processed_for_aod: Dict[str, str] = {}  # aod basename -> processed line

    # Track mapping failure reasons
    missing_p_json: List[Tuple[str, str, str]] = []
    missing_key_in_json: List[Tuple[str, str, str]] = []

    for orig, p, n in parsed:
        if p not in mappings:
            unmapped_processed_lines.append(orig)
            missing_p_json.append((orig, p, n))
            continue

        mp = mappings[p]
        if n not in mp:
            unmapped_processed_lines.append(orig)
            missing_key_in_json.append((orig, p, n))
            continue

        src = mp[n]
        aod_base = basename(src)
        mapped_aod_for_processed[orig] = aod_base

        if aod_base not in first_seen_processed_for_aod:
            first_seen_processed_for_aod[aod_base] = orig
            keep_processed_lines.append(orig)
        else:
            dup_processed_lines.append(orig)

    # Any unparseable lines are also unmapped outputs
    unmapped_processed_lines.extend(unparseable_processed)

    # Compare needed vs mapped
    mapped_aod_set = set(first_seen_processed_for_aod.keys())  # unique AOD basenames mapped at least once
    missing_needed = sorted(needed_set - mapped_aod_set)
    extra_mapped_not_needed = sorted(mapped_aod_set - needed_set)

    # Summary counts
    n_needed_unique = len(needed_set)
    n_mapped_unique = len(mapped_aod_set)
    n_missing_needed = len(missing_needed)

    n_keep = len(keep_processed_lines)
    n_dups = len(dup_processed_lines)
    n_unmapped = len(unmapped_processed_lines)

    expected = len(needed_set)
    actual = len(keep_processed_lines)

    one_to_one_ok = (n_missing_needed == 0 and actual == expected)

    print(f"One-to-one check (keep == needed_unique AND no missing): {one_to_one_ok}")
    print(f"  needed_unique = {expected}")
    print(f"  keep_lines    = {actual}")

    # Write outputs
    out_keep = Path("processed_keep_first.txt")
    out_dups = Path("processed_duplicates_to_remove.txt")
    out_unmapped = Path("processed_unmapped.txt")
    out_missing_needed = Path("missing_needed_aod_basenames.txt")
    out_extra = Path("mapped_not_in_needed_aod_basenames.txt")
    out_report = Path("summary_report.txt")

    out_keep.write_text("\n".join(keep_processed_lines) + ("\n" if keep_processed_lines else ""), encoding="utf-8")
    out_dups.write_text("\n".join(dup_processed_lines) + ("\n" if dup_processed_lines else ""), encoding="utf-8")
    out_unmapped.write_text("\n".join(unmapped_processed_lines) + ("\n" if unmapped_processed_lines else ""), encoding="utf-8")
    out_missing_needed.write_text("\n".join(missing_needed) + ("\n" if missing_needed else ""), encoding="utf-8")
    out_extra.write_text("\n".join(extra_mapped_not_needed) + ("\n" if extra_mapped_not_needed else ""), encoding="utf-8")

    # Build an explicit statement
    all_found = (n_missing_needed == 0)
    verdict = "ALL needed files were found (processed)." if all_found else "NOT all needed files were found (missing exist)."

    summary_lines = []
    summary_lines.append("===============================================")
    summary_lines.append("Processed vs Needed check + dedupe report")
    summary_lines.append("===============================================")
    summary_lines.append(f"Needed list:    {needed_txt} (unique basenames: {n_needed_unique}, total lines: {len(needed_lines)})")
    summary_lines.append(f"Processed list: {processed_txt} (total lines: {len(processed_lines)}, parseable: {len(parsed)})")
    summary_lines.append(f"JSON tag/dir:   {mettag} in {json_dir}")
    if missing_jsons:
        summary_lines.append(f"WARNING: Missing JSON files: {len(missing_jsons)} (see summary_report.txt details)")
    summary_lines.append("")
    summary_lines.append("VERDICT:")
    summary_lines.append(f"  {verdict}")
    summary_lines.append("")
    summary_lines.append("COUNTS:")
    summary_lines.append(f"  Unique needed AOD files:           {n_needed_unique}")
    summary_lines.append(f"  Unique mapped AOD files processed: {n_mapped_unique}")
    summary_lines.append(f"  Missing needed AOD files:          {n_missing_needed}")
    summary_lines.append(f"  Keep processed outputs (1st each): {n_keep}")
    summary_lines.append(f"  Duplicate processed outputs:       {n_dups}")
    summary_lines.append(f"  Unmapped processed outputs:        {n_unmapped}")
    summary_lines.append(f"  Mapped-but-not-needed AOD files:   {len(extra_mapped_not_needed)}")
    summary_lines.append("")
    summary_lines.append("OUTPUT FILES WRITTEN:")
    summary_lines.append(f"  {out_keep}  (deduped processed list; keep these)")
    summary_lines.append(f"  {out_dups}  (duplicates that map to already-seen needed AOD; remove these)")
    summary_lines.append(f"  {out_unmapped}  (processed lines that failed parse/mapping)")
    summary_lines.append(f"  {out_missing_needed}  (missing needed AOD basenames)")
    summary_lines.append(f"  {out_extra}  (mapped AOD basenames not in needed list)")
    summary_lines.append("")

    # Add mapping issue detail (trim in stdout but full in report)
    # We'll write full detail into summary_report.txt
    report_lines = list(summary_lines)
    if missing_jsons:
        report_lines.append("Missing JSON mapping files:")
        report_lines.extend([f"  - {x}" for x in missing_jsons])
        report_lines.append("")
    if unparseable_processed:
        report_lines.append(f"Unparseable processed lines: {len(unparseable_processed)} (first 20 shown)")
        for x in unparseable_processed[:20]:
            report_lines.append(f"  {x}")
        report_lines.append("")
    if missing_p_json:
        report_lines.append(f"Processed entries with missing pX JSON: {len(missing_p_json)} (first 20 shown)")
        for orig, p, n in missing_p_json[:20]:
            report_lines.append(f"  p={p} n={n}  line={orig}")
        report_lines.append("")
    if missing_key_in_json:
        report_lines.append(f"Processed entries with missing key N in JSON: {len(missing_key_in_json)} (first 20 shown)")
        for orig, p, n in missing_key_in_json[:20]:
            report_lines.append(f"  p={p} n={n}  line={orig}")
        report_lines.append("")

    out_report.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    # Print concise stdout summary
    print("\n".join(summary_lines))


if __name__ == "__main__":
    main()

