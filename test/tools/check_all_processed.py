#!/usr/bin/env python3
import json
import re
import sys
from pathlib import Path
from collections import Counter, defaultdict
from typing import Optional, Tuple


# --------- helpers ---------
def read_nonblank_lines(path: Path) -> list[str]:
    out = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s:
                out.append(s)
    return out


def basename(p: str) -> str:
    # works for /store/... and relative paths
    return p.rsplit("/", 1)[-1]


# processed line example:
# kucmsntuple_MET_R17_SVIPM100_v31_MET_AOD_Run2017F/p9/260107_163403/0002/kucmsntuple_2604.root
RE_P_AND_N = re.compile(r"/(p[1-9])/.*/kucmsntuple_(\d+)\.root$")


#def parse_processed_line(line: str) -> tuple[str, str] | None:
def parse_processed_line(line: str) -> Optional[Tuple[str, str]]:

    """
    Returns (pX, N_as_string) or None if not parseable.
    """
    m = RE_P_AND_N.search(line)
    if not m:
        return None
    return (m.group(1), m.group(2))


def load_mapping_json(path: Path) -> dict[str, str]:
    """
    Load a mapping json like {"19":[srcpath], "15":[srcpath], ...}
    Return dict: key N (string) -> srcpath (string)
    """
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    out: dict[str, str] = {}
    for k, v in data.items():
        if isinstance(v, list) and len(v) >= 1 and isinstance(v[0], str):
            out[str(k)] = v[0]
        elif isinstance(v, str):
            out[str(k)] = v
        else:
            # unexpected structure; ignore but keep going
            continue
    return out


# --------- main ---------
def main():
    if len(sys.argv) != 5:
        print(
            "Usage:\n"
            f"  {sys.argv[0]} NEEDED.txt PROCESSED.txt JSON_DIR METTAG\n\n"
            "Where:\n"
            "  NEEDED.txt    = lines like /store/data/.../XXXX.root\n"
            "  PROCESSED.txt = lines like .../p9/.../kucmsntuple_19.root\n"
            "  JSON_DIR      = directory containing Met2017F_files_p1.json ... p9.json\n"
            "  METTAG        = json basename prefix, e.g. Met2017F_files  (so files are Met2017F_files_p9.json)\n"
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

    # 1) Load needed AOD basenames
    needed_lines = read_nonblank_lines(needed_txt)
    needed_basenames = [basename(x) for x in needed_lines]
    needed_set = set(needed_basenames)
    needed_counts = Counter(needed_basenames)

    # 2) Load processed lines and extract (pX, N)
    processed_lines = read_nonblank_lines(processed_txt)
    parsed = []
    unparseable_processed = []
    for line in processed_lines:
        pn = parse_processed_line(line)
        if pn is None:
            unparseable_processed.append(line)
        else:
            parsed.append((line, pn[0], pn[1]))  # (original, pX, N)

    # 3) Load json mappings p1..p9
    mappings: dict[str, dict[str, str]] = {}
    missing_jsons = []
    for i in range(1, 10):
        p = f"p{i}"
        jpath = json_dir / f"{mettag}_{p}.json"
        if not jpath.exists():
            missing_jsons.append(str(jpath))
            continue
        mappings[p] = load_mapping_json(jpath)

    # 4) Convert processed outputs -> source AOD basenames (via mapping)
    mapped_sources = []  # list of source basenames derived from processed
    mapping_missing = []  # processed lines where json didn't contain N
    mapping_missing_p = []  # processed lines where pX json missing entirely
    mapped_source_by_processed = {}  # processed line -> source path

    for orig, p, n in parsed:
        if p not in mappings:
            mapping_missing_p.append((orig, p, n))
            continue
        mp = mappings[p]
        if n not in mp:
            mapping_missing.append((orig, p, n))
            continue
        src = mp[n]
        mapped_source_by_processed[orig] = src
        mapped_sources.append(basename(src))

    processed_source_counts = Counter(mapped_sources)
    processed_source_set = set(processed_source_counts.keys())

    # 5) Compare
    missing_needed = sorted(needed_set - processed_source_set)
    extra_processed = sorted(processed_source_set - needed_set)

    # Duplicates: same source processed multiple times
    dup_sources = sorted([k for k, c in processed_source_counts.items() if c > 1])

    # Needed duplicates (rare, but report)
    dup_needed = sorted([k for k, c in needed_counts.items() if c > 1])

    # 6) Report
    print("===============================================")
    print("Check: were all needed AOD files processed?")
    print("===============================================")
    print(f"Needed list:    {needed_txt}  (lines: {len(needed_lines)} unique basenames: {len(needed_set)})")
    print(f"Processed list: {processed_txt}  (lines: {len(processed_lines)} parseable: {len(parsed)})")
    print(f"JSON dir:       {json_dir}")
    print(f"JSON tag:       {mettag}  (expecting {mettag}_p1.json ... {mettag}_p9.json)")
    print("")

    if missing_jsons:
        print("WARNING: Missing JSON mapping files:")
        for x in missing_jsons:
            print(f"  - {x}")
        print("")

    if unparseable_processed:
        print(f"WARNING: {len(unparseable_processed)} processed lines did not match expected pattern '/pX/.../kucmsntuple_N.root'")
        print("  First 10 examples:")
        for x in unparseable_processed[:10]:
            print(f"    {x}")
        print("")

    if mapping_missing_p:
        print(f"WARNING: {len(mapping_missing_p)} processed entries referenced a pX with no JSON loaded.")
        print("  First 10 examples:")
        for orig, p, n in mapping_missing_p[:10]:
            print(f"    p={p} n={n}  line={orig}")
        print("")

    if mapping_missing:
        print(f"WARNING: {len(mapping_missing)} processed entries had no key N in the corresponding JSON mapping.")
        print("  First 10 examples:")
        for orig, p, n in mapping_missing[:10]:
            print(f"    p={p} n={n}  line={orig}")
        print("")

    # Core verdict
    if not missing_needed and not mapping_missing and not mapping_missing_p and not missing_jsons:
        print("✅ RESULT: All needed files appear to have been processed (based on JSON mappings).")
    else:
        print("⚠️ RESULT: Not all needed files are confirmed processed.")
    print("")

    print(f"Missing needed (by basename): {len(missing_needed)}")
    print(f"Extra processed (not in needed list): {len(extra_processed)}")
    print(f"Processed sources derived from mappings: {len(processed_source_set)} unique")
    print("")

    if dup_sources:
        print(f"NOTE: {len(dup_sources)} source files appear processed more than once (same AOD basename mapped multiple times).")
        print("  First 20:")
        for x in dup_sources[:20]:
            print(f"    {x}  (count {processed_source_counts[x]})")
        print("")

    if dup_needed:
        print(f"NOTE: {len(dup_needed)} needed basenames appear more than once in the needed list.")
        print("  First 20:")
        for x in dup_needed[:20]:
            print(f"    {x}  (count {needed_counts[x]})")
        print("")

    # Print missing + extra (trim)
    if missing_needed:
        print("---- Missing needed files (first 50 basenames) ----")
        for x in missing_needed[:50]:
            print(x)
        if len(missing_needed) > 50:
            print(f"... ({len(missing_needed) - 50} more)")
        print("")

    if extra_processed:
        print("---- Extra processed files (first 50 basenames) ----")
        for x in extra_processed[:50]:
            print(x)
        if len(extra_processed) > 50:
            print(f"... ({len(extra_processed) - 50} more)")
        print("")

    # Optional: write detailed outputs alongside current directory
    out_missing = Path("missing_needed.txt")
    out_extra = Path("extra_processed.txt")
    out_mapping_issues = Path("mapping_issues.txt")

    out_missing.write_text("\n".join(missing_needed) + ("\n" if missing_needed else ""), encoding="utf-8")
    out_extra.write_text("\n".join(extra_processed) + ("\n" if extra_processed else ""), encoding="utf-8")

    with out_mapping_issues.open("w", encoding="utf-8") as f:
        if missing_jsons:
            f.write("Missing JSON mapping files:\n")
            for x in missing_jsons:
                f.write(f"{x}\n")
            f.write("\n")
        if unparseable_processed:
            f.write("Unparseable processed lines:\n")
            for x in unparseable_processed:
                f.write(f"{x}\n")
            f.write("\n")
        if mapping_missing_p:
            f.write("Processed entries with missing pX JSON:\n")
            for orig, p, n in mapping_missing_p:
                f.write(f"p={p} n={n}  line={orig}\n")
            f.write("\n")
        if mapping_missing:
            f.write("Processed entries with missing key N in JSON:\n")
            for orig, p, n in mapping_missing:
                f.write(f"p={p} n={n}  line={orig}\n")
            f.write("\n")

    print("Wrote:")
    print(f"  {out_missing}  (missing needed basenames)")
    print(f"  {out_extra}    (extra processed basenames)")
    print(f"  {out_mapping_issues} (parse/mapping problems)")
    print("")


if __name__ == "__main__":
    main()

