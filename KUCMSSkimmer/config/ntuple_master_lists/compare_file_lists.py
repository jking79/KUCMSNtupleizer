#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from collections import defaultdict

KEY_RE = re.compile(r"(mGl-\d+_mN2-\d+_mN1-\d+_ct(?:0p)?\d+)")
ROOT_RE = re.compile(r"kucmsntuple_(\d+)\.root$")
PATH_KEY_RE = re.compile(r"(mGl-\d+_mN2-\d+_mN1-\d+_ct(?:0p)?\d+)")

def normalize_key(key):
    # filenames use ct0p1 / ct0p5, paths may use ct1 / ct5
    key = key.replace("_ct0p1", "_ct1")
    key = key.replace("_ct0p5", "_ct5")
    return key


def get_key_from_name(name):
    m = KEY_RE.search(name)
    if not m:
        return None
    return normalize_key(m.group(1))

def read_root_indices(list_file, expected_key=None):
    indices = set()
    bad_lines = []
    wrong_key_lines = []

    with open(list_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            km = PATH_KEY_RE.search(line)
            if expected_key is not None:
                if not km:
                    wrong_key_lines.append(("NO_KEY", line))
                else:
                    line_key = normalize_key(km.group(1))
                    if line_key != expected_key:
                        wrong_key_lines.append((line_key, line))

            m = ROOT_RE.search(line)
            if not m:
                bad_lines.append(line)
                continue

            indices.add(int(m.group(1)))

    return indices, bad_lines, wrong_key_lines


def collect_lists(directory, pattern):
    out = {}

    for path in sorted(Path(directory).glob(pattern)):
        key = get_key_from_name(path.name)
        if key is None:
            print(f"[WARN] Could not parse key from filename: {path}")
            continue

        if key in out:
            print(f"[WARN] Duplicate key {key}")
            print(f"       existing: {out[key]}")
            print(f"       new:      {path}")

        out[key] = path

    return out


def main():
    parser = argparse.ArgumentParser(
        description="Compare two sets of SMS file lists by mass/ctau key and kucmsntuple index."
    )
    parser.add_argument("--dir-a", default=".", help="Directory containing first set of list files")
    parser.add_argument("--dir-b", default=".", help="Directory containing second set of list files")
    parser.add_argument("--pattern-a", default="SMS_SVIPM100_v34_SMS-GlGl-GZ_FULLMINI_*.txt")
    #parser.add_argument("--pattern-a", default="SMS_gogoGZ_SVHPM100_FullMini_v34_SMS-GlGl-GZ_FULLMINI_*.txt")
    parser.add_argument("--pattern-b", default="SMS_SVHPM100_v34p1_gogoGZ_FULLMINI_*.txt")
    parser.add_argument("--verbose", action="store_true", help="Print detailed missing indices")

    args = parser.parse_args()

    lists_a = collect_lists(args.dir_a, args.pattern_a)
    lists_b = collect_lists(args.dir_b, args.pattern_b)

    keys_a = set(lists_a)
    keys_b = set(lists_b)

    only_a = sorted(keys_a - keys_b)
    only_b = sorted(keys_b - keys_a)
    common = sorted(keys_a & keys_b)

    print("========================================")
    print("List-file pairing summary")
    print("========================================")
    print(f"Set A files: {len(lists_a)}")
    print(f"Set B files: {len(lists_b)}")
    print(f"Matched keys: {len(common)}")
    print(f"Only in A: {len(only_a)}")
    print(f"Only in B: {len(only_b)}")

    if only_a:
        print("\nOnly in A:")
        for k in only_a:
            print(f"  {k}  {lists_a[k]}")

    if only_b:
        print("\nOnly in B:")
        for k in only_b:
            print(f"  {k}  {lists_b[k]}")

    print("\n========================================")
    print("Per-key ROOT index comparison")
    print("========================================")

    n_bad = 0

    for key in common:

        idx_a, bad_a, wrong_key_a = read_root_indices(lists_a[key], key)
        idx_b, bad_b, wrong_key_b = read_root_indices(lists_b[key], key)

        missing_in_b = sorted(idx_a - idx_b)
        missing_in_a = sorted(idx_b - idx_a)

        ok = (
            len(idx_a) == len(idx_b)
            and not missing_in_a
            and not missing_in_b
            and not bad_a
            and not bad_b
            and not wrong_key_a
            and not wrong_key_b
        )

        status = "OK" if ok else "BAD"
        print(f"{status:3s} {key:40s}  A={len(idx_a):4d}  B={len(idx_b):4d}")

        if not ok:
            n_bad += 1

            print(f"    A file: {lists_a[key]}")
            print(f"    B file: {lists_b[key]}")

            if bad_a:
                print(f"    bad/unparsed lines in A: {len(bad_a)}")
            if bad_b:
                print(f"    bad/unparsed lines in B: {len(bad_b)}")

            if missing_in_b:
                print(f"    in A not B: {len(missing_in_b)}")
                if args.verbose:
                    print("      " + " ".join(map(str, missing_in_b)))

            if missing_in_a:
                print(f"    in B not A: {len(missing_in_a)}")
                if args.verbose:
                    print("      " + " ".join(map(str, missing_in_a)))

            if wrong_key_a:
                print(f"    wrong mass/ct key lines in A: {len(wrong_key_a)}")
                if args.verbose:
                    for found_key, line in wrong_key_a[:20]:
                        print(f"      found {found_key}")
                        print(f"        {line}")
            
            if wrong_key_b:
                print(f"    wrong mass/ct key lines in B: {len(wrong_key_b)}")
                if args.verbose:
                    for found_key, line in wrong_key_b[:20]:
                        print(f"      found {found_key}")
                        print(f"        {line}")


    print("\n========================================")
    print("Final summary")
    print("========================================")

    if n_bad == 0 and not only_a and not only_b:
        print("All matched list files have identical kucmsntuple_XXX.root index sets.")
    else:
        print(f"Found issues in {n_bad} matched keys.")
        print(f"Unmatched list files: only A={len(only_a)}, only B={len(only_b)}")


if __name__ == "__main__":
    main()
