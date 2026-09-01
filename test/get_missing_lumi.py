#!/usr/bin/env python3

import json
import argparse


def ranges_to_set(ranges):
    """
    Convert JSON lumi ranges such as:
        [[1, 10], [15, 20]]
    into a set:
        {1, 2, ..., 10, 15, ..., 20}
    """
    lumis = set()

    for start, end in ranges:
        lumis.update(range(start, end + 1))

    return lumis


def set_to_ranges(values):
    """
    Convert a set of integers into compact inclusive ranges.

    Example:
        {1, 2, 3, 7, 8, 10}

    becomes:
        [[1, 3], [7, 8], [10, 10]]
    """
    if not values:
        return []

    values = sorted(values)

    ranges = []
    start = values[0]
    previous = values[0]

    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue

        ranges.append([start, previous])
        start = value
        previous = value

    ranges.append([start, previous])

    return ranges


def find_missing(subset_data, reference_data):
    """
    Return runs/lumisections present in reference_data but missing
    from subset_data.
    """
    missing = {}

    for run, reference_ranges in reference_data.items():

        reference_lumis = ranges_to_set(reference_ranges)

        if run in subset_data:
            subset_lumis = ranges_to_set(subset_data[run])
        else:
            subset_lumis = set()

        missing_lumis = reference_lumis - subset_lumis

        if missing_lumis:
            missing[run] = set_to_ranges(missing_lumis)

    return missing


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare two CMS-style JSON run/lumisection files and "
            "write lumisections present in the reference file but "
            "missing from the subset file."
        )
    )

    parser.add_argument(
        "subset",
        help="JSON file containing the subset of runs/lumisections"
    )

    parser.add_argument(
        "reference",
        help="JSON file containing the complete/reference runs/lumisections"
    )

    parser.add_argument(
        "output",
        help="Output JSON file containing missing runs/lumisections"
    )

    args = parser.parse_args()

    with open(args.subset, "r") as f:
        subset_data = json.load(f)

    with open(args.reference, "r") as f:
        reference_data = json.load(f)

    missing = find_missing(subset_data, reference_data)

    with open(args.output, "w") as f:
        json.dump(missing, f, indent=2, sort_keys=True)

    print(f"Wrote missing lumisections to: {args.output}")

    total_missing = sum(
        len(ranges_to_set(ranges))
        for ranges in missing.values()
    )

    print(f"Runs with missing lumisections: {len(missing)}")
    print(f"Total missing lumisections:     {total_missing}")


if __name__ == "__main__":
    main()
