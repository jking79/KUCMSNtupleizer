#!/usr/bin/env python3
"""Jointly fit an NSC timing-resolution curve to AllGains MC CSV inputs.

The script combines SRO/DRO and EB/EE measurements for each matching MC
dataset and calibration state.  The fit shares one resolution curve across
the four categories.  It fits nonnegative variance contributions, allowing
the stochastic contribution (and therefore S) to be exactly zero.

By default, a first fit identifies a unique category point in a bin whose
pull is more than 2.5 units from the median pull of the other three categories.
That point is excluded from the final fit; bins without one unique outlier
are kept intact.  The cut can be changed or disabled on the command line.

By default, inputs are read from macros/mc_res_info and outputs are written
to macros/mc_res_info/joint_nsc_fit.  Cali and NoCali groups are fit
separately.
"""

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares, nnls


CATEGORIES = (("SRO", "EB"), ("SRO", "EE"), ("DRO", "EB"), ("DRO", "EE"))
ERROR_COLUMNS = {
    "fitted-scatter": "sigma_fitted_scatter_nsc_error",
    "statistical": "sigma_statistical_error",
}


def number(row, name):
    try:
        value = float(row.get(name, ""))
        return value if math.isfinite(value) else math.nan
    except (TypeError, ValueError):
        return math.nan


def read_category_csv(path, error_column, min_x, max_x):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"No rows in {path}")

    metadata_names = ("run_start", "run_end", "pair_type", "dataset_tag",
                      "gain_selection", "calibration_state", "region")
    metadata = {name: rows[0].get(name, "") for name in metadata_names}
    for row in rows:
        if any(row.get(name, "") != metadata[name] for name in metadata_names):
            raise ValueError(f"Mixed histogram metadata in {path}")

    points = []
    rejected = 0
    for row in rows:
        x = number(row, "x_center")
        y = number(row, "sigma_fitted_scatter_nsc")
        error = number(row, error_column)
        valid = number(row, "gaussian_fit_valid")
        if (not math.isfinite(x) or x <= 0 or x < min_x or x > max_x or
                not math.isfinite(y) or y <= 0 or
                not math.isfinite(error) or error <= 0 or
                not math.isfinite(valid) or valid <= 0):
            rejected += 1
            continue
        points.append({
            "x": x,
            "x_low": number(row, "x_low"),
            "x_high": number(row, "x_high"),
            "sigma": y,
            "error": error,
            "occupancy": number(row, "occupancy"),
            "x_bin": row.get("x_bin", ""),
            "source_csv": str(path),
        })
    return metadata, points, rejected


def collect_groups(input_dir, dataset_filter, calibration_state, error_column,
                   min_x, max_x):
    groups = defaultdict(dict)
    rejected_by_file = {}
    candidates = sorted(input_dir.glob("*_resolution_bin_results.csv"))
    for path in candidates:
        with path.open(newline="") as stream:
            first = next(csv.DictReader(stream), None)
        if not first:
            continue
        if first.get("gain_selection") != "AllGains":
            continue
        pair = first.get("pair_type", "")
        region = first.get("region", "")
        if (pair, region) not in CATEGORIES:
            continue
        tag = first.get("dataset_tag", "")
        state = first.get("calibration_state", "")
        if dataset_filter and dataset_filter not in tag:
            continue
        if calibration_state != "both" and state.lower() != calibration_state.lower():
            continue

        metadata, points, rejected = read_category_csv(
            path, error_column, min_x, max_x)
        key = (metadata["run_start"], metadata["run_end"], tag, state)
        category = (pair, region)
        if category in groups[key]:
            raise ValueError(f"Duplicate {category} CSV for {key}: {path}")
        groups[key][category] = points
        rejected_by_file[path.name] = rejected

    if not groups:
        raise ValueError(f"No matching AllGains SRO/DRO EB/EE CSVs found in {input_dir}")

    for key, categories in groups.items():
        missing = set(CATEGORIES) - set(categories)
        if missing:
            missing_text = ", ".join(f"{pair}-{region}" for pair, region in sorted(missing))
            raise ValueError(f"Incomplete four-category group {key}; missing {missing_text}")
        for category in CATEGORIES:
            if not categories[category]:
                raise ValueError(f"No usable points for {key}, {category}")
    return groups, rejected_by_file


def solve_nonnegative(design, y, errors):
    """Weighted least squares in sigma, with nonnegative variance terms."""
    variance_error = np.maximum(2.0 * y * errors, 1e-12)
    weighted_design = design / variance_error[:, None]
    weighted_y = (y * y) / variance_error
    linear_seed, _ = nnls(weighted_design, weighted_y)
    scale = max(float(np.median(y * y)), 1e-8)
    seeds = [
        np.maximum(linear_seed, scale * 1e-5),
        np.array([scale / 3.0] * design.shape[1]),
    ]
    if design.shape[1] == 3:
        seeds.extend((
            np.array([0.85 * scale, 0.10 * scale, 0.05 * scale]),
            np.array([0.05 * scale, 0.85 * scale, 0.10 * scale]),
            np.array([0.05 * scale, 0.10 * scale, 0.85 * scale]),
        ))

    def residuals(terms):
        prediction = np.sqrt(np.maximum(design @ terms, 0.0))
        return (prediction - y) / errors

    results = [least_squares(
        residuals,
        seed,
        bounds=(np.zeros(design.shape[1]), np.full(design.shape[1], np.inf)),
        method="trf",
        x_scale="jac",
        max_nfev=20000,
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    ) for seed in seeds]
    return min(results, key=lambda result: float(np.dot(result.fun, result.fun)))


def fit_joint(points):
    if len(points) <= 3:
        raise ValueError(f"NSC fit needs at least 4 valid points; got {len(points)}")
    x = np.asarray([point["x"] for point in points], dtype=float)
    y = np.asarray([point["sigma"] for point in points], dtype=float)
    errors = np.asarray([point["error"] for point in points], dtype=float)
    pivot = float(np.exp(np.mean(np.log(x))))
    z = pivot / x

    # At x=pivot, these coefficients are the N, S, and C variance
    # contributions.  This scaling keeps their numerical magnitudes similar.
    nsc_design = np.column_stack((z * z, z, np.ones_like(z)))
    nsc = solve_nonnegative(nsc_design, y, errors)
    n_terms = np.maximum(nsc.x, 0.0)
    nsc_chi2 = float(np.dot(nsc.fun, nsc.fun))
    nsc_cov = np.linalg.pinv(nsc.jac.T @ nsc.jac)
    boundary_tolerance = max(float(np.max(n_terms)), 1e-12) * 1e-7
    s_at_zero = bool(n_terms[1] <= boundary_tolerance)

    noise = pivot * math.sqrt(n_terms[0])
    stochastic = math.sqrt(n_terms[1] * pivot)
    constant = math.sqrt(n_terms[2] / 2.0)

    def transformed_error(index, derivative):
        variance = float(nsc_cov[index, index])
        return abs(derivative) * math.sqrt(max(variance, 0.0))

    noise_error = (transformed_error(0, pivot / (2.0 * math.sqrt(n_terms[0])))
                   if n_terms[0] > boundary_tolerance else math.nan)
    stochastic_error = (math.nan if s_at_zero else
                        transformed_error(1, math.sqrt(pivot) / (2.0 * math.sqrt(n_terms[1]))))
    constant_error = (transformed_error(2, 1.0 / (2.0 * math.sqrt(2.0 * n_terms[2])))
                      if n_terms[2] > boundary_tolerance else math.nan)

    nc_design = np.column_stack((z * z, np.ones_like(z)))
    nc = solve_nonnegative(nc_design, y, errors)
    nc_terms = np.maximum(nc.x, 0.0)
    nc_chi2 = float(np.dot(nc.fun, nc.fun))

    predictions = np.sqrt(np.maximum(nsc_design @ n_terms, 0.0))
    points_out = []
    for point, prediction in zip(points, predictions):
        item = dict(point)
        item["fit_sigma"] = float(prediction)
        item["residual"] = point["sigma"] - float(prediction)
        item["pull"] = item["residual"] / point["error"]
        points_out.append(item)

    total_at_pivot = float(np.sum(n_terms))
    fractions = n_terms / total_at_pivot if total_at_pivot > 0 else np.full(3, math.nan)
    return {
        "pivot_x": pivot,
        "n_points": len(points),
        "ndof_nsc": len(points) - 3,
        "n_terms_at_pivot": float(n_terms[0]),
        "s_terms_at_pivot": float(n_terms[1]),
        "c_terms_at_pivot": float(n_terms[2]),
        "n_fraction_at_pivot": float(fractions[0]),
        "s_fraction_at_pivot": float(fractions[1]),
        "c_fraction_at_pivot": float(fractions[2]),
        "N": noise,
        "N_error": noise_error,
        "S": stochastic,
        "S_error": stochastic_error,
        "S_at_zero": s_at_zero,
        "C": constant,
        "C_error": constant_error,
        "chi2_nsc": nsc_chi2,
        "chi2_ndof_nsc": nsc_chi2 / (len(points) - 3),
        "fit_success": bool(nsc.success),
        "fit_message": str(nsc.message),
        "chi2_nc": nc_chi2,
        "ndof_nc": len(points) - 2,
        "delta_chi2_nsc_vs_nc": max(nc_chi2 - nsc_chi2, 0.0),
        "N_nc": pivot * math.sqrt(nc_terms[0]),
        "C_nc": math.sqrt(nc_terms[1] / 2.0),
        "points": points_out,
        "curve_design": nsc_design,
        "curve_terms": n_terms,
    }


def point_key(point):
    return (point["source_csv"], str(point["x_bin"]))


def add_peer_pull_filter(categories, initial_fit, threshold):
    """Exclude a unique pull outlier when the other three categories agree."""
    initial_by_key = {point_key(point): point for point in initial_fit["points"]}
    bins = defaultdict(list)
    annotated = {category: [] for category in CATEGORIES}
    for category in CATEGORIES:
        for source in categories[category]:
            point = dict(source)
            point["pair_type"], point["region"] = category
            initial = initial_by_key[point_key(point)]
            point["initial_fit_sigma"] = initial["fit_sigma"]
            point["initial_pull"] = initial["pull"]
            point["peer_pull_median"] = math.nan
            point["peer_pull_delta"] = math.nan
            point["included_in_fit"] = True
            point["filter_reason"] = ""
            bin_key = (round(point["x_low"], 8), round(point["x_high"], 8))
            bins[bin_key].append(point)
            annotated[category].append(point)

    excluded = []
    for bin_points in bins.values():
        if len(bin_points) != len(CATEGORIES):
            continue
        candidates = []
        for point in bin_points:
            peer_pulls = [other["initial_pull"] for other in bin_points if other is not point]
            peer_center = median(peer_pulls)
            delta = point["initial_pull"] - peer_center
            point["peer_pull_median"] = peer_center
            point["peer_pull_delta"] = delta
            if abs(delta) > threshold:
                candidates.append(point)
        # Keep a bin intact if it does not have one clearly isolated category.
        if len(candidates) == 1:
            point = candidates[0]
            point["included_in_fit"] = False
            point["filter_reason"] = (
                f"initial pull differs from other-category median by "
                f"{point['peer_pull_delta']:.3g} (> {threshold:g})")
            excluded.append(point)
    return annotated, excluded


def model_prediction(x, fit):
    z = fit["pivot_x"] / x
    a, b, c = fit["curve_terms"]
    return math.sqrt(max(a * z * z + b * z + c, 0.0))


def summarize_pulls(fit):
    pulls = np.asarray([point["pull"] for point in fit["points"]], dtype=float)
    center = float(np.mean(pulls))
    rms = float(np.sqrt(np.mean(pulls * pulls)))
    centered_rms = float(np.std(pulls, ddof=1)) if len(pulls) > 1 else 0.0
    skew = (float(np.mean((pulls-center)**3) / centered_rms**3)
            if centered_rms > 0 else math.nan)
    fit["pull_mean"] = center
    fit["pull_median"] = float(np.median(pulls))
    fit["pull_rms"] = rms
    fit["pull_skew"] = skew
    fit["fraction_pulls_within_1"] = float(np.mean(np.abs(pulls) <= 1.0))
    fit["fraction_pulls_within_2"] = float(np.mean(np.abs(pulls) <= 2.0))


def safe_stem(parts):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", "_".join(parts)).strip("_")


def write_outputs(output_dir, group_key, categories, fit, error_name,
                  filter_enabled, filter_threshold):
    run_start, run_end, dataset_tag, state = group_key
    stem = safe_stem((run_start, run_end, dataset_tag, state, "AllGains_joint_NSC"))
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "joint_resolution_fit_summary.csv"
    point_path = output_dir / "joint_resolution_fit_points.csv"
    summary_fields = [
        "run_start", "run_end", "dataset_tag", "calibration_state", "gain_selection",
        "categories", "error_column", "peer_pull_filter_enabled", "peer_pull_cut",
        "n_points_input", "n_points_excluded", "chi2_initial_nsc", "n_points", "ndof_nsc", "N", "N_error",
        "S", "S_error", "S_at_zero", "C", "C_error", "n_terms_at_pivot",
        "s_terms_at_pivot", "c_terms_at_pivot", "n_fraction_at_pivot",
        "s_fraction_at_pivot", "c_fraction_at_pivot", "chi2_nsc", "chi2_ndof_nsc",
        "chi2_nc", "ndof_nc", "delta_chi2_nsc_vs_nc", "N_nc", "C_nc", "fit_success",
        "fit_message", "pull_mean", "pull_median", "pull_rms", "pull_skew",
        "fraction_pulls_within_1", "fraction_pulls_within_2", "plot_file",
    ]
    summary_row = {
        "run_start": run_start, "run_end": run_end, "dataset_tag": dataset_tag,
        "calibration_state": state, "gain_selection": "AllGains",
        "categories": "SRO-EB;SRO-EE;DRO-EB;DRO-EE", "error_column": error_name,
        "peer_pull_filter_enabled": filter_enabled,
        "peer_pull_cut": filter_threshold if filter_enabled else "",
        "plot_file": stem + ".png",
        **{name: fit[name] for name in summary_fields if name in fit},
    }
    write_header = not summary_path.exists()
    with summary_path.open("a", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=summary_fields)
        if write_header:
            writer.writeheader()
        writer.writerow(summary_row)

    point_fields = [
        "run_start", "run_end", "dataset_tag", "calibration_state", "pair_type", "region",
        "gain_selection", "x_bin", "x_low", "x_high", "x", "occupancy", "sigma",
        "error", "initial_pull", "peer_pull_median", "peer_pull_delta",
        "included_in_fit", "filter_reason", "fit_sigma", "residual", "pull", "source_csv",
    ]
    write_header = not point_path.exists()
    with point_path.open("a", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=point_fields)
        if write_header:
            writer.writeheader()
        for category in CATEGORIES:
            for point in categories[category]:
                prediction = model_prediction(point["x"], fit)
                residual = point["sigma"] - prediction
                writer.writerow({
                    "run_start": run_start, "run_end": run_end,
                    "dataset_tag": dataset_tag, "calibration_state": state,
                    "pair_type": category[0], "region": category[1],
                    "gain_selection": "AllGains", "x_bin": point["x_bin"],
                    "x_low": point["x_low"], "x_high": point["x_high"],
                    "x": point["x"], "occupancy": point["occupancy"],
                    "sigma": point["sigma"], "error": point["error"],
                    "initial_pull": point["initial_pull"],
                    "peer_pull_median": point["peer_pull_median"],
                    "peer_pull_delta": point["peer_pull_delta"],
                    "included_in_fit": point["included_in_fit"],
                    "filter_reason": point["filter_reason"],
                    "fit_sigma": prediction, "residual": residual,
                    "pull": residual / point["error"],
                    "source_csv": point["source_csv"],
                })

    plot_group(output_dir / (stem + ".png"), group_key, categories, fit)
    return summary_path, point_path, output_dir / (stem + ".png")


def plot_group(path, group_key, categories, fit):
    colors = {("SRO", "EB"): "tab:blue", ("SRO", "EE"): "tab:orange",
              ("DRO", "EB"): "tab:green", ("DRO", "EE"): "tab:red"}
    markers = {("SRO", "EB"): "o", ("SRO", "EE"): "s",
               ("DRO", "EB"): "^", ("DRO", "EE"): "D"}
    fig, (axis, residual_axis) = plt.subplots(
        2, 1, figsize=(10, 8), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]}, constrained_layout=True)
    all_x = []
    excluded_labeled = False
    for category in CATEGORIES:
        points = categories[category]
        included = [point for point in points if point["included_in_fit"]]
        excluded = [point for point in points if not point["included_in_fit"]]
        x = np.asarray([point["x"] for point in included])
        y = np.asarray([point["sigma"] for point in included])
        error = np.asarray([point["error"] for point in included])
        all_x.extend([point["x"] for point in points])
        label = f"{category[0]} {category[1]}"
        if included:
            axis.errorbar(x, y, yerr=error, fmt=markers[category], color=colors[category],
                          label=label, capsize=2, markersize=5, linestyle="none")
        if excluded:
            axis.errorbar([point["x"] for point in excluded],
                          [point["sigma"] for point in excluded],
                          yerr=[point["error"] for point in excluded], fmt="x",
                          color="black", label=("Excluded peer-pull outlier" if not excluded_labeled else "_nolegend_"),
                          capsize=2, markersize=7, linestyle="none")
            excluded_labeled = True
        for point in points:
            final_pull = (point["sigma"] - model_prediction(point["x"], fit)) / point["error"]
            residual_axis.plot(point["x"], final_pull,
                               marker=(markers[category] if point["included_in_fit"] else "x"),
                               color=(colors[category] if point["included_in_fit"] else "black"),
                               linestyle="none", markersize=4)

    x_curve = np.geomspace(min(all_x), max(all_x), 500)
    z = fit["pivot_x"] / x_curve
    a, b, c = fit["curve_terms"]
    y_curve = np.sqrt(a * z * z + b * z + c)
    axis.plot(x_curve, y_curve, color="black", linewidth=2, label="Shared NSC fit")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_ylabel("Timing resolution σ [ns]")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(ncol=3, fontsize=9)
    run_start, run_end, dataset_tag, state = group_key
    axis.set_title(f"AllGains joint MC resolution: {dataset_tag} ({state}), runs {run_start}–{run_end}\n"
                   f"N={fit['N']:.4g}, S={fit['S']:.4g}, C={fit['C']:.4g}; "
                   f"χ²/ndof={fit['chi2_nsc']:.2f}/{fit['ndof_nsc']}; "
                   f"excluded={fit.get('n_points_excluded', 0)}")
    residual_axis.axhline(0, color="black", linewidth=1)
    residual_axis.axhline(2, color="gray", linewidth=1, linestyle="--")
    residual_axis.axhline(-2, color="gray", linewidth=1, linestyle="--")
    residual_axis.set_xscale("log")
    residual_axis.set_xlabel(r"$A_{eff}/\sigma_n$")
    residual_axis.set_ylabel("Pull")
    residual_axis.grid(True, which="both", alpha=0.25)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    default_input = Path(__file__).resolve().parent / "mc_res_info"
    parser.add_argument("--input-dir", type=Path, default=default_input,
                        help="directory containing the per-histogram resolution-bin CSVs")
    parser.add_argument("--output-dir", type=Path,
                        help="output directory (default: INPUT_DIR/joint_nsc_fit)")
    parser.add_argument("--dataset-tag", default="",
                        help="optional substring filter for dataset_tag")
    parser.add_argument("--calibration-state", choices=("both", "Cali", "NoCali"),
                        default="both", help="fit Cali and NoCali separately")
    parser.add_argument("--error-source", choices=tuple(ERROR_COLUMNS),
                        default="fitted-scatter",
                        help="per-bin error column used for weighting")
    parser.add_argument("--min-x", type=float, default=0.0)
    parser.add_argument("--max-x", type=float, default=math.inf)
    parser.add_argument("--peer-pull-cut", type=float, default=2.5,
                        help="exclude a unique point farther than this many pull units from its bin peers")
    parser.add_argument("--no-peer-pull-filter", action="store_true",
                        help="fit all valid points and only report their peer-pull diagnostics")
    args = parser.parse_args(argv)
    input_dir = args.input_dir.expanduser().resolve()
    output_dir = (args.output_dir or input_dir / "joint_nsc_fit").expanduser().resolve()
    if not input_dir.is_dir():
        parser.error(f"input directory does not exist: {input_dir}")
    if args.min_x < 0 or args.max_x <= args.min_x:
        parser.error("require 0 <= --min-x < --max-x")
    if args.peer_pull_cut <= 0:
        parser.error("--peer-pull-cut must be positive")

    try:
        groups, rejected = collect_groups(
            input_dir, args.dataset_tag, args.calibration_state,
            ERROR_COLUMNS[args.error_source], args.min_x, args.max_x)
        output_dir.mkdir(parents=True, exist_ok=True)
        for generated_csv in ("joint_resolution_fit_summary.csv",
                              "joint_resolution_fit_points.csv"):
            path = output_dir / generated_csv
            if path.exists():
                path.unlink()
        for key, categories in sorted(groups.items()):
            points = []
            for category in CATEGORIES:
                for source in categories[category]:
                    points.append({**source, "pair_type": category[0], "region": category[1]})
            initial_fit = fit_joint(points)
            filter_threshold = math.inf if args.no_peer_pull_filter else args.peer_pull_cut
            annotated_categories, excluded = add_peer_pull_filter(
                categories, initial_fit, filter_threshold)
            accepted_points = [point for category in CATEGORIES
                               for point in annotated_categories[category]
                               if point["included_in_fit"]]
            fit = fit_joint(accepted_points)
            fit["n_points_input"] = len(points)
            fit["n_points_excluded"] = len(excluded)
            fit["chi2_initial_nsc"] = initial_fit["chi2_nsc"]
            fit["peer_pull_cut"] = args.peer_pull_cut
            summarize_pulls(fit)
            summary_path, point_path, plot_path = write_outputs(
                output_dir, key, annotated_categories, fit,
                ERROR_COLUMNS[args.error_source], not args.no_peer_pull_filter,
                args.peer_pull_cut)
            print(f"{key[2]} {key[3]}: kept {fit['n_points']}/{fit['n_points_input']} points; "
                  f"N={fit['N']:.6g}, S={fit['S']:.6g}{' (at zero)' if fit['S_at_zero'] else ''}, "
                  f"C={fit['C']:.6g}, chi2/ndof={fit['chi2_nsc']:.3f}/{fit['ndof_nsc']}; "
                  f"pull mean/median/rms/skew={fit['pull_mean']:.2f}/"
                  f"{fit['pull_median']:.2f}/{fit['pull_rms']:.2f}/{fit['pull_skew']:.2f}")
            for point in excluded:
                print(f"  excluded {point['pair_type']} {point['region']} "
                      f"bin {point['x_low']:g}-{point['x_high']:g}: "
                      f"peer-pull delta={point['peer_pull_delta']:.2f}")
            print(f"  wrote {summary_path}, {point_path}, {plot_path}")
        if rejected:
            total_rejected = sum(rejected.values())
            print(f"Skipped {total_rejected} invalid, empty, or out-of-range bins across "
                  f"{len(rejected)} input CSVs.")
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
