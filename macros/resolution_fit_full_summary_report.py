#!/usr/bin/env python3
"""Create a browsable HTML/plot bundle for a resolution-fit batch summary."""

import argparse
import csv
import html
import math
import os
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    from ROOT import TF1, TGraphErrors, gROOT
    HAS_ROOT = True
    gROOT.SetBatch(True)
except ImportError:
    HAS_ROOT = False


def number(row, name):
    try:
        value = float(row.get(name, ""))
        return value if math.isfinite(value) else math.nan
    except (TypeError, ValueError):
        return math.nan


def display(value, digits=6):
    return "" if not math.isfinite(value) else f"{value:.{digits}g}"


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_table(output, headers, rows, table_class=""):
    class_attr = f" class='{table_class}'" if table_class else ""
    output.write(f"<table{class_attr}><thead><tr>")
    for header in headers:
        output.write(f"<th>{html.escape(str(header))}</th>")
    output.write("</tr></thead><tbody>")
    for row in rows:
        output.write("<tr>")
        for value in row:
            output.write(f"<td>{value}</td>")
        output.write("</tr>")
    output.write("</tbody></table>")


def resolve_artifact(value, summary_dir, output_dir):
    """Resolve paths written by the batch job from either working directory."""
    if not value:
        return Path()
    path = Path(value)
    candidates = (path, summary_dir / path, summary_dir.parent / path,
                  summary_dir / path.name, summary_dir.parent / path.name,
                  output_dir / path, output_dir / path.name)
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    return path


def make_bin_plot(rows, path, title):
    if not HAS_MATPLOTLIB:
        return False
    x = [number(row, "x_center") for row in rows]
    xerr = [(number(row, "x_high") - number(row, "x_low"))/2 for row in rows]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    ax = axes[0, 0]
    for values_name, errors_name, label, color in (
        ("sigma_statistical", "sigma_statistical_error", "statistical", "tab:gray"),
        ("sigma_model_aware_nsc", "sigma_model_aware_nsc_error", "model-aware NSC", "tab:blue"),
        ("sigma_fitted_scatter_nsc", "sigma_fitted_scatter_nsc_error", "fitted-scatter NSC", "tab:red"),
    ):
        values = [number(row, values_name) for row in rows]
        errors = [number(row, errors_name) for row in rows]
        mask = [math.isfinite(xx) and math.isfinite(yy) for xx, yy in zip(x, values)]
        ax.errorbar([xx for xx, keep in zip(x, mask) if keep], [yy for yy, keep in zip(values, mask) if keep],
                    xerr=[ee for ee, keep in zip(xerr, mask) if keep],
                    yerr=[ee for ee, keep in zip(errors, mask) if keep], fmt="o-", ms=3, capsize=2,
                    label=label, color=color)
    ax.set_xscale("log"); ax.set_xlabel("Energy/amplitude-bin center"); ax.set_ylabel("Time resolution")
    ax.set_title("Resolution versus energy"); ax.grid(True, which="both", alpha=.25); ax.legend(fontsize=8)
    ax = axes[0, 1]
    occupancy = [number(row, "occupancy") for row in rows]
    mask = [math.isfinite(xx) and math.isfinite(oo) and oo > 0 for xx, oo in zip(x, occupancy)]
    ax.errorbar([xx for xx, keep in zip(x, mask) if keep], [oo for oo, keep in zip(occupancy, mask) if keep],
                xerr=[ee for ee, keep in zip(xerr, mask) if keep], fmt="o-")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("Energy/amplitude-bin center"); ax.set_ylabel("Occupancy")
    ax.set_title("Bin occupancy"); ax.grid(True, which="both", alpha=.25)
    ax = axes[1, 0]
    valid = [number(row, "gaussian_fit_valid") for row in rows]
    status = [number(row, "gaussian_fit_status") for row in rows]
    cov = [number(row, "gaussian_covariance_status") for row in rows]
    ax.step(x, valid, where="mid", label="valid")
    ax.step(x, [1 if math.isfinite(v) and v == 0 else 0 for v in status], where="mid", label="status = 0")
    ax.step(x, [1 if math.isfinite(v) and v >= 2 else 0 for v in cov], where="mid", label="covariance >= 2")
    ax.set_xscale("log"); ax.set_ylim(-.1, 1.1); ax.set_xlabel("Energy/amplitude-bin center"); ax.set_ylabel("Flag")
    ax.set_title("Gaussian-fit quality"); ax.grid(True, which="both", alpha=.25); ax.legend(fontsize=8)
    ax = axes[1, 1]
    ax.plot(x, [number(row, "sigma_model_error") for row in rows], "o-", ms=3, label="model error")
    ax.plot(x, [number(row, "sigma_model_variations") for row in rows], "s--", ms=3, label="model variations")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("Energy/amplitude-bin center"); ax.set_ylabel("Value")
    ax.set_title("Gaussian model diagnostics"); ax.grid(True, which="both", alpha=.25); ax.legend(fontsize=8)
    fig.suptitle(title, fontsize=10)
    fig.savefig(path, dpi=160); plt.close(fig)
    return True


def make_resolution_plot(row, bin_rows, path):
    """Plot fitted-scatter NSC points and the saved NSC resolution model."""
    if not HAS_MATPLOTLIB or not bin_rows:
        return False
    points = []
    for item in bin_rows:
        x = number(item, "x_center")
        y = number(item, "sigma_fitted_scatter_nsc")
        error = number(item, "sigma_fitted_scatter_nsc_error")
        if math.isfinite(x) and x > 0 and math.isfinite(y) and y > 0 and math.isfinite(error) and error >= 0:
            points.append((x, y, error))
    if not points:
        return False

    fig, axis = plt.subplots(figsize=(9, 6), constrained_layout=True)
    x_values, y_values, y_errors = zip(*points)
    axis.errorbar(x_values, y_values, yerr=y_errors, fmt="o", ms=4, capsize=2,
                  color="tab:blue", label="Fitted-scatter NSC bins")
    noise = number(row, "fitted_scatter_nsc_noise")
    constant = number(row, "fitted_scatter_nsc_constant")
    stochastic = number(row, "fitted_scatter_nsc_stochastic")
    if all(math.isfinite(value) for value in (noise, constant, stochastic)):
        low, high = min(x_values), max(x_values)
        if high > low > 0:
            log_low, log_high = math.log10(low), math.log10(high)
            curve_x = [10 ** (log_low + index * (log_high - log_low) / 199) for index in range(200)]
            curve_y = [math.sqrt((noise * noise) / (x * x) + 2 * constant * constant +
                                 (stochastic * stochastic) / x) for x in curve_x]
            axis.plot(curve_x, curve_y, color="tab:red", linewidth=2,
                      label="Fitted-scatter NSC model")
            scatter = 100 * number(row, "fitted_scatter_nsc_fractional_scatter")
            parameter_text = (f"N = {display(noise)} ns·GeV\n"
                              f"C = {display(constant)} ns\n"
                              f"S = {display(stochastic)} ns·√GeV\n"
                              f"f = {display(scatter)}%")
            axis.text(0.03, 0.97, parameter_text, transform=axis.transAxes,
                      va="top", ha="left", fontsize=9,
                      bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "0.7"})
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("E_eff [GeV]")
    axis.set_ylabel("σ(t₁ − t₂) [ns]")
    axis.set_title(f"{row.get('hist_name', 'Resolution')} — fitted-scatter NSC", fontsize=8)
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(fontsize=9)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return True


def _resolution_value(parameters, x):
    noise, constant, stochastic = parameters
    return math.sqrt((noise * noise) / (x * x) + 2 * constant * constant +
                     (stochastic * stochastic) / x)


def _body_chi2(parameters, points):
    total = 0.0
    for x, y, error in points:
        sigma = max(error, 1e-9)
        total += ((y - _resolution_value(parameters, x)) / sigma) ** 2
    return total


def _fit_body_points(points, initial):
    """Refit the NSC model with ROOT's weighted least-squares fitter."""
    if not HAS_ROOT or len(points) < 4:
        return None
    graph = TGraphErrors(len(points))
    for index, (x, y, error) in enumerate(points):
        graph.SetPoint(index, x, y)
        graph.SetPointError(index, 0.0, max(error, 1e-9))
    low = min(point[0] for point in points)
    high = max(point[0] for point in points)
    if not high > low > 0:
        return None
    function = TF1("body_resolution_model", "sqrt((([0]*[0])/(x*x))+(2*[1]*[1])+([2]*[2]/x))", low, high)
    noise, constant, stochastic = initial
    noise = max(noise, 1e-3); constant = max(constant, 1e-4); stochastic = max(stochastic, 1e-4)
    function.SetParameters(noise, constant, stochastic)
    upper_limits = (max(100.0, 4.0 * noise), max(2.0, 4.0 * constant), max(10.0, 4.0 * stochastic))
    for index, upper in enumerate(upper_limits):
        function.SetParLimits(index, 0.0, upper)
    try:
        fit_result = graph.Fit(function, "QSN")
        status = int(fit_result)
        covariance_status = int(fit_result.Get().CovMatrixStatus())
    except Exception:
        return None
    parameters = tuple(float(function.GetParameter(index)) for index in range(3))
    if status != 0 or not all(math.isfinite(value) for value in parameters):
        return None
    at_limit = False
    for value, high_limit in zip(parameters, upper_limits):
        tolerance = 1e-5 * max(1.0, abs(value), abs(high_limit))
        if value <= tolerance or value >= high_limit - tolerance:
            at_limit = True
    return parameters, _body_chi2(parameters, points), status, covariance_status, at_limit


def _tail_residuals(parameters, points):
    return [(y - _resolution_value(parameters, x)) / max(error, 1e-9)
            for x, y, error in points]


def make_body_resolution_plot(row, bin_rows, path):
    """Test excluding up to three tail bins and plot the selected body refit."""
    info = {"tested": 0, "excluded": 0, "accepted": False, "message": "not tested"}
    if not HAS_MATPLOTLIB or not HAS_ROOT:
        info["message"] = "body refit unavailable (requires matplotlib and PyROOT)"
        return False, info
    points = []
    for item in bin_rows:
        x = number(item, "x_center")
        y = number(item, "sigma_fitted_scatter_nsc")
        error = number(item, "sigma_fitted_scatter_nsc_error")
        if math.isfinite(x) and x > 0 and math.isfinite(y) and y > 0 and math.isfinite(error) and error > 0:
            points.append((x, y, error))
    minimum_body_bins = 6
    if len(points) < minimum_body_bins + 1:
        info["message"] = f"fewer than {minimum_body_bins + 1} valid bins; cannot test exclusions"
        return False, info
    inclusive = (number(row, "fitted_scatter_nsc_noise"),
                 number(row, "fitted_scatter_nsc_constant"),
                 number(row, "fitted_scatter_nsc_stochastic"))
    if not all(math.isfinite(value) for value in inclusive):
        info["message"] = "inclusive NSC parameters unavailable"
        return False, info

    selected = None
    # Require a well-supported body refit, a material improvement, and a
    # coherent, >=3-sigma tail departure before excluding any bins.
    diagnostics = []
    for excluded in range(1, min(3, len(points) - minimum_body_bins) + 1):
        body = points[:-excluded]
        inclusive_chi2 = _body_chi2(inclusive, body)
        candidate = _fit_body_points(body, inclusive)
        info["tested"] = excluded
        if candidate is None:
            diagnostics.append(f"{excluded}: fit failed")
            continue
        parameters, candidate_chi2, _, covariance_status, at_limit = candidate
        improvement = inclusive_chi2 - candidate_chi2
        tail_pulls = _tail_residuals(parameters, points[-excluded:])
        inclusive_tail_pulls = _tail_residuals(inclusive, points[-excluded:])
        inconsistent_with_inclusive = all(abs(pull) > 1.0 for pull in inclusive_tail_pulls)
        coherent_tail = (bool(tail_pulls) and
                         (all(pull > 0 for pull in tail_pulls) or all(pull < 0 for pull in tail_pulls)) and
                         sum(pull * pull for pull in tail_pulls) >= 9.0)
        improved_body = (inclusive_chi2 > 0 and improvement >= max(2.0, 0.10 * inclusive_chi2))
        accepted = (improved_body and inconsistent_with_inclusive and coherent_tail and
                    covariance_status >= 2 and not at_limit)
        diagnostics.append(f"{excluded}: Δχ²={improvement:.2f}, inclusive pulls={','.join(f'{p:.2f}' for p in inclusive_tail_pulls)}, body pulls={','.join(f'{p:.2f}' for p in tail_pulls)}, cov={covariance_status}, atLimit={int(at_limit)}")
        if accepted:
            selected = (excluded, body, parameters, inclusive_chi2, candidate_chi2, diagnostics[-1])
            break

    if selected is None:
        excluded, body, parameters = 0, points, inclusive
        inclusive_chi2 = candidate_chi2 = _body_chi2(inclusive, body)
        info["message"] = "inclusive model retained; " + ("; ".join(diagnostics) if diagnostics else "no exclusions testable")
    else:
        excluded, body, parameters, inclusive_chi2, candidate_chi2, decision = selected
        info.update({"excluded": excluded, "accepted": True,
                     "message": f"excluded final {excluded} bin(s); body χ² {inclusive_chi2:.3g} -> {candidate_chi2:.3g}; {decision}"})
    info.update({"parameters": parameters, "body_chi2": candidate_chi2,
                 "inclusive_body_chi2": inclusive_chi2})

    fig, axis = plt.subplots(figsize=(9, 6), constrained_layout=True)
    body_x, body_y, body_error = zip(*body)
    axis.errorbar(body_x, body_y, yerr=body_error, fmt="o", ms=4, capsize=2,
                  color="tab:blue", label="Body bins")
    if excluded:
        tail = points[-excluded:]
        tail_x, tail_y, tail_error = zip(*tail)
        axis.errorbar(tail_x, tail_y, yerr=tail_error, fmt="x", ms=6, capsize=2,
                      color="tab:gray", label=f"Excluded tail ({excluded})")
    low, high = min(point[0] for point in points), max(point[0] for point in points)
    curve_x = [10 ** (math.log10(low) + index * (math.log10(high) - math.log10(low)) / 199) for index in range(200)]
    axis.plot(curve_x, [_resolution_value(inclusive, x) for x in curve_x], "k--", linewidth=1.5,
              label="Inclusive NSC model")
    axis.plot(curve_x, [_resolution_value(parameters, x) for x in curve_x], color="tab:red", linewidth=2,
              label="Body-only NSC refit" if excluded else "Inclusive NSC model")
    axis.set_xscale("log"); axis.set_yscale("log")
    axis.set_xlabel("E_eff [GeV]"); axis.set_ylabel("σ(t₁ − t₂) [ns]")
    axis.set_title(f"{row.get('hist_name', 'Resolution')} — tail-influence test", fontsize=8)
    axis.grid(True, which="both", alpha=0.25); axis.legend(fontsize=9)
    parameter_text = (f"N = {display(parameters[0])} ns·GeV\n"
                      f"C = {display(parameters[1])} ns\n"
                      f"S = {display(parameters[2])} ns·√GeV")
    axis.text(0.03, 0.97, parameter_text, transform=axis.transAxes,
              va="top", ha="left", fontsize=9,
              bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "0.7"})
    axis.text(0.03, 0.03, info["message"], transform=axis.transAxes, fontsize=9,
              bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "0.7"})
    fig.savefig(path, dpi=160); plt.close(fig)
    return True, info


def comparison_key(row):
    return tuple(row.get(name, "") for name in
                 ("run_start", "run_end", "pair_type", "dataset_tag", "calibration_state", "region"))


def matched_comparison_bins(gain1_rows, all_rows):
    def key(item):
        return (round(number(item, "x_low"), 8), round(number(item, "x_high"), 8))
    gain1 = {key(item): item for item in gain1_rows if math.isfinite(number(item, "x_low"))}
    all_gain = {key(item): item for item in all_rows if math.isfinite(number(item, "x_low"))}
    common = []
    for bin_key in sorted(set(gain1) & set(all_gain), key=lambda value: value[0]):
        common.append((gain1[bin_key], all_gain[bin_key]))
    return common, len(gain1), len(all_gain)


def selected_fit_parameters(record):
    info = record.get("body_info", {})
    parameters = info.get("parameters")
    if parameters and len(parameters) == 3:
        return tuple(parameters)
    return tuple(number(record["row"], name) for name in
                 ("fitted_scatter_nsc_noise", "fitted_scatter_nsc_constant", "fitted_scatter_nsc_stochastic"))


def make_comparison_plot(gain1_record, all_record, common, path):
    if not HAS_MATPLOTLIB or not common:
        return False
    points = []
    for g1, ag in common:
        g1_values = (number(g1, "x_center"), number(g1, "sigma_fitted_scatter_nsc"), number(g1, "sigma_fitted_scatter_nsc_error"))
        ag_values = (number(ag, "x_center"), number(ag, "sigma_fitted_scatter_nsc"), number(ag, "sigma_fitted_scatter_nsc_error"))
        if all(math.isfinite(value) and value > 0 for value in (*g1_values, *ag_values)):
            points.append((g1_values, ag_values, number(g1, "occupancy"), number(ag, "occupancy")))
    if not points:
        return False
    x = [point[0][0] for point in points]
    g1_y = [point[0][1] for point in points]; g1_e = [point[0][2] for point in points]
    ag_y = [point[1][1] for point in points]; ag_e = [point[1][2] for point in points]
    g1_params = selected_fit_parameters(gain1_record)
    ag_params = selected_fit_parameters(all_record)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    ax = axes[0, 0]
    ax.errorbar(x, g1_y, yerr=g1_e, fmt="o", ms=4, capsize=2, label="Gain1Only bins", color="tab:orange")
    ax.errorbar(x, ag_y, yerr=ag_e, fmt="s", ms=4, capsize=2, label="AllGains bins", color="tab:blue")
    low, high = min(x), max(x)
    curve_x = [10 ** (math.log10(low) + index * (math.log10(high) - math.log10(low)) / 199) for index in range(200)]
    ax.plot(curve_x, [_resolution_value(g1_params, value) for value in curve_x], color="tab:orange", linestyle="--", label="Gain1Only selected fit")
    ax.plot(curve_x, [_resolution_value(ag_params, value) for value in curve_x], color="tab:blue", linestyle="--", label="AllGains selected fit")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("E_eff [GeV]"); ax.set_ylabel("σ(t₁ − t₂) [ns]")
    ax.set_title("Measured bins and selected fits"); ax.grid(True, which="both", alpha=.25); ax.legend(fontsize=8)
    delta = [ag_y[index] - g1_y[index] for index in range(len(x))]
    delta_error = [math.hypot(ag_e[index], g1_e[index]) for index in range(len(x))]
    ax = axes[0, 1]
    ax.errorbar(x, delta, yerr=delta_error, fmt="o-", ms=4, capsize=2, color="tab:purple")
    ax.axhline(0.0, color="black", linewidth=1)
    ax.set_xscale("log"); ax.set_xlabel("E_eff [GeV]"); ax.set_ylabel("AllGains − Gain1Only [ns]")
    ax.set_title("Measured bin difference"); ax.grid(True, which="both", alpha=.25)
    g1_residual = [g1_y[index] - _resolution_value(g1_params, x[index]) for index in range(len(x))]
    ag_residual = [ag_y[index] - _resolution_value(ag_params, x[index]) for index in range(len(x))]
    g1_pull = [g1_residual[index] / g1_e[index] for index in range(len(x))]
    ag_pull = [ag_residual[index] / ag_e[index] for index in range(len(x))]
    ax = axes[1, 0]
    ax.axhline(0.0, color="black", linewidth=1)
    ax.plot(x, g1_pull, "o-", ms=4, color="tab:orange", label="Gain1Only")
    ax.plot(x, ag_pull, "s-", ms=4, color="tab:blue", label="AllGains")
    ax.set_xscale("log"); ax.set_xlabel("E_eff [GeV]"); ax.set_ylabel("(bin − fit) / bin error")
    ax.set_title("Residuals from selected fits"); ax.grid(True, which="both", alpha=.25); ax.legend(fontsize=8)
    ax = axes[1, 1]
    occupancy_g1 = [point[2] for point in points]
    occupancy_ag = [point[3] for point in points]
    ratio = [a / b if math.isfinite(a) and math.isfinite(b) and b > 0 else math.nan for a, b in zip(occupancy_ag, occupancy_g1)]
    ax.axhline(1.0, color="black", linewidth=1); ax.plot(x, ratio, "o-", color="tab:green")
    ax.set_xscale("log"); ax.set_xlabel("E_eff [GeV]"); ax.set_ylabel("AllGains / Gain1Only occupancy")
    ax.set_title("Occupancy ratio"); ax.grid(True, which="both", alpha=.25)
    row = gain1_record["row"]
    fig.suptitle(f"{row.get('hist_name', '')} — Gain1Only vs AllGains", fontsize=10)
    fig.savefig(path, dpi=160); plt.close(fig)
    return True


def make_comparison_pages(records, html_dir, plot_dir):
    grouped = {}
    for record in records.values():
        grouped.setdefault(comparison_key(record["row"]), {})[record["row"].get("gain_selection", "")] = record
    links = {}
    for key, pair in grouped.items():
        gain1 = pair.get("Gain1Only"); all_gain = pair.get("AllGains")
        if not gain1 or not all_gain:
            continue
        common, n_gain1, n_all = matched_comparison_bins(gain1["bin_rows"], all_gain["bin_rows"])
        if not common:
            continue
        hist = gain1["row"].get("hist_name", "comparison")
        page_name = hist + "_gain_compare.html"
        plot_name = hist + "_gain_compare.png"
        make_comparison_plot(gain1, all_gain, common, plot_dir / plot_name)
        metrics = []
        for index, (g1, ag) in enumerate(common, 1):
            g1_sigma = number(g1, "sigma_fitted_scatter_nsc"); ag_sigma = number(ag, "sigma_fitted_scatter_nsc")
            g1_error = number(g1, "sigma_fitted_scatter_nsc_error"); ag_error = number(ag, "sigma_fitted_scatter_nsc_error")
            difference = ag_sigma - g1_sigma; combined = math.hypot(ag_error, g1_error)
            pull = difference / combined if combined > 0 else math.nan
            relative = 100 * difference / g1_sigma if g1_sigma > 0 else math.nan
            metrics.append({"index": index, "g1": g1, "ag": ag, "difference": difference,
                            "combined": combined, "pull": pull, "relative": relative,
                            "interval_overlap": abs(difference) <= ag_error + g1_error})
        tail = metrics[-min(3, len(metrics)):]
        tail_pulls = [abs(item["pull"]) for item in tail if math.isfinite(item["pull"])]
        tail_signs = [1 if item["difference"] > 0 else -1 for item in tail if item["difference"] != 0]
        sign_consistent = bool(tail_signs) and len(set(tail_signs)) == 1
        growing = len(tail_pulls) >= 2 and tail_pulls[-1] >= tail_pulls[0]
        significant_tail = sum(pull >= 1 for pull in tail_pulls) >= min(2, len(tail_pulls)) and sign_consistent and growing
        first_pull = next((item for item in metrics if math.isfinite(item["pull"]) and abs(item["pull"]) >= 1), None)
        first_nonoverlap = next((item for item in metrics if not item["interval_overlap"]), None)
        largest = max(metrics, key=lambda item: abs(item["difference"]))
        selected1 = gain1["body_info"].get("message", "")
        selected2 = all_gain["body_info"].get("message", "")
        with (html_dir / page_name).open("w") as output:
            output.write("<!doctype html><html><head><meta charset='utf-8'><title>Gain comparison</title>")
            output.write("<style>body{font-family:sans-serif}table{border-collapse:collapse;font-size:12px}th,td{border:1px solid #bbb;padding:4px}th{background:#eee;position:sticky;top:0}.flag{background:#ffd6d6}</style></head><body>")
            output.write("<p><a href='resolution_fit_summary.html'>Back to summary</a></p>")
            output.write(f"<h1>{html.escape(hist)} — Gain1Only vs AllGains</h1>")
            output.write(f"<p>Gain1Only fit: {html.escape(selected1)}<br>AllGains fit: {html.escape(selected2)}</p>")
            output.write(f"<p>Matched bins: {len(common)} (Gain1Only {n_gain1}, AllGains {n_all})</p>")
            output.write(f"<p><b>Tail assessment:</b> {'significant growing difference' if significant_tail else 'no significant growing difference'}; tail uses the final {len(tail)} matched bins.</p>")
            output.write(f"<p>First |pull| ≥ 1: {first_pull['index'] if first_pull else 'none'}; first non-overlapping error bars: {first_nonoverlap['index'] if first_nonoverlap else 'none'}; largest absolute difference: {display(largest['difference'])} ns in bin {largest['index']}.</p>")
            if (plot_dir / plot_name).is_file():
                output.write(f"<img src='../plots/{html.escape(plot_name)}' width='1100'>")
            output.write("<h2>Bin comparison</h2><table><thead><tr><th>Bin</th><th>Range</th><th>Gain1Only σ</th><th>AllGains σ</th><th>Δσ [ns]</th><th>Combined error</th><th>Pull</th><th>Δ [%]</th><th>Interval overlap</th><th>Tail exclusion</th></tr></thead><tbody>")
            gain1_excluded = gain1["body_info"].get("excluded", 0); all_excluded = all_gain["body_info"].get("excluded", 0)
            for item in metrics:
                g1 = item["g1"]; ag = item["ag"]
                tail_label = []
                if gain1_excluded and item["index"] > len(metrics) - gain1_excluded: tail_label.append("Gain1Only")
                if all_excluded and item["index"] > len(metrics) - all_excluded: tail_label.append("AllGains")
                css = " class='flag'" if abs(item["pull"]) >= 1 or not item["interval_overlap"] else ""
                output.write(f"<tr{css}><td>{item['index']}</td><td>{display(number(g1,'x_low'))}-{display(number(g1,'x_high'))}</td><td>{display(number(g1,'sigma_fitted_scatter_nsc'))} ± {display(number(g1,'sigma_fitted_scatter_nsc_error'))}</td><td>{display(number(ag,'sigma_fitted_scatter_nsc'))} ± {display(number(ag,'sigma_fitted_scatter_nsc_error'))}</td><td>{display(item['difference'])}</td><td>{display(item['combined'])}</td><td>{display(item['pull'])}</td><td>{display(item['relative'])}</td><td>{'yes' if item['interval_overlap'] else 'no'}</td><td>{html.escape(', '.join(tail_label))}</td></tr>")
            output.write("</tbody></table></body></html>\n")
        links[gain1["row"].get("hist_name", "")] = page_name
    return links


def make_subset_plots(rows, plot_dir):
    plots = {}
    if not HAS_MATPLOTLIB:
        return plots
    gain_styles = {"AllGains": ("tab:blue", "o"), "Gain1Only": ("tab:orange", "s")}
    fit_fields = (
        ("fitted_scatter_nsc_noise", "fitted_scatter_nsc_noise_error", "Noise term"),
        ("fitted_scatter_nsc_constant", "fitted_scatter_nsc_constant_error", "Constant term"),
        ("fitted_scatter_nsc_stochastic", "fitted_scatter_nsc_stochastic_error", "Stochastic term"),
        ("fitted_scatter_nsc_fractional_scatter", "fitted_scatter_nsc_fractional_scatter_error", "Fractional scatter (%)"),
    )
    for region in ("EB", "EE"):
        for pair in ("SRO", "DRO", "ZEE"):
            for cali in ("Cali", "NoCali"):
                selected = [row for row in rows if row.get("region") == region and row.get("pair_type") == pair and row.get("calibration_state") == cali]
                if not selected:
                    continue
                filename = f"summary_{region}_{pair}_{cali}.png"
                fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
                for gain, (color, marker) in gain_styles.items():
                    gain_rows = [row for row in selected if row.get("gain_selection") == gain]
                    if not gain_rows:
                        continue
                    run = [number(row, "run_start") for row in gain_rows]
                    for axis, value_name, error_name, title in zip(axes.flat, *zip(*fit_fields)):
                        values = [number(row, value_name) for row in gain_rows]
                        errors = [number(row, error_name) for row in gain_rows]
                        if "fractional_scatter" in value_name:
                            values = [100*v for v in values]; errors = [100*v for v in errors]
                        mask = [math.isfinite(rr) and math.isfinite(vv) for rr, vv in zip(run, values)]
                        axis.errorbar([rr for rr, keep in zip(run, mask) if keep], [vv for vv, keep in zip(values, mask) if keep],
                                      yerr=[ee for ee, keep in zip(errors, mask) if keep], fmt=marker+"-", ms=4,
                                      capsize=2, alpha=.8, color=color, label=gain)
                        axis.set_title(title); axis.set_xlabel("Starting run"); axis.set_ylabel("Percent" if "fractional" in value_name else "Fit value")
                        axis.grid(True, alpha=.25); axis.legend(fontsize=8)
                fig.suptitle(f"{region} {pair} {cali}", fontsize=12)
                fig.savefig(plot_dir/filename, dpi=160); plt.close(fig)
                plots[(region, pair, cali)] = filename
    return plots


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary_csv", type=Path)
    parser.add_argument("--output-dir", type=Path,
                        help="Report bundle directory; defaults to the summary CSV directory")
    args = parser.parse_args()
    summary_csv = args.summary_csv.resolve()
    output_dir = (args.output_dir or summary_csv.parent).resolve()
    html_dir = output_dir / "html"
    plot_dir = output_dir / "plots"
    html_dir.mkdir(parents=True, exist_ok=True); plot_dir.mkdir(parents=True, exist_ok=True)
    rows = read_csv(summary_csv)
    if not rows:
        raise SystemExit(f"No rows found in {summary_csv}")
    subset_plots = make_subset_plots(rows, plot_dir)
    page_rows = []
    records = {}
    for row in rows:
        hist = row.get("hist_name", "histogram")
        safe = hist + "_report.html"
        bin_csv = resolve_artifact(row.get("bin_csv_file", ""), summary_csv.parent, output_dir)
        bin_rows = read_csv(bin_csv) if bin_csv.is_file() else []
        bin_plot_name = hist + "_bin_diagnostics.png"
        bin_plot_path = plot_dir / bin_plot_name
        if bin_rows: make_bin_plot(bin_rows, bin_plot_path, hist)
        resolution_plot_name = hist + "_resolution.png"
        resolution_plot_path = plot_dir / resolution_plot_name
        if not make_resolution_plot(row, bin_rows, resolution_plot_path):
            resolution_plot_name = ""
        body_resolution_plot_name = hist + "_body_resolution.png"
        body_resolution_plot_path = plot_dir / body_resolution_plot_name
        body_plot_ok, body_fit_info = make_body_resolution_plot(row, bin_rows, body_resolution_plot_path)
        if not body_plot_ok:
            body_resolution_plot_name = ""
        root_file = resolve_artifact(row.get("output_root_file", ""), summary_csv.parent, output_dir)
        page_rows.append((row, safe, root_file, bin_csv, bin_rows,
                          bin_plot_name if bin_plot_path.exists() else "",
                          resolution_plot_name if resolution_plot_path.exists() else "",
                          body_resolution_plot_name if body_resolution_plot_path.exists() else "",
                          body_fit_info))
        records[hist] = {"row": row, "bin_rows": bin_rows, "body_info": body_fit_info}
    for row, page_name, root_file, bin_csv, bin_rows, bin_plot_name, resolution_plot_name, body_resolution_plot_name, body_fit_info in page_rows:
        hist = row.get("hist_name", "histogram")
        root_href = os.path.relpath(root_file, html_dir) if root_file.is_file() else root_file.name
        csv_href = os.path.relpath(bin_csv, html_dir) if bin_csv.is_file() else bin_csv.name
        with (html_dir/page_name).open("w") as output:
            output.write("<!doctype html><html><head><meta charset='utf-8'><title>Resolution histogram</title>")
            output.write("<style>body{font-family:sans-serif}table{border-collapse:collapse;font-size:12px}th,td{border:1px solid #bbb;padding:4px}th{background:#eee;position:sticky;top:0}.compact td{white-space:nowrap}.ok{background:#e5f5e0}.failed{background:#ffd6d6}</style></head><body>")
            output.write(f"<p><a href='resolution_fit_summary.html'>Back to summary</a></p><h1>{html.escape(hist)}</h1>")
            output.write(f"<p>Run range: {html.escape(row.get('run_start',''))}-{html.escape(row.get('run_end',''))}; {html.escape(row.get('region',''))}; {html.escape(row.get('pair_type',''))}; {html.escape(row.get('gain_selection',''))}; {html.escape(row.get('calibration_state',''))}</p>")
            output.write(f"<p><a href='{html.escape(root_href)}'>ROOT fit file</a> &nbsp; <a href='{html.escape(csv_href)}'>bin CSV</a></p>")
            if resolution_plot_name:
                output.write(f"<h2>Resolution fit</h2><img src='../plots/{html.escape(resolution_plot_name)}' width='900'>")
            if body_resolution_plot_name:
                output.write(f"<h2>Tail-influence test</h2><p>{html.escape(body_fit_info['message'])}</p><img src='../plots/{html.escape(body_resolution_plot_name)}' width='900'>")
            if bin_plot_name: output.write(f"<img src='../plots/{html.escape(bin_plot_name)}' width='1000'>")
            families = ("model_aware_nsc", "model_aware_nc", "stat_only_nsc", "stat_only_nc", "fitted_scatter_nsc", "fitted_scatter_nc")
            summary_rows = []
            for family in families:
                summary_rows.append([html.escape(family), display(number(row, family+"_noise")), display(number(row, family+"_stochastic")), display(number(row, family+"_constant")), display(100*number(row, family+"_fractional_scatter")), display(number(row, family+"_probability")), display(number(row, family+"_valid"), 0)])
            output.write("<h2>Fit results</h2>")
            write_table(output, ["Fit", "Noise", "Stochastic", "Constant", "Scatter (%)", "Probability", "Valid"], summary_rows, "compact")
            if bin_rows:
                compact = []
                for index, item in enumerate(bin_rows, 1):
                    compact.append([str(index), display(number(item,"x_low")), display(number(item,"x_high")), display(number(item,"occupancy")), display(number(item,"sigma_statistical")), display(number(item,"sigma_statistical_error")), display(number(item,"sigma_fitted_scatter_nsc")), display(number(item,"sigma_fitted_scatter_nsc_error")), display(number(item,"gaussian_fit_valid"),0)])
                output.write("<h2>Resolution-bin results</h2>")
                write_table(output, ["Bin","Low","High","Occupancy","Sigma","Stat. error","FS sigma","FS error","Gaussian valid"], compact, "compact")
                output.write("<details><summary>Full CSV columns</summary>")
                write_table(output, list(bin_rows[0]), [[html.escape(str(item.get(key,""))) for key in bin_rows[0]] for item in bin_rows])
                output.write("</details>")
            output.write("</body></html>\n")
    comparison_links = make_comparison_pages(records, html_dir, plot_dir)
    overview = html_dir / "resolution_fit_summary.html"
    with overview.open("w") as output:
        output.write("<!doctype html><html><head><meta charset='utf-8'><title>Full resolution fit summary</title>")
        output.write("<style>body{font-family:sans-serif}table{border-collapse:collapse;font-size:12px}th,td{border:1px solid #bbb;padding:4px}th{background:#eee;position:sticky;top:0}.ok{background:#e5f5e0}.failed{background:#ffd6d6}</style></head><body><h1>Full resolution fit summary</h1>")
        output.write(f"<p>Input: <code>{html.escape(str(summary_csv))}</code></p>")
        if subset_plots:
            output.write("<h2>Subset plots</h2>")
            output.write("<div style='display:grid;grid-template-columns:repeat(2,1fr);gap:16px'>")
            for region in ("EB", "EE"):
                for pair in ("SRO", "DRO", "ZEE"):
                    names = [subset_plots.get((region, pair, cali)) for cali in ("Cali", "NoCali")]
                    if not any(names):
                        continue
                    output.write(f"<div><h3>{region} {pair}</h3><div style='display:flex;gap:8px'>")
                    for cali, name in zip(("Cali", "NoCali"), names):
                        if name:
                            output.write(f"<div style='width:50%'><h4>{cali}</h4><img src='../plots/{html.escape(name)}' width='100%'></div>")
                    output.write("</div></div>")
            output.write("</div>")
        else: output.write("<p>Plotting skipped: matplotlib is not installed.</p>")
        output.write("<h2>Histograms</h2><table><thead><tr><th>Histogram</th><th>Runs</th><th>Pair</th><th>Gain</th><th>Cali</th><th>Region</th><th>Status</th><th>FS N</th><th>FS C</th><th>FS S</th><th>FS scatter (%)</th><th>Tail test</th><th>Gain comparison</th><th>Page</th></tr></thead><tbody>")
        for row, page_name, _, _, _, _, _, _, body_fit_info in page_rows:
            status = row.get("processing_status", ""); css = "ok" if status == "success" else "failed"
            tail_status = (f"excluded {body_fit_info.get('excluded', 0)}"
                           if body_fit_info.get("accepted") else "inclusive retained")
            comparison = comparison_links.get(row.get("hist_name", ""), "")
            comparison_cell = f"<a href='{html.escape(comparison)}'>compare</a>" if comparison else ""
            output.write(f"<tr class='{css}'><td><code>{html.escape(row.get('hist_name',''))}</code></td><td>{html.escape(row.get('run_start',''))}-{html.escape(row.get('run_end',''))}</td><td>{html.escape(row.get('pair_type',''))}</td><td>{html.escape(row.get('gain_selection',''))}</td><td>{html.escape(row.get('calibration_state',''))}</td><td>{html.escape(row.get('region',''))}</td><td>{html.escape(status)}</td><td>{display(number(row,'fitted_scatter_nsc_noise'))}</td><td>{display(number(row,'fitted_scatter_nsc_constant'))}</td><td>{display(number(row,'fitted_scatter_nsc_stochastic'))}</td><td>{display(100*number(row,'fitted_scatter_nsc_fractional_scatter'))}</td><td>{html.escape(tail_status)}</td><td>{comparison_cell}</td><td><a href='{html.escape(page_name)}'>details</a></td></tr>")
        output.write("</tbody></table></body></html>\n")
    print(f"Wrote {overview}")
    print(f"Wrote {len(page_rows)} histogram pages, {len(comparison_links)} gain-comparison pages, and {len(subset_plots)} subset plots")


if __name__ == "__main__":
    main()
