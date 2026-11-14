#!/usr/bin/env python3
"""Create an HTML report and comparison plot from resolution_fit_summary.csv."""

import argparse
import csv
import html
import math
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def number(row, name):
    try:
        value = float(row.get(name, ""))
        return value if math.isfinite(value) else math.nan
    except (TypeError, ValueError):
        return math.nan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary_csv", type=Path)
    parser.add_argument("-o", "--output-prefix", type=Path,
                        help="Output prefix; defaults to the input filename without .csv")
    args = parser.parse_args()
    prefix = args.output_prefix or args.summary_csv.with_suffix("")
    with args.summary_csv.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise SystemExit(f"No rows found in {args.summary_csv}")

    plot_names = {}
    if HAS_MATPLOTLIB:
        gain_styles = {"AllGains": ("tab:blue", "o"), "Gain1Only": ("tab:orange", "s")}
        for region in ("EB", "EE"):
            for pair_type in ("SRO", "DRO", "ZEE"):
                for calibration_state in ("Cali", "NoCali"):
                    selected = [row for row in rows
                                if row.get("region") == region
                                and row.get("pair_type") == pair_type
                                and row.get("calibration_state") == calibration_state]
                    if not selected:
                        continue
                    plot_path = prefix.with_name(f"{prefix.name}_{region}_{pair_type}_{calibration_state}.png")
                    fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
                    for gain_selection in ("AllGains", "Gain1Only"):
                        gain_rows = [row for row in selected if row.get("gain_selection") == gain_selection]
                        if not gain_rows:
                            continue
                        run = [number(row, "run_start") for row in gain_rows]
                        color, marker = gain_styles[gain_selection]
                        for axis, value_name, error_name, title in (
                            (axes[0, 0], "fitted_scatter_nsc_noise", "fitted_scatter_nsc_noise_error", "Noise term"),
                            (axes[0, 1], "fitted_scatter_nsc_constant", "fitted_scatter_nsc_constant_error", "Constant term"),
                            (axes[1, 0], "fitted_scatter_nsc_stochastic", "fitted_scatter_nsc_stochastic_error", "Stochastic term"),
                            (axes[1, 1], "fitted_scatter_nsc_fractional_scatter", "fitted_scatter_nsc_fractional_scatter_error", "Fractional scatter (%)"),
                        ):
                            values = [number(row, value_name) for row in gain_rows]
                            errors = [number(row, error_name) for row in gain_rows]
                            if "fractional_scatter" in value_name:
                                values = [100 * value for value in values]
                                errors = [100 * error for error in errors]
                            mask = [math.isfinite(rr) and math.isfinite(vv) for rr, vv in zip(run, values)]
                            axis.errorbar([rr for rr, keep in zip(run, mask) if keep],
                                          [vv for vv, keep in zip(values, mask) if keep],
                                          yerr=[ee for ee, keep in zip(errors, mask) if keep],
                                          fmt=marker + "-", ms=4, capsize=2, alpha=0.8,
                                          label=gain_selection, color=color)
                            axis.set_title(title); axis.set_xlabel("Starting run"); axis.grid(True, alpha=0.25)
                            axis.set_ylabel("Percent" if axis is axes[1, 1] else "Fit value")
                            axis.legend(fontsize=8)
                    fig.suptitle(f"{region} {pair_type} {calibration_state}", fontsize=12)
                    fig.savefig(plot_path, dpi=160); plt.close(fig)
                    plot_names[(region, pair_type, calibration_state)] = plot_path.name

    report = prefix.with_suffix(".html")
    with report.open("w") as output:
        output.write("<!doctype html><html><head><meta charset='utf-8'>")
        output.write("<title>Resolution fit summary</title>")
        output.write("<style>body{font-family:sans-serif} table{border-collapse:collapse;font-size:12px} th,td{border:1px solid #bbb;padding:4px} th{background:#eee;position:sticky;top:0} .failed{background:#ffd6d6} .ok{background:#e5f5e0}</style>")
        output.write("</head><body><h1>Resolution fit summary</h1>")
        output.write(f"<p>Input: <code>{html.escape(str(args.summary_csv))}</code></p>")
        if plot_names:
            output.write("<h2>Subset plots</h2><div style='display:grid;grid-template-columns:repeat(2,1fr);gap:12px'>")
            for key, plot_name in plot_names.items():
                output.write(f"<div><h3>{html.escape(' '.join(key))}</h3><img src='{html.escape(plot_name)}' width='600'></div>")
            output.write("</div>")
        else:
            output.write("<p>Plotting skipped: matplotlib is not installed in this Python environment.</p>")
        success = sum(row.get("processing_status") == "success" for row in rows)
        output.write(f"<p>Histograms: {len(rows)} &nbsp; Successful: {success} &nbsp; Failed: {len(rows)-success}</p>")
        output.write("<table><thead><tr><th>Histogram</th><th>Runs</th><th>Pair</th><th>Gain</th><th>Cali</th><th>Region</th><th>Status</th><th>FS NSC N</th><th>FS NSC S</th><th>FS NSC C</th><th>FS scatter (%)</th><th>Prob.</th><th>Bins</th><th>Subset plot</th><th>Bin CSV</th></tr></thead><tbody>")
        for row in rows:
            status = row.get("processing_status", "")
            css = "ok" if status == "success" else "failed"
            csv_path = Path(row.get("bin_csv_file", "")).name
            output.write(f"<tr class='{css}'>")
            output.write(f"<td><code>{html.escape(row.get('hist_name',''))}</code></td>")
            output.write(f"<td>{html.escape(row.get('run_start',''))}-{html.escape(row.get('run_end',''))}</td>")
            output.write(f"<td>{html.escape(row.get('pair_type',''))}</td><td>{html.escape(row.get('gain_selection',''))}</td>")
            output.write(f"<td>{html.escape(row.get('calibration_state',''))}</td><td>{html.escape(row.get('region',''))}</td><td>{html.escape(status)}</td>")
            output.write(f"<td>{number(row,'fitted_scatter_nsc_noise'):.5g}</td><td>{number(row,'fitted_scatter_nsc_stochastic'):.5g}</td><td>{number(row,'fitted_scatter_nsc_constant'):.5g}</td>")
            output.write(f"<td>{100*number(row,'fitted_scatter_nsc_fractional_scatter'):.4g}</td><td>{number(row,'fitted_scatter_nsc_probability'):.4g}</td><td>{number(row,'fitted_scatter_nsc_fit_points'):.0f}</td>")
            plot_name = plot_names.get((row.get("region"), row.get("pair_type"), row.get("calibration_state")), "")
            if plot_name:
                output.write(f"<td><a href='{html.escape(plot_name)}'>plot</a></td>")
            else:
                output.write("<td></td>")
            if csv_path:
                output.write(f"<td><a href='{html.escape(csv_path)}'>bin CSV</a></td>")
            else:
                output.write("<td></td>")
            output.write("</tr>")
        output.write("</tbody></table></body></html>\n")
    if plot_names:
        print(f"Wrote {len(plot_names)} subset plots with prefix {prefix.name}_<region>_<pair>_<calibration>.png")
    print(f"Wrote {report}")


if __name__ == "__main__":
    main()
