#!/usr/bin/env python3
"""Make a compact plot and Markdown table from one resolution-bin CSV."""

import argparse
import csv
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


def finite(values):
    return [math.isfinite(value) for value in values]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("-o", "--output-prefix", type=Path,
                        help="Output prefix; defaults to the input filename without .csv")
    args = parser.parse_args()

    prefix = args.output_prefix or args.csv_file.with_suffix("")
    with args.csv_file.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise SystemExit(f"No rows found in {args.csv_file}")

    x = [number(row, "x_center") for row in rows]
    xerr = [(number(row, "x_high") - number(row, "x_low")) / 2.0 for row in rows]
    sigma = [number(row, "sigma_statistical") for row in rows]
    sigma_err = [number(row, "sigma_statistical_error") for row in rows]
    model = [number(row, "sigma_model_aware_nsc") for row in rows]
    model_err = [number(row, "sigma_model_aware_nsc_error") for row in rows]
    scatter = [number(row, "sigma_fitted_scatter_nsc") for row in rows]
    scatter_err = [number(row, "sigma_fitted_scatter_nsc_error") for row in rows]
    occupancy = [number(row, "occupancy") for row in rows]
    probability = [number(row, "gaussian_fit_valid") for row in rows]
    fit_prob = [number(row, "sigma_model_aware_nsc") for row in rows]

    title = rows[0].get("hist_name", args.csv_file.stem)
    plot_name = ""
    if HAS_MATPLOTLIB:
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        ax = axes[0, 0]
        for values, errors, label, color in ((sigma, sigma_err, "statistical", "tab:gray"),
                                              (model, model_err, "model-aware NSC", "tab:blue"),
                                              (scatter, scatter_err, "fitted-scatter NSC", "tab:red")):
            mask = [math.isfinite(xx) and math.isfinite(yy) for xx, yy in zip(x, values)]
            ax.errorbar([xx for xx, keep in zip(x, mask) if keep], [yy for yy, keep in zip(values, mask) if keep],
                        yerr=[ee for ee, keep in zip(errors, mask) if keep], xerr=[ee for ee, keep in zip(xerr, mask) if keep],
                        fmt="o-", ms=3, lw=1, capsize=2, label=label, color=color)
        ax.set_xscale("log"); ax.set_xlabel("Energy/amplitude-bin center"); ax.set_ylabel("Time resolution")
        ax.set_title("Resolution versus energy"); ax.grid(True, which="both", alpha=0.25); ax.legend(fontsize=8)
        ax = axes[0, 1]
        mask = [math.isfinite(xx) and math.isfinite(oo) and oo > 0 for xx, oo in zip(x, occupancy)]
        ax.errorbar([xx for xx, keep in zip(x, mask) if keep], [oo for oo, keep in zip(occupancy, mask) if keep],
                    xerr=[ee for ee, keep in zip(xerr, mask) if keep], fmt="o-", ms=3)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("Energy/amplitude-bin center"); ax.set_ylabel("Occupancy")
        ax.set_title("Bin occupancy"); ax.grid(True, which="both", alpha=0.25)
        ax = axes[1, 0]
        valid = [number(row, "gaussian_fit_valid") for row in rows]
        statuses = [number(row, "gaussian_fit_status") for row in rows]
        cov = [number(row, "gaussian_covariance_status") for row in rows]
        ax.step(x, valid, where="mid", label="Gaussian fit valid")
        ax.step(x, [1 if math.isfinite(ss) and ss == 0 else 0 for ss in statuses], where="mid", label="fit status = 0")
        ax.step(x, [1 if math.isfinite(cc) and cc >= 2 else 0 for cc in cov], where="mid", label="covariance >= 2")
        ax.set_xscale("log"); ax.set_ylim(-0.1, 1.1); ax.set_xlabel("Energy/amplitude-bin center"); ax.set_ylabel("Flag")
        ax.set_title("Per-bin Gaussian-fit quality"); ax.grid(True, which="both", alpha=0.25); ax.legend(fontsize=8)
        ax = axes[1, 1]
        model_err_only = [number(row, "sigma_model_error") for row in rows]
        variations = [number(row, "sigma_model_variations") for row in rows]
        ax.plot(x, model_err_only, "o-", ms=3, label="model error"); ax.plot(x, variations, "s--", ms=3, label="model variations")
        ax.set_xscale("log"); ax.set_xlabel("Energy/amplitude-bin center"); ax.set_ylabel("Value")
        ax.set_title("Gaussian model diagnostics"); ax.grid(True, which="both", alpha=0.25); ax.legend(fontsize=8)
        fig.suptitle(title, fontsize=10)
        fig.savefig(prefix.with_suffix(".png"), dpi=160); plt.close(fig)
        plot_name = prefix.with_suffix('.png').name

    markdown = prefix.with_suffix(".md")
    with markdown.open("w") as output:
        output.write(f"# Resolution-bin report\n\n`{title}`\n\n")
        if plot_name:
            output.write(f"![Resolution diagnostics]({plot_name})\n\n")
        else:
            output.write("Plotting skipped: matplotlib is not installed in this Python environment.\n\n")
        output.write("| Bin | Range | Occupancy | Sigma | Stat. error | Fitted-scatter sigma | FS error | Gaussian valid |\n")
        output.write("|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for index, row in enumerate(rows, 1):
            output.write("| {0} | {1}-{2} | {3:.6g} | {4:.6g} | {5:.6g} | {6:.6g} | {7:.6g} | {8} |\n".format(
                index, number(row, "x_low"), number(row, "x_high"), number(row, "occupancy"),
                number(row, "sigma_statistical"), number(row, "sigma_statistical_error"),
                number(row, "sigma_fitted_scatter_nsc"), number(row, "sigma_fitted_scatter_nsc_error"),
                int(number(row, "gaussian_fit_valid")) if math.isfinite(number(row, "gaussian_fit_valid")) else ""))
    if plot_name:
        print(f"Wrote {prefix.with_suffix('.png')}")
    print(f"Wrote {markdown}")


if __name__ == "__main__":
    main()
