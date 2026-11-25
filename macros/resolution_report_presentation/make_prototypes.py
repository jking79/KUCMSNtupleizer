#!/usr/bin/env python3
"""Regenerate the three review plots from resolution_report and its selected-fit audit."""
import csv
import glob
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family": "Helvetica",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Helvetica",
    "mathtext.it": "Helvetica:italic",
    "mathtext.bf": "Helvetica:bold",
    "mathtext.sf": "Helvetica",
    "mathtext.default": "regular",
    "mathtext.fallback": "stixsans",
})
from matplotlib.ticker import LogLocator, FuncFormatter, FixedLocator

HERE = Path(__file__).resolve().parent
MACROS = HERE.parent
REPORT = MACROS / "resolution_report"
PLOTS = HERE / "plots"
AUDIT = HERE / "fit_selection_audit.csv"
XMIN, XMAX = 15, 800
YMIN, YMAX = 0.04, 4
COLORS = {"SRO_Cali": "#2025b2", "SRO_NoCali": "#087e7d",
          "DRO_Cali": "#169b2e", "DRO_NoCali": "#f06412",
          "SRO": "#2025b2", "DRO": "#087e7d", "ZEE": "#169b2e",
          "SRO_EB_Cali": "#2025b2", "SRO_EE_Cali": "#087e7d", "SRO_EB_NoCali": "#169b2e", "SRO_EE_NoCali": "#777400",
          "DRO_EB_Cali": "#169b2e", "DRO_EE_Cali": "#f06412",
          "SRO_Gain1Only": "#2025b2", "SRO_AllGains": "#087e7d",
          "DRO_Gain1Only": "#169b2e", "DRO_AllGains": "#f06412"}
MARKERS = ["o", "s", "^", "v", "D"]

with AUDIT.open(newline="") as stream:
    audit = {row["hist_name"]: row for row in csv.DictReader(stream)}


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def find_record(pair, region, cali, gain, year_token):
    choices = []
    for hist, audit_row in audit.items():
        if (audit_row["pair_type"] == pair and audit_row["region"] == region and
            audit_row["calibration_state"] == cali and audit_row["gain_selection"] == gain and
            year_token in audit_row["dataset_tag"]):
            choices.append((hist, audit_row))
    if len(choices) != 1:
        raise RuntimeError(f"Expected one {pair}/{region}/{cali}/{gain}/{year_token} row, found {len(choices)}")
    hist, fitrow = choices[0]
    bin_files = list(REPORT.glob(f"{hist}_NSC_resfit_resolution_bin_results.csv"))
    if not bin_files:
        # Some report versions keep the hist name as the bin CSV basename.
        bin_files = list(REPORT.glob(f"{hist}*resolution_bin_results.csv"))
    if len(bin_files) != 1:
        raise RuntimeError(f"Expected one bin CSV for {hist}, found {len(bin_files)}")
    return fitrow, read_rows(bin_files[0])


def values(record):
    fitrow, bins = record
    pts = []
    for row in bins:
        try:
            x, low, high = (float(row[k]) for k in ("x_center", "x_low", "x_high"))
            y = float(row["sigma_fitted_scatter_nsc"])
            err = float(row["sigma_fitted_scatter_nsc_error"])
        except (ValueError, TypeError):
            continue
        if all(math.isfinite(v) for v in (x, low, high, y, err)) and x > 0 and y > 0 and err > 0:
            pts.append((x, (high-low)/2, y, err))
    excluded = min(int(fitrow["excluded_tail_bins"]), max(0, len(pts)-1))
    return pts, pts[:-excluded] if excluded else pts, pts[-excluded:] if excluded else [], fitrow


def model(x, n, s, c):
    return np.sqrt((n/x)**2 + (s*s)/x + 2*c*c)

# Keep numpy import local to keep the record parsing dependency small.
import numpy as np

def label_formatter(val, _):
    # Label 1–5 subdivisions per decade; leave subdivisions 6–9 blank.
    if val <= 0:
        return ""
    exponent = math.floor(math.log10(val))
    mantissa = val / (10 ** exponent)
    nearest = round(mantissa)
    if abs(mantissa-nearest) > 1e-7 or nearest not in (1, 2, 3, 4, 5):
        return ""
    if abs(val-YMAX) < 1e-9:
        return f"{val:g}"
    if nearest == 1:
        return f"{val:g}"
    return f"{nearest}" if exponent == 0 else (f"{nearest}×10^{exponent}" if exponent > 0 else f"{nearest}×10$^{{{exponent}}}$")


def fmt_digits(value, digits):
    """Format with an exact count of printed digits (decimal point excluded)."""
    value = float(value)
    magnitude = abs(value)
    integer_digits = max(1, len(str(int(magnitude)))) if magnitude >= 1 else 1
    decimals = max(0, digits - integer_digits)
    while True:
        rendered = f"{value:.{decimals}f}"
        digit_count = sum(char.isdigit() for char in rendered)
        if digit_count <= digits or decimals == 0:
            return rendered
        decimals -= 1


def make_plot(filename, curves, subtitle, run_text, xmax_override=None):
    if xmax_override is not None:
        plot_xmax = xmax_override
    elif "ecal_dpg" in filename:
        plot_xmax = 800
    elif any(record[0]["gain_selection"] == "AllGains" for _, _, record in curves):
        plot_xmax = 3200
    else:
        plot_xmax = 1600
    # 4:3 legacy proportions with the plot pad widened 10% in x.
    fig = plt.figure(figsize=(14.667, 10), dpi=288)
    ax = fig.add_axes([0.15, 0.12, 0.80, 0.79])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(XMIN, plot_xmax); ax.set_ylim(YMIN, YMAX)
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1,)))
    if plot_xmax > 800:
        minor_ticks = [float(m * 10**e) for e in range(0, 5) for m in range(2, 10) if XMIN < m * 10**e < plot_xmax]
        minor_ticks.append(float(plot_xmax))
        ax.xaxis.set_minor_locator(FixedLocator(minor_ticks))
    else:
        ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    ax.yaxis.set_major_locator(FixedLocator([0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4]))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    ax.yaxis.set_major_formatter(FuncFormatter(label_formatter))
    ax.yaxis.set_minor_formatter(FuncFormatter(label_formatter))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda val, pos: "" if abs(val-XMIN) < 1e-8 else (f"{int(val)}" if val >= 1 else f"{val:g}")))
    x_labels = (20, 30, 50, 200, 300, 500, 800, 1600, 2000, 3200)
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda val, pos: (f"{int(val)}" if any(abs(val-x) < 1e-7 for x in x_labels) else "")))
    ax.tick_params(axis="both", which="major", labelsize=13.5, direction="in", top=True, right=True, length=8)
    ax.tick_params(axis="both", which="minor", labelsize=13.5, direction="in", top=True, right=True, length=5)
    ax.grid(True, which="both", color="#888888", linestyle=":", linewidth=0.75, alpha=.9)
    ax.set_xlabel(r"$A_{\mathrm{eff}}/\sigma_n$", fontsize=25, labelpad=16)
    ax.set_ylabel(r"$\sigma(t_1-t_2)$ [ns]", fontsize=25, labelpad=18)

    legend_handles, legend_labels = [], []
    for i, (key, label, record) in enumerate(curves):
        points, used, excluded, fitrow = values(record)
        color = COLORS[key]
        x, xerr, y, yerr = (np.array(col) for col in zip(*used))
        err = ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=MARKERS[i], color=color,
                          markersize=9, elinewidth=1.4, capsize=2.7, linewidth=0, zorder=4)
        legend_handles.append(err); legend_labels.append(label)
        if excluded:
            ex, exerr, ey, eyerr = (np.array(col) for col in zip(*excluded))
            ax.errorbar(ex, ey, xerr=exerr, yerr=eyerr, fmt="x", color=color,
                        markersize=9, mew=2, elinewidth=1.4, capsize=2.7, linewidth=0, zorder=5)
        n, c, s = (float(fitrow[f"fit_noise_N"]), float(fitrow["fit_constant_C"]), float(fitrow["fit_stochastic_S"]))
        visible_points = [point for point in points if XMIN <= point[0] <= plot_xmax]
        if visible_points:
            xx = np.geomspace(visible_points[0][0], visible_points[-1][0], 500)
            ax.plot(xx, model(xx, n, s, c), color=color, linewidth=2.3, zorder=2)

    legend = ax.legend(legend_handles, legend_labels, loc="upper right", bbox_to_anchor=(0.985, .985),
                       frameon=False, fontsize=21, handlelength=2.5, labelspacing=.8,
                       prop={"weight": "bold", "size": 21})
    # CMS-style header and in-pad subtitle.
    fig.text(.15, .965, r"$\bf{CMS}$ $\it{Preliminary}$", fontsize=30, ha="left", va="top")
    fig.text(.95, .965, run_text, fontsize=26.25, ha="right", va="top")
    ax.text(.10, .9625, subtitle, transform=ax.transAxes, ha="left", va="top", fontsize=21, linespacing=1.05)

    # Model expression and selected report fit parameters, placed in the plot body.
    parameter_y_top = 0.025 + (len(curves) - 1) * 0.036
    equation_y = parameter_y_top + (0.270 - parameter_y_top) / 3.0 + 0.036
    ax.text(.025, equation_y, r"$\sigma_i^2=\left(\frac{N}{A_{\mathrm{eff}}/\sigma_n}\right)^2 + \frac{S^2}{A_{\mathrm{eff}}/\sigma_n} + 2C^2$",
            transform=ax.transAxes, fontsize=22.75, ha="left", va="bottom")
    for i, (key, label, record) in enumerate(curves):
        _, _, _, fitrow = values(record)
        n, ne = float(fitrow["fit_noise_N"]), float(fitrow["fit_noise_N_error"])
        c, ce = float(fitrow["fit_constant_C"]), float(fitrow["fit_constant_C_error"])
        s, se = float(fitrow["fit_stochastic_S"]), float(fitrow["fit_stochastic_S_error"])
        ax.text(.025, parameter_y_top - i*.036, f"N : {fmt_digits(n, 4)} ± {fmt_digits(ne, 3)} [ns]   "
                f"S : {fmt_digits(s, 3)} ± {fmt_digits(se, 3)} [ns]   "
                f"C : {fmt_digits(c, 5)} ± {fmt_digits(ce, 5)} [ns]",
                transform=ax.transAxes, fontsize=19.25, color=COLORS[key], ha="left", va="bottom")

    if plot_xmax > 800:
        fig.canvas.draw()
        for tick in ax.xaxis.get_minor_ticks():
            if abs(tick.get_loc() - plot_xmax) < 1e-7:
                tick.label1.set_horizontalalignment("right")
    stem = PLOTS / filename
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=288, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)


# Run 2 2018 UL, EB: SRO/DRO, default uncalibrated and Chronus.
run2 = [(f"{p}_{c}", f"{p}" + (" Chronus" if c == "Cali" else ""), find_record(p, "EB", c, "Gain1Only", "18UL_miniv2"))
        for p,c in (("SRO","Cali"),("SRO","NoCali"),("DRO","Cali"),("DRO","NoCali"))]
make_plot("prototype_run2_2018_EB_SRO_DRO", run2, "2018 UL\nEB", "Run 2 (13 TeV)")

# Run 3 2025 Prompt: SRO EB/EE, uncalibrated and Chronus, AllGains.
run3 = [(f"SRO_{region}_{cali}", f"{region}" + (" Chronus" if cali == "Cali" else ""),
         find_record("SRO", region, cali, "AllGains", "25_prompt"))
        for region,cali in (("EB","Cali"),("EE","Cali"),("EB","NoCali"),("EE","NoCali"))]
make_plot("prototype_run3_2025_SRO_EB_EE", run3, "2025 Prompt · Uncorrected CC time\nAllGains · SRO Resolution", "Run 3 (13.6 TeV)")

# ECAL DPG review comparison: EB SRO/DRO/ZEE, uncalibrated Gain1Only.
dpg = [(p, p, find_record(p, "EB", "NoCali", "Gain1Only", "25_prompt")) for p in ("SRO","DRO","ZEE")]
make_plot("prototype_ecal_dpg_2025_EB_SRO_DRO_ZEE", dpg, "2025 Prompt · Uncorrected CC time\nEB", "Run 3 (13.6 TeV)")
