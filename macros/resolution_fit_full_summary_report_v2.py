#!/usr/bin/env python3
"""Versioned full resolution report with CMS-style era presentation plots.

The original resolution_fit_full_summary_report.py remains the report engine.
This wrapper runs it unchanged, then adds ROOT/TDR-style presentation plots.
"""

import argparse
import html
import math
import re
import sys
from pathlib import Path

import resolution_fit_full_summary_report as base

try:
    from ROOT import (TCanvas, TF1, TGraphErrors, TLatex, TLegend, TMultiGraph,
                      gROOT, gStyle, kBlue, kCyan, kGreen, kOrange, kRed,
                      kViolet, kGray)
    from jwk_tdr_style_py import setTDRStyle
    HAS_ROOT = True
    gROOT.SetBatch(True)
except ImportError:
    HAS_ROOT = False


COLORS = (kBlue + 2, kCyan + 3, kGreen + 2, kOrange + 7,
          kViolet + 7, kRed + 1, kGray + 2, kBlue + 4) if HAS_ROOT else ()


def selected_fit(row, bin_rows):
    """Reproduce the report's tail-selection rule for presentation curves."""
    points = []
    for item in bin_rows:
        x = base.number(item, "x_center")
        y = base.number(item, "sigma_fitted_scatter_nsc")
        error = base.number(item, "sigma_fitted_scatter_nsc_error")
        if math.isfinite(x) and x > 0 and math.isfinite(y) and y > 0 and math.isfinite(error) and error > 0:
            points.append((x, y, error))
    inclusive = tuple(base.number(row, name) for name in (
        "fitted_scatter_nsc_noise", "fitted_scatter_nsc_constant",
        "fitted_scatter_nsc_stochastic"))
    if len(points) < 5 or not all(math.isfinite(value) for value in inclusive):
        return inclusive, 0
    minimum_body_bins = 6
    for excluded in range(1, min(3, len(points) - minimum_body_bins) + 1):
        body = points[:-excluded]
        inclusive_chi2 = base._body_chi2(inclusive, body)
        candidate = base._fit_body_points(body, inclusive)
        if candidate is None:
            continue
        parameters, candidate_chi2, _, covariance_status, at_limit = candidate
        improvement = inclusive_chi2 - candidate_chi2
        tail_pulls = base._tail_residuals(parameters, points[-excluded:])
        inclusive_tail_pulls = base._tail_residuals(inclusive, points[-excluded:])
        inconsistent_with_inclusive = all(abs(pull) > 1.0 for pull in inclusive_tail_pulls)
        coherent_tail = (bool(tail_pulls) and
                         (all(pull > 0 for pull in tail_pulls) or all(pull < 0 for pull in tail_pulls)) and
                         sum(pull * pull for pull in tail_pulls) >= 9.0)
        improved_body = (inclusive_chi2 > 0 and improvement >= max(2.0, 0.10 * inclusive_chi2))
        if (improved_body and inconsistent_with_inclusive and coherent_tail and
                covariance_status >= 2 and not at_limit):
            return parameters, excluded
    return inclusive, 0


def era_label(row):
    tag = row.get("dataset_tag", "")
    match = re.search(r"eg_(16|17|18|22|23|24|25|26)", tag)
    year = match.group(1) if match else ""
    if year == "18":
        if "ABC" in tag:
            year += " ABC"
        elif "_D" in tag:
            year += " D"
    return f"20{year}"


def era_key(row):
    return (row.get("run_start", ""), row.get("run_end", ""), row.get("dataset_tag", ""))


def run_label(row):
    try:
        return "Run 2 (13 TeV)" if int(row.get("run_start", 0)) < 350000 else "Run 3 (13.6 TeV)"
    except ValueError:
        return "ECAL timing resolution"


def load_records(summary_csv, output_dir, tail_info=None):
    records = {}
    for row in base.read_csv(summary_csv):
        bin_csv = base.resolve_artifact(row.get("bin_csv_file", ""), summary_csv.parent, output_dir)
        bins = base.read_csv(bin_csv) if bin_csv.is_file() else []
        saved = (tail_info or {}).get(row.get("hist_name", ""), {})
        parameters = saved.get("parameters")
        excluded = saved.get("excluded", 0)
        if not parameters:
            parameters, excluded = selected_fit(row, bins)
        records[row.get("hist_name", "")] = {
            "row": row, "bins": bins, "parameters": parameters, "excluded": excluded,
        }
    return records


def points_for(record, zee_cut=False):
    points = []
    for item in record["bins"]:
        if zee_cut and base.number(item, "occupancy") < 5000:
            continue
        x = base.number(item, "x_center")
        low = base.number(item, "x_low")
        high = base.number(item, "x_high")
        y = base.number(item, "sigma_fitted_scatter_nsc")
        error = base.number(item, "sigma_fitted_scatter_nsc_error")
        if all(math.isfinite(value) for value in (x, low, high, y, error)) and x > 0 and y > 0 and error > 0:
            points.append((x, (high - low) / 2.0, y, error))
    return points


def make_root_plot(curves, output_path, title, subtitle, run_text):
    if not HAS_ROOT or not curves:
        return False
    setTDRStyle()
    gStyle.SetOptStat(0)
    gStyle.SetPadTickX(1); gStyle.SetPadTickY(1)
    gStyle.SetPadTopMargin(0.09)
    canvas = TCanvas("presentation_canvas", "presentation_canvas", 300, 50, 1200, 900)
    canvas.SetLogx(True); canvas.SetLogy(True); canvas.SetGridx(True); canvas.SetGridy(True)
    multigraph = TMultiGraph()
    graphs = []; functions = []; legend = TLegend(0.66, 0.64, 0.93, 0.90)
    legend.SetBorderSize(0); legend.SetTextSize(0.028)
    for index, curve in enumerate(curves):
        points = points_for(curve["record"], zee_cut=curve.get("zee", False))
        if not points:
            continue
        excluded_count = min(curve["record"].get("excluded", 0), len(points) - 1)
        fit_points = points[:-excluded_count] if excluded_count else points
        excluded_points = points[-excluded_count:] if excluded_count else []
        graph = TGraphErrors(len(fit_points))
        for point_index, (x, xerr, y, error) in enumerate(fit_points):
            graph.SetPoint(point_index, x, y)
            graph.SetPointError(point_index, xerr, error)
        color = COLORS[index % len(COLORS)]
        graph.SetMarkerColor(color); graph.SetLineColor(color); graph.SetMarkerStyle(20 + index % 8)
        graph.SetMarkerSize(1.35); graph.SetLineWidth(1)
        multigraph.Add(graph, "P"); graphs.append(graph)
        legend.AddEntry(graph, curve["label"], "ep")
        if excluded_points:
            excluded_graph = TGraphErrors(len(excluded_points))
            for point_index, (x, xerr, y, error) in enumerate(excluded_points):
                excluded_graph.SetPoint(point_index, x, y)
                excluded_graph.SetPointError(point_index, xerr, error)
            excluded_graph.SetMarkerColor(kGray + 2); excluded_graph.SetLineColor(kGray + 2)
            excluded_graph.SetMarkerStyle(24); excluded_graph.SetMarkerSize(1.35)
            multigraph.Add(excluded_graph, "P"); graphs.append(excluded_graph)
        params = curve["record"]["parameters"]
        if len(params) == 3 and all(math.isfinite(value) for value in params):
            function = TF1(f"presentation_fit_{index}",
                           "sqrt((([0]*[0])/(x*x))+(2*[1]*[1])+([2]*[2]/x))", 15.0, 1600.0)
            function.SetParameters(*params); function.SetLineColor(color); function.SetLineWidth(2)
            functions.append(function)
    if not graphs:
        canvas.Close(); return False
    multigraph.Draw("AP")
    multigraph.SetTitle("")
    multigraph.GetXaxis().SetTitle("A_{eff}/#sigma_{n}")
    multigraph.GetYaxis().SetTitle("#sigma(t_{1}-t_{2}) [ns]")
    multigraph.GetXaxis().SetLimits(15.0, 1600.0)
    multigraph.SetMinimum(0.01); multigraph.SetMaximum(10.0)
    multigraph.GetYaxis().SetRangeUser(0.01, 10.0)
    for function in functions:
        function.Draw("same")
    legend.Draw()
    latex = TLatex(); latex.SetNDC(); latex.SetTextFont(42)
    latex.SetTextSize(0.045); latex.DrawLatex(0.16, 0.93, "#bf{CMS} #it{Preliminary}")
    latex.SetTextSize(0.040); latex.DrawLatex(0.70, 0.93, run_text)
    latex.SetTextSize(0.030); latex.DrawLatex(0.16, 0.85, subtitle)
    latex.SetTextSize(0.023)
    latex.DrawLatex(0.175, 0.285, "#sigma_{i}^{2}=(N/(A_{eff}/#sigma_{n}))^{2}+S^{2}/(A_{eff}/#sigma_{n})+2C^{2}")
    latex.SetTextSize(0.022)
    for index, curve in enumerate(curves):
        params = curve["record"]["parameters"]
        if len(params) == 3 and all(math.isfinite(value) for value in params):
            text = (f"{curve['label']}: N={params[0]:.4g}, S={params[2]:.4g}, C={params[1]:.4g}")
            latex.DrawLatex(0.175, 0.24 - index * 0.027, text)
    canvas.Modified(); canvas.Update(); canvas.Print(str(output_path)); canvas.Close()
    return True


def make_presentation_plots(records, output_dir):
    if not HAS_ROOT:
        print("Presentation plots skipped: PyROOT is unavailable.")
        return {}
    presentation_dir = output_dir / "plots" / "presentation"
    presentation_dir.mkdir(parents=True, exist_ok=True)
    eras = {}
    for record in records.values():
        eras.setdefault(era_key(record["row"]), []).append(record)
    generated = {}
    for key, era_records in eras.items():
        by_cat = {(r["row"].get("pair_type"), r["row"].get("gain_selection"),
                   r["row"].get("calibration_state"), r["row"].get("region")): r for r in era_records}
        anchor = era_records[0]["row"]
        stem = "_".join(key)
        era_name = era_label(anchor)
        configs = []
        for region in ("EB", "EE"):
            configs.extend([
                (f"set1_uncal_gain1_{region}", f"{region}: SRO/DRO/ZEE, Gain1Only, NoCali",
                 [(pair, "Gain1Only", "NoCali", region) for pair in ("SRO", "DRO", "ZEE")]),
                (f"set2_gain1_cali_vs_uncali_{region}", f"{region}: SRO/DRO, Gain1Only, Cali vs NoCali",
                 [(pair, "Gain1Only", cali, region) for pair in ("SRO", "DRO") for cali in ("Cali", "NoCali")]),
                (f"set3_allgains_cali_vs_uncali_{region}", f"{region}: SRO/DRO, AllGains, Cali vs NoCali",
                 [(pair, "AllGains", cali, region) for pair in ("SRO", "DRO") for cali in ("Cali", "NoCali")]),
                (f"set4_cali_gain1_vs_allgains_{region}", f"{region}: SRO/DRO, Cali, Gain1Only vs AllGains",
                 [(pair, gain, "Cali", region) for pair in ("SRO", "DRO") for gain in ("Gain1Only", "AllGains")]),
            ])
        configs.append(("set5_cali_allgains_EB_EE", "SRO/DRO, calibrated AllGains, EB and EE",
                        [(pair, "AllGains", "Cali", region) for pair in ("SRO", "DRO") for region in ("EB", "EE")]))
        era_outputs = []
        for name, subtitle, selections in configs:
            curves = []
            for pair, gain, cali, region in selections:
                record = by_cat.get((pair, gain, cali, region))
                if record:
                    label = " ".join((pair, region if name.endswith("EB_EE") else "", gain if "set1" not in name else "", cali if "set1" not in name else "")).strip()
                    curves.append({"record": record, "label": label, "zee": pair == "ZEE"})
            if not curves:
                continue
            path = presentation_dir / f"{stem}_{name}.png"
            if make_root_plot(curves, path, stem, subtitle, run_label(anchor)):
                era_outputs.append((name, subtitle, path.name))
        generated[key] = era_outputs
        page_name = f"presentation_{stem}.html"
        with (output_dir / "html" / page_name).open("w") as output:
            output.write("<!doctype html><html><head><meta charset='utf-8'><title>Era presentation plots</title></head><body>")
            output.write("<p><a href='resolution_fit_summary.html'>Back to summary</a></p>")
            output.write(f"<h1>{html.escape(era_name)}: {html.escape(key[0])}-{html.escape(key[1])}</h1>")
            for _, subtitle, filename in era_outputs:
                output.write(f"<h2>{html.escape(subtitle)}</h2><img src='../plots/presentation/{html.escape(filename)}' width='1000'>")
            output.write("</body></html>\n")
    index = output_dir / "html" / "presentation_index.html"
    with index.open("w") as output:
        output.write("<!doctype html><html><head><meta charset='utf-8'><title>Era presentation plots</title></head><body><h1>Era presentation plots</h1>")
        for key, plots in generated.items():
            page_name = f"presentation_{'_'.join(key)}.html"
            output.write(f"<h2>{html.escape(era_label({'dataset_tag': key[2]}))}: {html.escape(key[0])}-{html.escape(key[1])}</h2><a href='{html.escape(page_name)}'>Open era page ({len(plots)} plots)</a>")
        output.write("</body></html>\n")
    overview = output_dir / "html" / "resolution_fit_summary.html"
    if overview.is_file():
        content = overview.read_text()
        marker = "<h2>Era presentation plots</h2>"
        section = f"{marker}<p><a href='presentation_index.html'>Open era presentation plot index</a></p>"
        if marker not in content:
            overview.write_text(content.replace("</body>", section + "</body>"))
    return generated


def main():
    presentation_only = "--presentation-only" in sys.argv
    tail_info = {}
    original_body_plot = base.make_body_resolution_plot

    def capture_body_plot(row, bin_rows, path):
        result, info = original_body_plot(row, bin_rows, path)
        tail_info[row.get("hist_name", "")] = info
        return result, info

    if not presentation_only:
        base.make_body_resolution_plot = capture_body_plot
        original_argv = sys.argv[:]
        sys.argv = [argument for argument in sys.argv if argument != "--presentation-only"]
        try:
            base.main()
        finally:
            sys.argv = original_argv
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("summary_csv", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--presentation-only", action="store_true")
    args, _ = parser.parse_known_args()
    summary_csv = args.summary_csv.resolve()
    output_dir = (args.output_dir or summary_csv.parent).resolve()
    records = load_records(summary_csv, output_dir, tail_info)
    generated = make_presentation_plots(records, output_dir)
    print(f"Wrote {sum(len(value) for value in generated.values())} era presentation plots")


if __name__ == "__main__":
    main()
