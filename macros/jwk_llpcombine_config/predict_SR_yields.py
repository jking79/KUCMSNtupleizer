#!/usr/bin/env python3
"""
Compute predicted SR yields from BFI output JSON using the shape transfer formula:

    predicted_SR_b11 = (buoy_anchor_yield / CR_anchor_yield) * CR_b11_yield

Where:
    CR_anchor_yield   = data yield in anchor channel b00 bin  (Ch1CRHad100)
    buoy_anchor_yield = data yield in SR channel's anchor bin (Ch#CRXxx00)
    CR_b11_yield      = data yield in anchor channel b11 bin  (Ch1CRHad111)

Usage:
    python3 predict_SR_yields.py SV_SR_SingleAnchor_v49.json
    python3 predict_SR_yields.py SV_SR_SingleAnchor_v49.json --pred-bin 11 --max-sig 5
    python3 predict_SR_yields.py SV_SR_SingleAnchor_v49.json --sig-proc gogoGZ_2500_2400_2300_10
"""

import json
import math
import sys
import argparse

SIG_PREFIXES = ["gogoGZ", "gogo", "stop", "squark"]


def is_signal(proc):
    return any(proc.startswith(p) for p in SIG_PREFIXES)


def total_data(j, binname):
    """Return data yield for a bin (weighted events)."""
    if binname not in j or "data" not in j[binname]:
        return 0.0
    return j[binname]["data"][1]


def all_sig_ranked(j, binname):
    """Return list of (proc_name, weighted_yield) for all signal processes,
    sorted descending by yield, skipping zeros."""
    if binname not in j:
        return []
    sigs = [
        (proc, vals[1])
        for proc, vals in j[binname].items()
        if is_signal(proc) and vals[1] > 0
    ]
    return sorted(sigs, key=lambda x: x[1], reverse=True)


def total_sig(j, binname, sig_proc=None):
    """Return signal yield — specific process if given, else sum all."""
    if binname not in j:
        return 0.0
    if sig_proc:
        return j[binname].get(sig_proc, [0, 0, 0])[1]
    return sum(vals[1] for proc, vals in j[binname].items() if is_signal(proc))


def find_sr_channels(j, anchor_bin_suffix, anchor_bin_name):
    """Auto-detect SR channels: bins ending in anchor_bin_suffix that contain 'CR'
    but are not the anchor channel itself. Returns {sr_ch_name: buoy_anchor_bin}."""
    sr_channels = {}
    for binname in j.keys():
        if not binname.endswith(anchor_bin_suffix):
            continue
        if "CR" not in binname or binname == anchor_bin_name:
            continue
        ch_name    = binname[: -len(anchor_bin_suffix)]
        sr_ch_name = ch_name.replace("CR", "SR", 1)
        sr_channels[sr_ch_name] = binname
    return sr_channels


def main():
    parser = argparse.ArgumentParser(
        description="Predict SR yields via shape transfer from BFI JSON"
    )
    parser.add_argument("json_file", help="BFI output JSON")
    parser.add_argument("--anchor-ch",  default="Ch1CR",
                        help="Anchor channel name (default: Ch1CR)")
    parser.add_argument("--anchor-bin", default="00",
                        help="2-char suffix of the anchor bin (default: 00)")
    parser.add_argument("--pred-bin",   default="11",
                        help="2-char suffix of the bin to predict (default: 11)")
    parser.add_argument("--sig-proc",   default=None,
                        help="Show a single signal process instead of sum")
    parser.add_argument("--max-sig",    default=None, type=int,
                        help="Max number of signal processes to display per channel "
                             "(default: all non-zero)")
    args = parser.parse_args()

    with open(args.json_file) as f:
        j = json.load(f)

    anchor_ch       = args.anchor_ch
    anchor_bin_name = anchor_ch + args.anchor_bin   # e.g. Ch1CRHad100
    cr_pred_bin     = anchor_ch + args.pred_bin     # e.g. Ch1CRHad111

    cr_anchor_data = total_data(j, anchor_bin_name)
    cr_pred_data   = total_data(j, cr_pred_bin)

    # List available signal processes
    first_bin  = next(iter(j))
    sig_procs  = sorted(p for p in j[first_bin] if is_signal(p))
    n_sig      = len(sig_procs)
    if n_sig:
        print(f"\nAvailable signal processes ({n_sig}):")
        for sp in sig_procs[:5]:
            print(f"  {sp}")
        if n_sig > 5:
            print(f"  ... ({n_sig - 5} more)")
        if args.sig_proc:
            print(f"\nShowing single signal: {args.sig_proc}")
        else:
            print(f"\nShowing sum of all {n_sig} signal processes")

    # Anchor channel summary
    print(f"\nAnchor channel : {anchor_ch}")
    print(f"  Anchor bin ({args.anchor_bin}) : {anchor_bin_name:<28} data = {cr_anchor_data:.1f}")
    print(f"  Pred   bin ({args.pred_bin}) : {cr_pred_bin:<28} data = {cr_pred_data:.1f}")

    zero_anchor = cr_anchor_data <= 0
    if zero_anchor:
        print("\nWARNING: anchor bin data yield is zero — TF and predicted bkg will be zero.")

    sr_channels = find_sr_channels(j, args.anchor_bin, anchor_bin_name)
    if not sr_channels:
        print("\nNo SR channels found in JSON.")
        sys.exit(1)

    col = max(len(k) for k in sr_channels) + 2
    sep = "-" * (col + 10 + 14 + 14 + 10)

    hdr = (
        f"\n{'SR Channel':<{col}}"
        f"{'TF':>10}"
        f"{'Pred b'+args.pred_bin:>14}"
        f"{'Sig b'+args.pred_bin:>14}"
        f"{'S/sqrtB':>10}"
    )
    print(hdr)

    for sr_ch, buoy_anchor_bin in sorted(sr_channels.items()):
        sr_pred_bin = sr_ch + args.pred_bin

        buoy_data = total_data(j, buoy_anchor_bin)
        tf        = buoy_data / cr_anchor_data if not zero_anchor else 0.0
        pred_bkg  = tf * cr_pred_data

        mc_sig    = total_sig(j, sr_pred_bin, args.sig_proc)
        ranked    = all_sig_ranked(j, sr_pred_bin)
        if args.max_sig is not None:
            ranked = ranked[: args.max_sig]

        # Skip if nothing to show
        if pred_bkg == 0 and mc_sig == 0 and not ranked:
            continue

        ssqrtb = mc_sig / math.sqrt(pred_bkg) if pred_bkg > 0 else 0.0

        # Channel row
        print(sep)
        print(
            f"{sr_ch:<{col}}"
            f"{tf:>10.4f}"
            f"{pred_bkg:>14.4f}"
            f"{mc_sig:>14.4f}"
            f"{ssqrtb:>10.4f}"
        )

        # Yield calculation breakdown
        print(
            f"  {'':>{col-2}}"
            f"  {buoy_anchor_bin} ({buoy_data:.1f}) / "
            f"{anchor_bin_name} ({cr_anchor_data:.1f}) "
            f"* {cr_pred_bin} ({cr_pred_data:.1f})"
        )

        # Ranked signal sub-rows
        if ranked:
            sig_hdr_indent = " " * 4
            print(f"{sig_hdr_indent}{'Signal process':<{col + 6}}{'Yield':>14}{'S/sqrtB':>10}")
            for rank, (proc, yld) in enumerate(ranked, start=1):
                ssqrtb_proc = yld / math.sqrt(pred_bkg) if pred_bkg > 0 else 0.0
                print(f"  [{rank:>2}] {proc:<{col + 4}}{yld:>14.4f}{ssqrtb_proc:>10.4f}")

    print(sep)


if __name__ == "__main__":
    main()
