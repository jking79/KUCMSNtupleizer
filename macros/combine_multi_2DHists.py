#!/usr/bin/env python3

from ROOT import *
import argparse
import os
import glob


def combine_hists_from_files(
    input_dir,
    input_file_glob,
    input_tags,
    output_tag,
    output_file_name,
):

    # ----------------------------------------
    # Find input ROOT files
    # ----------------------------------------
    input_files = sorted(glob.glob(os.path.join(input_dir, input_file_glob)))

    if not input_files:
        raise RuntimeError(f"No input files found matching: {input_dir}/{input_file_glob}")

    print("[INFO] Input files:")
    for fname in input_files:
        print(f"  {fname}")

    combined_hists = {}

    # ----------------------------------------
    # Loop over files
    # ----------------------------------------
    for input_file_name in input_files:

        infile = TFile.Open(input_file_name, "READ")

        if not infile or infile.IsZombie():
            print(f"[WARN] Could not open {input_file_name}")
            continue

        print(f"[INFO] Reading {input_file_name}")

        # ----------------------------------------
        # Loop over all keys in file
        # ----------------------------------------
        for key in infile.GetListOfKeys():

            obj = key.ReadObj()

            if not obj.InheritsFrom("TH2"):
                continue

            input_hist_name = obj.GetName()
            output_hist_name = input_hist_name

            # Replace sample-specific tag with combined tag
            for tag in input_tags:
                output_hist_name = output_hist_name.replace(tag, output_tag)

            print(f"[INFO]  {input_hist_name}  -->  {output_hist_name}")

            # ----------------------------------------
            # First instance: clone histogram
            # ----------------------------------------
            if output_hist_name not in combined_hists:

                hout = obj.Clone(output_hist_name)
                hout.SetDirectory(0)
                hout.Sumw2()

                combined_hists[output_hist_name] = hout

            # ----------------------------------------
            # Later instances: add histogram
            # ----------------------------------------
            else:

                combined_hists[output_hist_name].Add(obj)

        infile.Close()

    if not combined_hists:
        raise RuntimeError("No TH2 histograms found in input files")

    # ----------------------------------------
    # Write combined histograms
    # ----------------------------------------
    outfile = TFile.Open(output_file_name, "UPDATE")

    if not outfile or outfile.IsZombie():
        raise RuntimeError(f"Could not open output file: {output_file_name}")

    outfile.cd()

    for hname, hist in combined_hists.items():
        print(f"[WRITE] {hname}")
        hist.Write(hname, TObject.kOverwrite)

    outfile.Close()

    print(f"[DONE] Wrote {len(combined_hists)} histograms to {output_file_name}")


def parse_args():

    parser = argparse.ArgumentParser(
        description="Combine matching TH2 histograms from multiple ROOT files"
    )

    parser.add_argument(
        "-d",
        "--inputDir",
        required=True,
        help="Directory containing input ROOT files"
    )

    parser.add_argument(
        "-g",
        "--glob",
        default="*.root",
        help="Input ROOT file glob pattern"
    )

    parser.add_argument(
        "-t",
        "--tags",
        nargs="+",
        required=True,
        help="Input sample tags to replace, e.g. eg_25Cp eg_25Dp eg_25Ep"
    )

    parser.add_argument(
        "-T",
        "--outputTag",
        required=True,
        help="Combined output tag, e.g. eg_25p"
    )

    parser.add_argument(
        "-o",
        "--outputFile",
        required=True,
        help="Output ROOT file"
    )

    return parser.parse_args()


def main():

    # --------------------------------------------------
    # OPTION 1:
    # Hard-coded settings
    # --------------------------------------------------
    use_manual_inputs = True

    if use_manual_inputs:

        #input_dir = "res2droot_eg_25p_cc_g1g2_nocali_xa_pm24b1200_v0422"
        input_dir = "res2droot_eg_25p_cc_g1g2_uncorr_xa_pm24b1200_v0422"

        input_file_glob = "res2dPlots_eg_25*p_xa_pm24b1200_v0422.root"

        input_tags = [
            "eg_25Cp",
            "eg_25Dp",
            "eg_25Ep",
            "eg_25Fp",
            "eg_25Gp",
        ]

        output_tag = "eg_25p"

        output_file_name = "../KUCMSSkimmer/ecal_config/res2dPlotsTFile.root"

    # --------------------------------------------------
    # OPTION 2:
    # Command-line arguments
    # --------------------------------------------------
    else:

        args = parse_args()

        input_dir = args.inputDir
        input_file_glob = args.glob
        input_tags = args.tags
        output_tag = args.outputTag
        output_file_name = args.outputFile

    combine_hists_from_files(
        input_dir=input_dir,
        input_file_glob=input_file_glob,
        input_tags=input_tags,
        output_tag=output_tag,
        output_file_name=output_file_name,
    )


if __name__ == "__main__":
    main()
