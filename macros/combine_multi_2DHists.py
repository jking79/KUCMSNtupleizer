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
        
            # Default scale, in case no tag matches
            hist_scale = 1.0
            matched_tag = None
        
            # Replace sample-specific tag with combined tag
            # and get the scale associated with that tag
            for tag, scale in input_tags:
                if tag in input_hist_name:
                    output_hist_name = output_hist_name.replace(tag, output_tag)
                    hist_scale = scale
                    matched_tag = tag
                    break
        
            if matched_tag is None:
                print(f"[WARN] No input tag matched histogram: {input_hist_name}")
                continue
        
            print(
                f"[INFO]  {input_hist_name}  -->  {output_hist_name}"
                f"   scale = {hist_scale}"
            )
        
            # ----------------------------------------
            # First instance: clone and scale histogram
            # ----------------------------------------
            if output_hist_name not in combined_hists:
        
                hout = obj.Clone(output_hist_name)
                hout.SetDirectory(0)
                hout.Sumw2()
                hout.Scale(hist_scale)
        
                combined_hists[output_hist_name] = hout
        
            # ----------------------------------------
            # Later instances: add scaled histogram
            # ----------------------------------------
            else:
        
                combined_hists[output_hist_name].Add(obj, hist_scale)

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
        #input_dir = "res2droot_eg_25p_cc_g1g2_uncorr_xa_pm24b1200_v0422"

        #input_dir = "res2droot_eg_25p_cc_hg_nocali_xa_pm24b1200_v0501"
        #input_dir = "res2droot_eg_25p_cc_lg_nocali_xa_pm24b1200_v0501"
        #input_dir = "res2droot_eg_25p_cc_hg_xa_pm24b1200_v0501"
        #input_dir = "res2droot_eg_25p_cc_lg_xa_pm24b1200_v0501"
        #input_dir = "res2droot_eg_25p_unc_hg_xa_pm24b1200_v0501"
        input_dir = "res2droot_eg_25p_unc_lg_xa_pm24b1200_v0501"

        #input_file_glob = "res2dPlots_eg_25*p_xa_pm24b1200_v0422.root"

        #input_file_glob = "res2dPlots_hg_nocali_eg_25*p_xa_pm24b1200_v0501.root"
        #input_file_glob = "res2dPlots_lg_nocali_cc_eg_25*p_xa_pm24b1200_v0501.root"
        #input_file_glob = "res2dPlots_hg_eg_25*p_xa_pm24b1200_v0501.root"
        #input_file_glob = "res2dPlots_lg_cc_eg_25*p_xa_pm24b1200_v0501.root"
        #input_file_glob = "res2dPlots_hg_unc_eg_25*p_xa_pm24b1200_v0501.root"
        input_file_glob = "res2dPlots_lg_unc_eg_25*p_xa_pm24b1200_v0501.root"

# For era-split data histograms, use scale = 1.0.
# The per-era TH2s preserve their observed event occupancy, so the correct
# relative era contribution is already present. Non-unit scales are kept for
# future use, e.g. MC normalization or explicitly shape-normalized inputs.


        input_tags = [
            ("eg_25Cp",1.0), #21.63),
            ("eg_25Dp",1.0), #25.52),
            ("eg_25Ep",1.0), #14.15),
            ("eg_25Fp",1.0), #26.89),
            ("eg_25Gp",1.0), #22.40),
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
