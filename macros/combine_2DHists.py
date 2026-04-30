#!/usr/bin/env python3

from ROOT import *
import argparse


def combine_hists(
    input_file_name,
    input_hist_names,
    output_hist_name,
):

    # ----------------------------------------
    # Open ROOT file
    # ----------------------------------------
    infile = TFile.Open(input_file_name, "UPDATE")

    if not infile or infile.IsZombie():
        raise RuntimeError(f"Could not open {input_file_name}")


    # ----------------------------------------
    # Get first histogram
    # ----------------------------------------
    h0 = infile.Get(input_hist_names[0])

    if not h0:
        raise RuntimeError(f"Could not find {input_hist_names[0]}")

    if not h0.InheritsFrom("TH2"):
        raise RuntimeError(f"{input_hist_names[0]} is not a TH2")


    # ----------------------------------------
    # Clone first histogram
    # ----------------------------------------
    hout = h0.Clone(output_hist_name)

    # Prevent ownership problems
    hout.SetDirectory(0)

    # Ensure error propagation
    hout.Sumw2()


    # ----------------------------------------
    # Add remaining histograms
    # ----------------------------------------
    for hname in input_hist_names[1:]:

        hin = infile.Get(hname)

        if not hin:
            print(f"[WARN] Missing histogram: {hname}")
            continue

        if not hin.InheritsFrom("TH2"):
            print(f"[WARN] {hname} is not a TH2")
            continue

        print(f"[INFO] Adding {hname}")

        hout.Add(hin)


    # ----------------------------------------
    # Write histogram back to file
    # ----------------------------------------
    infile.cd()

    hout.Write(output_hist_name, TObject.kOverwrite)

    infile.Close()

    print(f"[DONE] Wrote '{output_hist_name}' to {input_file_name}")


def parse_args():

    parser = argparse.ArgumentParser(
        description="Combine multiple TH2 histograms into one histogram"
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input ROOT file"
    )

    parser.add_argument(
        "-o",
        "--outputHist",
        required=True,
        help="Name of output histogram"
    )

    parser.add_argument(
        "-H",
        "--hists",
        nargs="+",
        required=True,
        help="Input histogram names"
    )

    return parser.parse_args()


def main():

    # --------------------------------------------------
    # OPTION 1:
    # Hard-coded settings
    # --------------------------------------------------
    #use_manual_inputs = False
    use_manual_inputs = True

    if use_manual_inputs:

        input_file_name = "/uscms/home/jaking/nobackup/el9/CMSSW_13_3_3/src/KUCMSNtupleizer/KUCMSNtupleizer/KUCMSSkimmer/ecal_config/res2dPlotsTFile.root"

        input_hist_names = [
            "ResMap_387742_398903_SRO_Data_Hist_eg_25Cp_xa_pm24b1200_v0422_EB",
            "ResMap_387742_398903_SRO_Data_Hist_eg_25Dp_xa_pm24b1200_v0422_EB",
            "ResMap_387742_398903_SRO_Data_Hist_eg_25Ep_xa_pm24b1200_v0422_EB",
            "ResMap_387742_398903_SRO_Data_Hist_eg_25Fp_xa_pm24b1200_v0422_EB",
            "ResMap_387742_398903_SRO_Data_Hist_eg_25Gp_xa_pm24b1200_v0422_EB",
        ]

        output_hist_name = "ResMap_387742_398903_SRO_Data_Hist_eg_25p_xa_pm24b1200_v0422_EB"

    # --------------------------------------------------
    # OPTION 2:
    # Command-line arguments
    # --------------------------------------------------
    else:

        args = parse_args()

        input_file_name = args.input
        input_hist_names = args.hists
        output_hist_name = args.outputHist


    # --------------------------------------------------
    # Run
    # --------------------------------------------------
    combine_hists(
        input_file_name,
        input_hist_names,
        output_hist_name,
    )


if __name__ == "__main__":
    main()

# python3 combine_2DHists.py -i myFile.root -o combinedHist -H histA histB histC histD
