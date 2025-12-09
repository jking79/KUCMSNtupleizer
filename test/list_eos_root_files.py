#!/usr/bin/env python3
import subprocess
import os
import argparse

EOS_PREFIX = "root://cmseos.fnal.gov/"
EOS_LS_CMD = "eos root://cmseos.fnal.gov ls {}"

def list_eos_directory(path, recursive=True, extension=".root"):
    """
    List all files with the given extension in an EOS directory recursively.
    Returns a list of full paths.
    """
    files = []

    try:
        result = subprocess.run(EOS_LS_CMD.format(path),
                                shell=True,
                                check=True,
                                stdout=subprocess.PIPE,
                                universal_newlines=True)
        entries = result.stdout.splitlines()
        #print( entries )
    except subprocess.CalledProcessError:
        print(f"⚠️ Could not list directory: {path}")
        return files

    for entry in entries:
        if not entry.endswith(extension):  # directory
            if recursive:
                files.extend(list_eos_directory(path+entry+'/', recursive=True, extension=extension))
        else:  # file
            if entry.endswith(extension):
                files.append(path+entry)

    return files

def save_files_by_index(file_paths, input_dir, ppath, index=1, output_dir="."):
    """
    Split files by a given path index and save each group into a separate txt file.
    """
    os.makedirs(output_dir, exist_ok=True)
    current_group = None
    outf = None

    for path in file_paths:
        cpath = path.replace(input_dir,"")
        parts = cpath.strip("/").split("/")
        key = parts[index] if len(parts) > index else "unknown"

        if key != current_group:
            if outf:
                outf.close()
            current_group = key
            outfile_name = os.path.join(output_dir, f"{key}.txt")
            outf = open(outfile_name, "w")
            print(f"Writing: {outfile_name}")

        opath = ppath + cpath
        outf.write(opath + "\n")

    if outf:
        outf.close()

def main():

    #parser = argparse.ArgumentParser()
    #parser.add_argument("--crabDir",'-d',help='path to crab submission directory',required=True)


    #args = parser.parse_args()
    #if args.crabDir[-1] == "/":
    #    crabdir = args.crabDir[:-1]
    #crabdir = os.path.split(crabdir)[-1]
    #match = "kucmsntuple_"
    #proc = crabdir[crabdir.find(match)+len(match):crabdir.find("_R")]

    #odir = crabdir[:crabdir.rfind(proc)]
    #parts = odir.split("_")
    #odir = "_".join(parts[1:])[:-1] #remove 'crab_' and ending '_'

    # EOS directory to scan
    eos_dir = "/store/user/lpcsusylep/jaking/KUCMSNtuple/"
    #eos_dir = "/store/user/lpcsusylep/jaking/"

    #input_dir = eos_dir + "kucmsntuple_MET_R18_SVIPM100_v31/MET/"
    #ppath = 'MET/'

    #input_dir = eos_dir + "kucmsntuple_EGamma_R18_InvMetPho30_NoSV_v31/EGamma/"
    #ppath = 'EGamma/'

    #input_dir = eos_dir + "kuncali/gammares_cali/DoubleEG/"
    #input_dir = eos_dir + "gammares_ul17/DoubleEG/"
    #ppath = 'DoubleEG/'

    #input_dir = eos_dir + "kucmsntuple_JetHT_R18_InvMET100_nolumimask_v31/JetHT/"
    #ppath = 'JetHT/'

    #input_dir = eos_dir + "kucmsntuple_MET_R16_SVIPM100_v31/MET/"
    #input_dir = eos_dir + "kucmsntuple_MET_R17_SVIPM100_v31/MET/"
    #input_dir = eos_dir + "/"+odir+"/"+proc+"/"
    #ppath = proc+'/'

    subdir = 'kucmsntuple_SMS_Sig_SVIPM100_v31_SMS-GlGlGZ_AODSIM_mGl-2500_mN2-2400_mN1-500_ct0p1/'
    input_dir = eos_dir + "kucmsntuple_SMS_Sig_SVIPM100_v31/CRAB_UserFiles/"+subdir
    ppath = 'kucmsntuple_SMS_Sig_SVIPM100_v31_SMS-GlGlGZ_AODSIM_mGl-2500_mN2-2400_mN1-500_ct0p1'

    print(f"Scanning EOS directory for .root files: {input_dir} ...")
    files = list_eos_directory(input_dir, recursive=True, extension=".root")
    print(f"✅ Found {len(files)} .root files.")

    # Optional: split files by index and save lists
    save_files_by_index(files, input_dir, ppath, index=0, output_dir="./")

if __name__ == "__main__":
    main()

