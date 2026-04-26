#!/usr/bin/env python3
import json
import subprocess
from typing import Any, Dict, List, Tuple
import sys
import csv

writer = csv.writer(sys.stdout, delimiter=",", lineterminator="\n")

def bashout(cmd: str) -> str:
    """Run shell command and return stdout as text."""
    return subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.STDOUT)

def das_json(query: str) -> List[Dict[str, Any]]:
    """Run dasgoclient --json and parse the returned JSON."""
    cmd = f"dasgoclient --json -query='{query}'"
    out = bashout(cmd).strip()
    if not out:
        return []
    return json.loads(out)

def safe_int(x) -> int:
    try:
        if x is None:
            return 0
        s = str(x).strip()
        return int(s) if s else 0
    except Exception:
        return 0

def get_summary(dataset: str) -> Tuple[int, int, int, int]:
    """
    Returns: (size_bytes, nfiles, nlumis, nevents)
    Uses: dasgoclient --json -query='dataset dataset=...'
    Selects the dbs3:filesummaries record (the one that contains size/nevents/nfiles/nlumis).
    """
    data = das_json(f"dataset dataset={dataset}")

    for rec in data:
        services = rec.get("das", {}).get("services", [])
        if any("dbs3:filesummaries" == s for s in services):
            d = (rec.get("dataset") or [{}])[0]
            size_bytes = safe_int(d.get("size", 0))
            nfiles     = safe_int(d.get("nfiles", d.get("num_file", 0)))
            nlumis     = safe_int(d.get("nlumis", d.get("num_lumi", 0)))
            nevents    = safe_int(d.get("nevents", 0))
            return size_bytes, nfiles, nlumis, nevents

    # If we didn't find filesummaries, return zeros (or raise)
    return 0, 0, 0, 0

def get_disk_sites(dataset: str) -> Tuple[bool, str]:
    """
    Returns: (on_disk_anywhere, disk_sites_csv)
    Uses: dasgoclient site dataset=...
    We treat "on disk" as: any DISK site with block_completion == 100%.
    """
    data = das_json(f"site dataset={dataset}")
    disk_sites = []
    for rec in data:
        # site query returns rec["site"] (usually), sometimes nested in "site":[{...}]
        sites = rec.get("site", [])
        for s in sites:
            name = s.get("name") or s.get("site") or ""
            se   = (s.get("se") or "").upper()
            bc   = s.get("block_completion") or s.get("block_completion_percent") or ""
            # bc often like "100.00%" (string)
            bc_str = str(bc)
            is_complete = bc_str.startswith("100")
            is_disk = ("DISK" in se) or ("MSS" not in se and "TAPE" not in se and "DISK" in se) or ( "TAPE" not in se )
            if name and is_disk and is_complete:
                disk_sites.append(name)

    # de-dup and stable sort
    disk_sites = sorted(set(disk_sites))
    return (len(disk_sites) > 0), ",".join(disk_sites)

def bytes_to_gb(size_bytes: int) -> float:
    # Use decimal GB for spreadsheet friendliness; change to / (1024**3) if you prefer GiB
    return size_bytes / 1e12

# -------------------------
# Your dataset lists
# -------------------------
dl2018MetMini = [
    "/MET/Run2018A-15Feb2022_UL2018-v1/MINIAOD",
    "/MET/Run2018B-15Feb2022_UL2018-v1/MINIAOD",
    "/MET/Run2018C-15Feb2022_UL2018-v1/MINIAOD",
    "/MET/Run2018D-15Feb2022_UL2018-v1/MINIAOD",
]

dlR2Met = [

    '/MET/Run2018A-15Feb2022_UL2018-v1/MINIAOD',
    '/MET/Run2018B-15Feb2022_UL2018-v1/MINIAOD',
    '/MET/Run2018C-15Feb2022_UL2018-v1/MINIAOD',
    '/MET/Run2018D-15Feb2022_UL2018-v1/MINIAOD',
    '/MET/Run2017A-09Aug2019_UL2017_rsb-v1/MINIAOD',
    '/MET/Run2017B-09Aug2019_UL2017_rsb-v1/MINIAOD',
    '/MET/Run2017C-09Aug2019_UL2017_rsb-v1/MINIAOD',
    '/MET/Run2017D-09Aug2019_UL2017_rsb-v1/MINIAOD',
    '/MET/Run2017E-09Aug2019_UL2017_rsb-v1/MINIAOD',
    '/MET/Run2017F-09Aug2019_UL2017_rsb-v1/MINIAOD',
    '/MET/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD',
    '/MET/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD',
    '/MET/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD',
    '/MET/Run2016G-21Feb2020_UL2016-v1/MINIAOD',
    '/MET/Run2016H-21Feb2020_UL2016-v2/MINIAOD',

]

dl2024Mini = [

    "/JetMET0/Run2024A-PromptReco-v1/MINIAOD",
    "/JetMET0/Run2024B-PromptReco-v1/MINIAOD",
    "/JetMET0/Run2024C-PromptReco-v1/MINIAOD",
    "/JetMET0/Run2024D-PromptReco-v1/MINIAOD",
    "/JetMET0/Run2024E-PromptReco-v1/MINIAOD",
    "/JetMET0/Run2024E-PromptReco-v2/MINIAOD",
    "/JetMET0/Run2024F-PromptReco-v1/MINIAOD",
    "/JetMET0/Run2024G-PromptReco-v1/MINIAOD",
    "/JetMET0/Run2024H-PromptReco-v1/MINIAOD",
    "/JetMET0/Run2024I-PromptReco-v1/MINIAOD",
    "/JetMET0/Run2024I-PromptReco-v2/MINIAOD",
    "/JetMET0/Run2024J-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2024A-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2024B-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2024C-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2024D-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2024E-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2024E-PromptReco-v2/MINIAOD",
    "/JetMET1/Run2024F-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2024G-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2024H-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2024I-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2024I-PromptReco-v2/MINIAOD",
    "/JetMET1/Run2024J-PromptReco-v1/MINIAOD",

]

dl2023Mini = [

    "/JetMET0/Run2023B-19Dec2023-v1/MINIAOD",
    "/JetMET0/Run2023C-19Dec2023-v1/MINIAOD",
    "/JetMET0/Run2023D-19Dec2023-v1/MINIAOD",
    "/JetMET1/Run2023B-19Dec2023-v1/MINIAOD",
    "/JetMET1/Run2023C-19Dec2023-v1/MINIAOD",
    "/JetMET1/Run2023D-19Dec2023-v1/MINIAOD",

]

dl2022Mini = [

    "/JetMET/Run2022C-19Dec2023-v1/MINIAOD",
    "/JetMET/Run2022D-19Dec2023-v1/MINIAOD",
    "/JetMET/Run2022E-19Dec2023-v1/MINIAOD",
    "/JetMET/Run2022F-19Dec2023-v2/MINIAOD",
    "/JetMET/Run2022G-19Dec2023-v1/MINIAOD",

]

dl2025Mini = [

    "/JetMET0/Run2025B-PromptReco-v1/MINIAOD",
    "/JetMET0/Run2025C-PromptReco-v1/MINIAOD",
    "/JetMET0/Run2025C-PromptReco-v2/MINIAOD",
    "/JetMET0/Run2025D-PromptReco-v1/MINIAOD",
    "/JetMET0/Run2025E-PromptReco-v1/MINIAOD",
    "/JetMET0/Run2025F-PromptReco-v1/MINIAOD",
    "/JetMET0/Run2025F-PromptReco-v2/MINIAOD",
    "/JetMET0/Run2025G-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2025B-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2025C-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2025C-PromptReco-v2/MINIAOD",
    "/JetMET1/Run2025D-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2025E-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2025F-PromptReco-v1/MINIAOD",
    "/JetMET1/Run2025F-PromptReco-v2/MINIAOD",
    "/JetMET1/Run2025G-PromptReco-v1/MINIAOD",

]

dl2023QCDmini = [

    "/QCD-4Jets_HT-1000to1200_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v4/MINIAODSIM",
    "/QCD-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v4/MINIAODSIM",
    "/QCD-4Jets_HT-1200to1500_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v4/MINIAODSIM",
    "/QCD-4Jets_HT-1500to2000_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3/MINIAODSIM",
    "/QCD-4Jets_HT-2000_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v4/MINIAODSIM",
    "/QCD-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3/MINIAODSIM",
    "/QCD-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v4/MINIAODSIM",
    "/QCD-4Jets_HT-40to70_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3/MINIAODSIM",
    "/QCD-4Jets_HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3/MINIAODSIM",
    "/QCD-4Jets_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3/MINIAODSIM",
    "/QCD-4Jets_HT-800to1000_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v4/MINIAODSIM",

]

dl2023ZJetsmini = [ 

    "/Zto2Nu-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v2/MINIAODSIM",
    "/Zto2Nu-4Jets_HT-1500to2500_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v2/MINIAODSIM",
    "/Zto2Nu-4Jets_HT-1500to2500_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v6_ext1-v2/MINIAODSIM",
    "/Zto2Nu-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v2/MINIAODSIM",
    "/Zto2Nu-4Jets_HT-2500_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v3/MINIAODSIM",
    "/Zto2Nu-4Jets_HT-2500_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v6_ext1-v2/MINIAODSIM",
    "/Zto2Nu-4Jets_HT-400to800_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v2/MINIAODSIM",
    "/Zto2Nu-4Jets_HT-400to800_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v6_ext1-v2/MINIAODSIM",
    "/Zto2Nu-4Jets_HT-800to1500_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v2/MINIAODSIM",
    "/Zto2Nu-4Jets_HT-800to1500_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v6_ext1-v2/MINIAODSIM",

]

datalist = dl2023ZJetsmini    # <-- set the list you actually want

# -------------------------
# TSV output (paste into Sheets/Excel)
# -------------------------
writer.writerow([
    "dataset", "size_TB", "nfiles", "nlumis", "nevents", "on_disk", "disk_sites"
])

for ds in datalist:
    size_bytes, nfiles, nlumis, nevents = get_summary(ds)
    on_disk, disk_sites = get_disk_sites(ds)

    writer.writerow([
        ds.strip(),
        f"{bytes_to_gb(size_bytes):.3f}",
        int(nfiles),
        int(nlumis),
        int(nevents),
        "YES" if on_disk else "NO",
        disk_sites.strip()
    ])

#    print("\t".join([
#        ds,
#        f"{bytes_to_gb(size_bytes):.3f}",
#        str(nfiles),
#        str(nlumis),
#        str(nevents),
#        "YES" if on_disk else "NO",
#        disk_sites
#    ]))

