#!/usr/bin/env python3
import os
import sys
import subprocess
import getpass
from datetime import datetime

# === User configuration ===
#/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/KUCMSNtupleizer/KUCMSNtupleizer/KUCMSSkimmer
CMSSW_BASE = f"/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3"
TARBALL = "CMSSW_13_3_3.tar.gz"
SUB_FILE = "run_skim.sub"
SCRIPT = "run_skim.sh"
LOGDIR = "logs"

def make_tarball():
    """Create or update a tarball of the full CMSSW area."""
    base_dir = os.path.dirname(CMSSW_BASE)
    tar_path = os.path.join(base_dir, TARBALL)
    print(f"Creating tarball: {tar_path}")

    if os.path.exists(tar_path):
        os.remove(tar_path)

    subprocess.run([
        "tar", "--exclude=.git", "--exclude=tmp", "-czf",
        tar_path, os.path.basename(CMSSW_BASE)
    ], cwd=base_dir, check=True)
    print("✔️  CMSSW tarball created successfully.\n")

def read_jobs(input_file):
    """Parse input file with 3 arguments per line."""
    jobs = []
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 3:
                raise ValueError(f"Line '{line}' does not have 3 space-delimited entries.")
            jobs.append(parts)
    print(f"Found {len(jobs)} jobs in {input_file}.")
    return jobs

def write_submission(jobs):
    """Write Condor submission file."""
    os.makedirs(LOGDIR, exist_ok=True)
    with open(SUB_FILE, "w") as f:
        f.write(f"""\
universe              = vanilla
executable            = {SCRIPT}
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files  = {TARBALL}
transfer_output_files = *.root

output                = {LOGDIR}/job_$(Cluster)_$(Process).out
error                 = {LOGDIR}/job_$(Cluster)_$(Process).err
log                   = {LOGDIR}/job_$(Cluster)_$(Process).log

+JobFlavour           = "longlunch"

queue ARG1, ARG2, ARG3 from (
""")
        for a1, a2, a3 in jobs:
            f.write(f"  {a1} {a2} {a3}\n")
        f.write(")\n")

def submit():
    """Submit Condor jobs."""
    subprocess.run(["condor_submit", SUB_FILE], check=True)
    print("✔️  Condor submission complete.\n")

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 run_condor.py <input_list.txt>")
        sys.exit(1)

    input_list = sys.argv[1]
    if not os.path.exists(input_list):
        print(f"❌ Error: input file '{input_list}' not found.")
        sys.exit(1)

    print(f"=== Condor job submission ===")
    print(f"Input file: {input_list}")
    print(f"CMSSW base: {CMSSW_BASE}")
    print(f"User: {getpass.getuser()}")
    print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    make_tarball()
    jobs = read_jobs(input_list)
    write_submission(jobs)
    submit()

if __name__ == "__main__":
    main()

