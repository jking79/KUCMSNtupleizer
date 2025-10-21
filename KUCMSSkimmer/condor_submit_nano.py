#!/usr/bin/env python3
"""
HTCondor Job Submission Script for NanoSV Processing
Creates a single submit file that queues multiple jobs using HTCondor's Queue mechanism.
Auto-creates config.tgz with macro and job-specific file lists.
Automatically submits jobs unless --dry-run is specified.
Organizes output files into logical subdirectories.
"""

import os
import sys
import argparse
import math
import tarfile
import shutil
import subprocess
from pathlib import Path

def create_output_structure(base_output_dir):
    """
    Create organized directory structure for job outputs.
    
    Args:
        base_output_dir (str): Base output directory path
    
    Returns:
        dict: Dictionary with paths to each subdirectory
    """
    base_path = Path(base_output_dir)
    
    # Define subdirectories
    subdirs = {
        'configs': base_path / 'configs',      # Config files (.tgz)
        'logs': base_path / 'logs',            # HTCondor logs (.log)
        'stdout': base_path / 'stdout',        # Standard output (.stdout)
        'stderr': base_path / 'stderr',        # Standard error (.stderr)
        'outputs': base_path / 'outputs',      # ROOT output files (.root)
        'scripts': base_path / 'scripts',      # Job scripts and submit files
        'temp': base_path / 'temp'             # Temporary files (cleaned up)
    }
    
    # Create all directories
    for name, path in subdirs.items():
        path.mkdir(parents=True, exist_ok=True)
        print(f"Created directory: {path}")
    
    return {name: str(path) for name, path in subdirs.items()}

def create_file_chunks(input_file_list, num_jobs=None, files_per_job=None):
    """
    Split the input file list into chunks for parallel processing.
    
    Args:
        input_file_list (str): Path to the text file containing list of files
        num_jobs (int, optional): Number of condor jobs to create
        files_per_job (int, optional): Number of files per job
    
    Returns:
        tuple: (list of file chunks, actual number of jobs)
    """
    with open(input_file_list, 'r') as f:
        files = [line.strip() for line in f if line.strip()]
    
    if not files:
        raise ValueError("Input file list is empty")
    
    if files_per_job:
        # Calculate number of jobs based on files per job
        actual_jobs = math.ceil(len(files) / files_per_job)
        chunk_size = files_per_job
    elif num_jobs:
        # Calculate files per job based on number of jobs
        actual_jobs = min(num_jobs, len(files))  # Don't create more jobs than files
        chunk_size = math.ceil(len(files) / actual_jobs)
    else:
        raise ValueError("Either num_jobs or files_per_job must be specified")
    
    # Split files into chunks
    chunks = []
    for i in range(0, len(files), chunk_size):
        chunk = files[i:i + chunk_size]
        chunks.append(chunk)
    
    return chunks, len(chunks)

def create_config_tgz(macro_path, file_chunk, job_id, dirs):
    """
    Create config.tgz for a specific job containing the macro and file list.
    
    Args:
        macro_path (str): Path to the nanoSVStudy.C macro
        file_chunk (list): List of files for this job
        job_id (int): Job identifier
        dirs (dict): Dictionary with directory paths
    
    Returns:
        str: Relative path to the created config.tgz file (from base output dir)
    """
    # Create temporary config directory
    config_dir = os.path.join(dirs['temp'], f"temp_config_job_{job_id}")
    os.makedirs(config_dir, exist_ok=True)
    
    try:
        # Copy macro to config directory
        macro_dest = os.path.join(config_dir, "nanoSVStudy.C")
        shutil.copy2(macro_path, macro_dest)
        
        # Create file list for this job
        file_list_path = os.path.join(config_dir, "job_files.txt")
        with open(file_list_path, 'w') as f:
            for file_path in file_chunk:
                f.write(f"{file_path}\n")
        
        # Create config.tgz in configs directory
        config_tgz_name = f"config_job_{job_id}.tgz"
        config_tgz_path = os.path.join(dirs['configs'], config_tgz_name)
        with tarfile.open(config_tgz_path, 'w:gz') as tar:
            # Add macro
            tar.add(macro_dest, arcname="nanoSVStudy.C")
            # Add file list
            tar.add(file_list_path, arcname="job_files.txt")
        
        # Return relative path from base output directory
        return f"configs/{config_tgz_name}"
    
    finally:
        # Clean up temporary directory
        shutil.rmtree(config_dir, ignore_errors=True)

def create_job_list_file(file_chunks, dirs):
    """
    Create a job list file for HTCondor Queue mechanism.
    Format: job_id config_file_name (just filename, not path)
    
    Args:
        file_chunks (list): List of file chunks
        dirs (dict): Directory paths
    
    Returns:
        str: Name of the created job list file (relative to scripts dir)
    """
    job_list_path = os.path.join(dirs['scripts'], "job_list.txt")
    with open(job_list_path, 'w') as f:
        for job_id, _ in enumerate(file_chunks):
            # Just use the filename - HTCondor will transfer it to job working directory
            config_name = f"config_job_{job_id}.tgz"
            f.write(f"{job_id} {config_name}\n")
    
    return "job_list.txt"

def create_bash_script_template(dirs):
    """
    Create a single bash script template that works for all jobs.
    
    Args:
        dirs (dict): Directory paths
    
    Returns:
        str: Name of the created bash script (relative to scripts dir)
    """
    bash_content = """#!/bin/bash

# Get job parameters from HTCondor
JOB_ID=$1
CONFIG_FILE=$2

echo "Starting job ${JOB_ID} on $(hostname) at $(date)"
echo "Config file: ${CONFIG_FILE}"
echo "Working directory: $(pwd)"
echo "Available files in working directory:"
ls -la

# Verify config file exists in working directory
if [ ! -f "${CONFIG_FILE}" ]; then
    echo "ERROR: Config file ${CONFIG_FILE} not found in working directory!"
    echo "Available files:"
    ls -la
    exit 1
fi

# Set up CMSSW environment
export SCRAM_ARCH=el8_amd64_gcc12
source /cvmfs/cms.cern.ch/cmsset_default.sh
scramv1 project CMSSW_13_3_0

# Move config to CMSSW directory
echo "Moving config file ${CONFIG_FILE} to CMSSW directory"
mv ${CONFIG_FILE} CMSSW_13_3_0/src/

cd CMSSW_13_3_0/src/
eval `scramv1 runtime -sh`

# Extract configuration (macro + file list)
CONFIG_BASENAME=$(basename ${CONFIG_FILE})
echo "Extracting ${CONFIG_BASENAME}"
tar -xzf ${CONFIG_BASENAME}

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to extract ${CONFIG_BASENAME}"
    exit 1
fi

echo "Working directory content after extraction:"
ls -la
echo "Files to process:"
if [ -f "job_files.txt" ]; then
    cat job_files.txt
else
    echo "ERROR: job_files.txt not found after extraction!"
    exit 1
fi

# Verify macro exists
if [ ! -f "nanoSVStudy.C" ]; then
    echo "ERROR: nanoSVStudy.C not found after extraction!"
    exit 1
fi

# Run ROOT analysis
echo "Processing files for job ${JOB_ID}"
root -l -b -q -e '.L nanoSVStudy.C+' -e 'nanoSVStudy("job_files.txt", "output_job_'${JOB_ID}'.root")'

if [ $? -ne 0 ]; then
    echo "ERROR: ROOT analysis failed!"
    exit 1
fi

# Verify output file was created
if [ ! -f "output_job_${JOB_ID}.root" ]; then
    echo "ERROR: Output file output_job_${JOB_ID}.root was not created!"
    exit 1
fi

# Transfer output back to the initial directory
echo "Transferring output file back"
mv output_job_${JOB_ID}.root ${_CONDOR_SCRATCH_DIR}

echo "Job ${JOB_ID} completed successfully at $(date)"
"""
    
    bash_filename = "nanoSVs_job.sh"
    bash_filepath = os.path.join(dirs['scripts'], bash_filename)
    
    with open(bash_filepath, 'w') as f:
        f.write(bash_content)
    
    # Make script executable
    os.chmod(bash_filepath, 0o755)
    
    return bash_filename

def create_condor_submit_file(job_list_file, bash_script, num_jobs, dirs, memory=2000, disk=1000000, cpus=1):
    """
    Create a single HTCondor submit file that queues multiple jobs.
    
    Args:
        job_list_file (str): Name of the job list file
        bash_script (str): Name of the bash script
        num_jobs (int): Number of jobs
        dirs (dict): Directory paths
        memory (int): Memory request in MB
        disk (int): Disk request in KB
        cpus (int): CPU request
    
    Returns:
        str: Name of the created submit file
    """
    submit_content = f"""# HTCondor Submit File for NanoSV Processing
# Submits {num_jobs} jobs using Queue mechanism
# Files are organized in subdirectories for better organization

universe = vanilla
executable = {bash_script}

# Transfer config files using absolute paths from submit directory
transfer_input_files = ../configs/$(config_file)
should_transfer_files = YES
when_to_transfer_output = ON_EXIT

# Output files organized in subdirectories
transfer_output_files = output_job_$(job_id).root
transfer_output_remaps = "output_job_$(job_id).root = ../outputs/output_job_$(job_id).root"

# Resource requests
request_memory = {memory}
request_disk = {disk}
request_cpus = {cpus}

# Job-specific arguments: job_id and just the config filename (transferred to job working dir)
arguments = $(job_id) $(config_file)

# Output files organized in subdirectories
output = ../stdout/nanoSV_job_$(job_id)_$(Process).stdout
error = ../stderr/nanoSV_job_$(job_id)_$(Process).stderr
log = ../logs/nanoSV_job_$(job_id)_$(Process).log

# Queue jobs from job list file
# Format: job_id config_file
queue job_id,config_file from {job_list_file}
"""
    
    submit_filename = "nanoSVs_submit.jdl"
    submit_filepath = os.path.join(dirs['scripts'], submit_filename)
    
    with open(submit_filepath, 'w') as f:
        f.write(submit_content)
    
    return submit_filename

def submit_jobs(submit_file, scripts_dir):
    """
    Submit the HTCondor jobs using condor_submit.
    
    Args:
        submit_file (str): Name of the submit file
        scripts_dir (str): Scripts directory path
    
    Returns:
        bool: True if submission successful, False otherwise
    """
    submit_path = os.path.join(scripts_dir, submit_file)
    
    # First, let's check if condor_submit is available and get its path
    try:
        # Check if condor_submit exists and get its location
        which_result = subprocess.run(
            ['which', 'condor_submit'],
            capture_output=True,
            text=True,
            check=True
        )
        condor_submit_path = which_result.stdout.strip()
        
        # Check if it's executable
        if not os.access(condor_submit_path, os.X_OK):
            print(f"✗ Error: {condor_submit_path} is not executable")
            return False
            
    except subprocess.CalledProcessError:
        print("✗ Error: condor_submit command not found in PATH.")
        print("Please ensure HTCondor is installed and in your PATH.")
        print("You can check with: which condor_submit")
        return False
    
    # Check if this is a Python script that needs special handling
    is_python_script = False
    try:
        with open(condor_submit_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('exec env LD_LIBRARY_PATH='):
                is_python_script = True
    except:
        pass
    
    try:
        print(f"\nSubmitting jobs to HTCondor...")
        print(f"Command: condor_submit {submit_file}")
        print(f"Working directory: {scripts_dir}")
        
        # Try different execution methods based on the script type
        if is_python_script:
            # For Python scripts with exec env, use shell=True
            print("Using shell execution for Python-based condor_submit")
            result = subprocess.run(
                f'condor_submit {submit_file}',
                cwd=scripts_dir,
                capture_output=True,
                text=True,
                shell=True,
                check=True
            )
        else:
            # Standard execution
            result = subprocess.run(
                [condor_submit_path, submit_file],
                cwd=scripts_dir,
                capture_output=True,
                text=True,
                check=True
            )
        
        print("✓ Jobs submitted successfully!")
        print(result.stdout)
        
        if result.stderr:
            print("Warnings:")
            print(result.stderr)
        
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"✗ Error submitting jobs:")
        print(f"Return code: {e.returncode}")
        print(f"Error output: {e.stderr}")
        if e.stdout:
            print(f"Standard output: {e.stdout}")
        
        # If the first method failed, try alternative approaches
        print("\nTrying alternative execution methods...")
        
        # Try with bash explicitly
        try:
            print("Attempting with bash -c...")
            result = subprocess.run(
                ['bash', '-c', f'condor_submit {submit_file}'],
                cwd=scripts_dir,
                capture_output=True,
                text=True,
                check=True
            )
            print("✓ Jobs submitted successfully with bash!")
            print(result.stdout)
            return True
        except subprocess.CalledProcessError as e2:
            print(f"Bash execution also failed: {e2}")
        
        # Try executing with python3 directly if it's a Python script
        if is_python_script:
            try:
                print("Attempting direct Python execution...")
                result = subprocess.run(
                    ['python3', condor_submit_path, submit_file],
                    cwd=scripts_dir,
                    capture_output=True,
                    text=True,
                    check=True
                )
                print("✓ Jobs submitted successfully with Python!")
                print(result.stdout)
                return True
            except subprocess.CalledProcessError as e3:
                print(f"Python execution also failed: {e3}")
        
        return False
        
    except OSError as e:
        print(f"✗ OS Error when trying to execute condor_submit:")
        print(f"Error: {e}")
        print("This might be an architecture mismatch or permission issue.")
        
        # Additional debugging information
        try:
            # Check file type
            file_result = subprocess.run(
                ['file', condor_submit_path],
                capture_output=True,
                text=True
            )
            print(f"File type: {file_result.stdout.strip()}")
            
            # Check if it's a script that might have wrong shebang
            with open(condor_submit_path, 'rb') as f:
                first_bytes = f.read(100)
                print(f"First bytes: {first_bytes}")
                
        except Exception as debug_e:
            print(f"Could not get additional debug info: {debug_e}")
        
        return False
    except FileNotFoundError:
        print("✗ Error: condor_submit command not found.")
        print("Please ensure HTCondor is installed and in your PATH.")
        return False

def create_submit_script(submit_file, scripts_dir, dry_run=False):
    """
    Create a simple script to submit the job.
    
    Args:
        submit_file (str): Name of the submit file
        scripts_dir (str): Scripts directory path
        dry_run (bool): If True, create script with dry-run flag
    """
    if dry_run:
        submit_content = f"""#!/bin/bash
# DRY RUN - Submit NanoSV processing jobs
echo "DRY RUN: Would submit jobs with:"
echo "condor_submit {submit_file}"
echo ""
echo "Running dry run..."
condor_submit --dry-run {submit_file}
echo ""
echo "To actually submit jobs, run:"
echo "condor_submit {submit_file}"
"""
    else:
        submit_content = f"""#!/bin/bash
# Submit NanoSV processing jobs
echo "Submitting all jobs..."
condor_submit {submit_file}
echo ""
echo "Jobs submitted! Check status with:"
echo "condor_q"
"""
    
    script_name = "submit_jobs.sh"
    script_path = os.path.join(scripts_dir, script_name)
    
    with open(script_path, 'w') as f:
        f.write(submit_content)
    
    os.chmod(script_path, 0o755)
    
    if dry_run:
        print(f"Created DRY RUN submit script: {script_path}")
    else:
        print(f"Created submit script: {script_path}")

def create_summary_file(file_chunks, dirs, total_files):
    """
    Create a summary file showing the job distribution and file organization.
    
    Args:
        file_chunks (list): List of file chunks
        dirs (dict): Directory paths
        total_files (int): Total number of files
    """
    summary_path = os.path.join(dirs['scripts'], "job_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("HTCondor Job Summary\n")
        f.write("===================\n\n")
        f.write(f"Total files: {total_files}\n")
        f.write(f"Total jobs: {len(file_chunks)}\n")
        f.write(f"Average files per job: {total_files / len(file_chunks):.1f}\n\n")
        
        f.write("Directory Structure:\n")
        f.write("├── configs/         # Configuration files (.tgz)\n")
        f.write("├── logs/            # HTCondor log files (.log)\n") 
        f.write("├── stdout/          # Standard output files (.stdout)\n")
        f.write("├── stderr/          # Standard error files (.stderr)\n")
        f.write("├── outputs/         # ROOT output files (.root)\n")
        f.write("├── scripts/         # Job scripts and submit files\n")
        f.write("└── temp/            # Temporary files (cleaned up)\n\n")
        
        f.write("Job Distribution:\n")
        for job_id, chunk in enumerate(file_chunks):
            f.write(f"Job {job_id:3d}: {len(chunk):4d} files -> configs/config_job_{job_id}.tgz\n")
        
        f.write(f"\nSubmit file: scripts/nanoSVs_submit.jdl\n")
        f.write(f"Output files: outputs/output_job_0.root, outputs/output_job_1.root, ..., outputs/output_job_{len(file_chunks)-1}.root\n")
        
        f.write(f"\nLog files:\n")
        f.write(f"  HTCondor logs: logs/nanoSV_job_*.log\n")
        f.write(f"  Standard output: stdout/nanoSV_job_*.stdout\n")
        f.write(f"  Standard error: stderr/nanoSV_job_*.stderr\n")
    
    print(f"Created job summary: {summary_path}")

def cleanup_temp_files(dirs):
    """
    Clean up temporary files and directories.
    
    Args:
        dirs (dict): Directory paths
    """
    temp_dir = Path(dirs['temp'])
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
        print(f"Cleaned up temporary directory: {temp_dir}")

def main():
    parser = argparse.ArgumentParser(
        description="Generate and submit HTCondor jobs using Queue mechanism for NanoSV processing with organized output structure",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create and submit 10 jobs
  python3 condor_submit_nano.py -f file_list.txt -j 10 -m nanoSVStudy.C -o jobs_output
  
  # Create and submit jobs with 50 files each
  python3 condor_submit_nano.py -f file_list.txt --files-per-job 50 -m nanoSVStudy.C -o jobs_output
  
  # Dry run to test (no submission)
  python3 condor_submit_nano.py -f file_list.txt -j 20 -m nanoSVStudy.C -o jobs_output --dry-run
  
  # Manual submission after dry-run
  cd jobs_output/scripts && condor_submit nanoSVs_submit.jdl

Output Structure:
  jobs_output/
  ├── configs/         # Configuration files (.tgz)
  ├── logs/            # HTCondor log files (.log)
  ├── stdout/          # Standard output files (.stdout)
  ├── stderr/          # Standard error files (.stderr)
  ├── outputs/         # ROOT output files (.root)
  ├── scripts/         # Job scripts and submit files
  └── temp/            # Temporary files (cleaned up)
        """
    )
    
    parser.add_argument('-f', '--file-list', required=True,
                        help='Text file containing list of input files to process')
    
    # Job specification - either number of jobs OR files per job
    job_group = parser.add_mutually_exclusive_group(required=True)
    job_group.add_argument('-j', '--jobs', type=int,
                          help='Number of HTCondor jobs to create')
    job_group.add_argument('--files-per-job', type=int,
                          help='Number of files to process per job')
    
    parser.add_argument('-m', '--macro', required=True,
                        help='Path to nanoSVStudy.C macro file')
    parser.add_argument('-o', '--output-dir', required=True,
                        help='Base output directory for organized job files and results')
    parser.add_argument('--memory', type=int, default=2000,
                        help='Memory request per job in MB (default: 2000)')
    parser.add_argument('--disk', type=int, default=1000000,
                        help='Disk request per job in KB (default: 1000000)')
    parser.add_argument('--cpus', type=int, default=1,
                        help='CPU request per job (default: 1)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Create all job files but do not submit jobs (test mode)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.file_list):
        print(f"Error: File list '{args.file_list}' not found")
        sys.exit(1)
    
    if not os.path.exists(args.macro):
        print(f"Error: Macro file '{args.macro}' not found")
        sys.exit(1)
    
    if args.jobs and args.jobs <= 0:
        print("Error: Number of jobs must be positive")
        sys.exit(1)
    
    if args.files_per_job and args.files_per_job <= 0:
        print("Error: Files per job must be positive")
        sys.exit(1)
    
    # Create organized output directory structure
    print("Creating organized directory structure...")
    dirs = create_output_structure(args.output_dir)
    
    # Read total number of files for reporting
    with open(args.file_list, 'r') as f:
        total_files = len([line.strip() for line in f if line.strip()])
    
    if args.jobs:
        print(f"\nCreating {args.jobs} HTCondor jobs for {total_files} files...")
        print(f"Approximately {total_files / args.jobs:.1f} files per job")
    else:
        estimated_jobs = math.ceil(total_files / args.files_per_job)
        print(f"\nCreating jobs with {args.files_per_job} files each for {total_files} total files...")
        print(f"Estimated {estimated_jobs} jobs will be created")
    
    print(f"Input file list: {args.file_list}")
    print(f"Macro file: {args.macro}")
    print(f"Base output directory: {args.output_dir}")
    print(f"Using HTCondor Queue mechanism with organized file structure")
    
    if args.dry_run:
        print("DRY RUN MODE: Jobs will be prepared but not submitted")
    else:
        print("AUTOMATIC SUBMISSION: Jobs will be submitted after creation")
    
    try:
        # Split files into chunks
        file_chunks, actual_jobs = create_file_chunks(
            args.file_list, 
            num_jobs=args.jobs, 
            files_per_job=args.files_per_job
        )
        
        print(f"\nCreated {actual_jobs} jobs with file distribution:")
        for i, chunk in enumerate(file_chunks):
            print(f"  Job {i}: {len(chunk)} files")
        
        # Create config.tgz files for each job
        print(f"\nCreating config files in configs/ directory...")
        for job_id, chunk in enumerate(file_chunks):
            config_path = create_config_tgz(args.macro, chunk, job_id, dirs)
            print(f"  Created {config_path}")
        
        # Create job list file for Queue mechanism
        job_list_file = create_job_list_file(file_chunks, dirs)
        print(f"Created job list: scripts/{job_list_file}")
        
        # Create single bash script template
        bash_script = create_bash_script_template(dirs)
        print(f"Created bash script: scripts/{bash_script}")
        
        # Create single submit file
        submit_file = create_condor_submit_file(
            job_list_file, bash_script, actual_jobs, dirs,
            args.memory, args.disk, args.cpus
        )
        print(f"Created submit file: scripts/{submit_file}")
        
        # Create submit script
        create_submit_script(submit_file, dirs['scripts'], args.dry_run)
        
        # Create summary file
        create_summary_file(file_chunks, dirs, total_files)
        
        # Clean up temporary files
        cleanup_temp_files(dirs)
        
        if args.dry_run:
            print(f"\nDRY RUN: Successfully created {actual_jobs} HTCondor jobs!")
            print("No jobs were submitted. Job files are ready for review.")
            print(f"\nTo submit jobs manually:")
            print(f"  cd {args.output_dir}/scripts")
            print(f"  condor_submit {submit_file}")
            print(f"  # OR")
            print(f"  ./submit_jobs.sh")
            print(f"\nOutput files will be organized in:")
            print(f"  Config files: {args.output_dir}/configs/")
            print(f"  ROOT outputs: {args.output_dir}/outputs/")
            print(f"  Log files: {args.output_dir}/logs/")
            print(f"  Stdout: {args.output_dir}/stdout/")
            print(f"  Stderr: {args.output_dir}/stderr/")
        else:
            # AUTOMATICALLY SUBMIT JOBS
            submission_success = submit_jobs(submit_file, dirs['scripts'])
            
            if submission_success:
                print(f"\n✓ Successfully created and submitted {actual_jobs} HTCondor jobs!")
                print(f"Output files will be organized in:")
                print(f"  ROOT outputs: {args.output_dir}outputs/output_job_0.root, ..., output_job_{actual_jobs-1}.root")
                print(f"  Log files: {args.output_dir}logs/")
                print(f"  Stdout: {args.output_dir}stdout/")
                print(f"  Stderr: {args.output_dir}stderr/")
            else:
                print(f"\n✗ Job files created successfully, but automatic submission failed.")
                print(f"You can manually submit the jobs:")
                print(f"  cd {args.output_dir}/scripts")
                print(f"  condor_submit {submit_file}")
                sys.exit(1)
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
