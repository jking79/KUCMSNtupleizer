#!/usr/bin/env python3
import subprocess

def get_empty_argument_jobs(user=None):
    """Return a list of ClusterIds of Condor jobs with empty arguments."""
    cmd = ["condor_q", "-long"]
    if user:
        cmd.append(user)
    result = subprocess.run(["bash",cmd], capture_output=True, text=True, check=True)

    cluster_id = None
    empty_jobs = []

    for line in result.stdout.splitlines():
        line = line.strip()
        if line.startswith("ClusterId ="):
            cluster_id = line.split('=')[1].strip()
        elif line.startswith("Arguments =") and cluster_id is not None:
            args = line.split('=', 1)[1].strip().strip('"')
            if args == "":
                empty_jobs.append(cluster_id)
            cluster_id = None  # reset for next job

    return empty_jobs

def remove_jobs(cluster_ids):
    """Remove the given Condor jobs."""
    if not cluster_ids:
        print("No jobs with empty arguments found.")
        return

    print("Jobs with empty arguments found:")
    for cid in cluster_ids:
        print(f"  {cid}")

    confirm = input("Do you want to remove these jobs? [y/N]: ").strip().lower()
    if confirm != 'y':
        print("Aborting. No jobs removed.")
        return

    for cid in cluster_ids:
        subprocess.run(["condor_rm", cid], check=False)
    print("Done. Jobs removed.")

def main():
    import getpass
    user = getpass.getuser()  # default to current user
    empty_jobs = get_empty_argument_jobs(user)
    remove_jobs(empty_jobs)

if __name__ == "__main__":
    main()

