import os
import argparse

def create_crab_status_script(folder_path):
    script_lines = []
    
    for folder_name in os.listdir(folder_path):
        full_path = os.path.join(folder_path, folder_name)
        if os.path.isdir(full_path):
            crab_status_line = f'crab status -d {full_path}'
            script_lines.append(crab_status_line)

    with open("crab_status_check.sh", "w") as script_file:
        script_file.write("\n".join(script_lines))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create crab status check script")
    parser.add_argument("folder", help="Path to the folder containing subfolders")

    args = parser.parse_args()
    create_crab_status_script(args.folder)
