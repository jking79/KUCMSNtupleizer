def create_crab_resubmit_script(input_file="crab_status_check.sh", output_file="crab_resubmit.sh"):
    with open(input_file, "r") as input_script:
        lines = input_script.readlines()

    resubmit_lines = [line.replace("status", "resubmit") for line in lines]

    with open(output_file, "w") as output_script:
        output_script.write("".join(resubmit_lines))

if __name__ == "__main__":
    create_crab_resubmit_script()
