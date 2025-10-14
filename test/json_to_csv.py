#!/usr/bin/env python3
"""
Convert a CMS-style lumi JSON file into a CSV file for Excel or Numbers.

Example:
    python3 lumi_json_to_csv.py Cert_Example.json lumi_list.csv
"""

import json
import csv
import sys
import os

def convert_lumi_json_to_csv(input_json, output_csv):
    # Read JSON
    with open(input_json, 'r') as f:
        data = json.load(f)

    # Write CSV
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Run", "LumiStart", "LumiEnd"])
        
        for run, lumiranges in data.items():
            for start, end in lumiranges:
                writer.writerow([run, start, end])

    print(f"✅ Conversion complete: {output_csv}")


def main():
    #if len(sys.argv) != 3:
    #    print("Usage: python3 lumi_json_to_csv.py input.json output.csv")
    #    sys.exit(1)
    
    #input_json = sys.argv[1]
    #output_csv = sys.argv[2]

    input_json = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt' 
    output_csv = 'UL2018_314472-325175_13TeV_GoldJson.csv'

    if not os.path.exists(input_json):
        print(f"❌ Error: Input file '{input_json}' not found.")
        sys.exit(1)

    convert_lumi_json_to_csv(input_json, output_csv)


if __name__ == "__main__":
    main()

