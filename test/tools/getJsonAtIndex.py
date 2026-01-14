#!/usr/bin/env python3
import json
import sys

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 extract_json_index.py <json_file> <index_key> [output.txt]")
        sys.exit(1)

    json_file = sys.argv[1]
    index_key = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else None

    # --- Load the JSON file ---
    try:
        with open(json_file, "r") as f:
            data = json.load(f)
    except Exception as e:
        print(f"❌ Error reading {json_file}: {e}")
        sys.exit(1)

    # --- Check for the requested key ---
    if index_key not in data:
        print(f"⚠️ Key '{index_key}' not found in {json_file}")
        print(f"Available keys: {list(data.keys())[:10]}{'...' if len(data) > 10 else ''}")
        sys.exit(1)

    value = data[index_key]

    # --- Print or write to file ---
    if output_file:
        with open(output_file, "w") as out:
            json.dump(value, out, indent=2)
        print(f"✅ Wrote contents of key '{index_key}' to {output_file}")
    else:
        print(f"✅ Contents of key '{index_key}':\n")
        print(json.dumps(value, indent=2))

if __name__ == "__main__":
    main()

