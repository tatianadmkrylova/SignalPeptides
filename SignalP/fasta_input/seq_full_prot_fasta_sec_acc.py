#!/usr/bin/env python3

import sys
import csv
import requests
import time

input_file = sys.argv[1]   # File with UniProt metadata (.csv)
output_file = sys.argv[2]  # Output FASTA file

with open(input_file, newline="", encoding="utf-8") as f, open(output_file, "w") as out:
    reader = csv.DictReader(f)

    for row in reader:
        input_acc = row["PrimaryAccession"].strip()

        url = f"https://rest.uniprot.org/uniprotkb/{input_acc}.json"
        r = requests.get(url)
        time.sleep(0.2)

        if r.status_code != 200:
            print(f"{input_acc}: not found")
            continue

        data = r.json()

        if "sequence" not in data:
            print(f"{input_acc}: no sequence in UniProt")
            continue

        current_primary_acc = data["primaryAccession"]
        seq = data["sequence"]["value"]

        out.write(f">{current_primary_acc}\n{seq}\n")

print(f"Saved to {output_file}")







