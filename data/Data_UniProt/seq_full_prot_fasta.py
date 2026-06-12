#!/usr/bin/env python3
# To download a full sequence of the protein via UniProt in fasta format
import sys
import csv
import requests
import time

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, newline="", encoding="utf-8") as f, open(output_file, "w") as out:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        acc = row["Accession Number"].strip()
        

        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        r = requests.get(url)

        time.sleep(0.2)

        if r.status_code != 200:
            print(f"{acc}: not found")
            continue

        data = r.json()

        if "sequence" not in data:
            print(f"{acc}: no sequence in UniProt")
            continue

        seq = data["sequence"]["value"]

        out.write(f">{acc}\n{seq}\n")

print(f"Saved to {output_file}")
