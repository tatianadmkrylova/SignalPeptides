#!/usr/bin/env python3

import sys
import csv
import requests
import time

input_file = sys.argv[1]
output_file = "signal_peptides_full.fasta"

with open(input_file, newline="", encoding="utf-8") as f, open(output_file, "w") as out:
    reader = csv.DictReader(f)

    for row in reader:
        acc = row["Accession Number"].strip()
        length_str = row["Length"].strip()

        if not length_str.isdigit():
            print(f"{acc}: invalid length ({length_str})")
            continue

        sp_len = int(length_str)

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

        signal_peptide = seq[:sp_len]

        out.write(f">{acc}|SP={sp_len}\n{signal_peptide}\n")

print(f"Saved to {output_file}")