#!/usr/bin/env python3

import sys
import csv
import requests
import time
import re
import os

input_file = sys.argv[1]
output_file = "signal_peptides.fasta"

with open(input_file) as f, open(output_file, "w") as out:
    reader = csv.reader(f, delimiter=",")

    for row in reader:
        if not row: # if row is empty, continue
            continue

        acc = row[0].strip() # strip() to hide the spaces

        if acc.lower() in ["accession", "id", "uniprot"]: # pass the headers
            continue

        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json" #access to uniprot db
        r = requests.get(url)

        time.sleep(0.2)

        if r.status_code != 200: # fault (no protein with the accession number)
            print(f"{acc}: not found")
            continue

        data = r.json()
        seq = data["sequence"]["value"] # full amino acids sequences

        found = False # flag to find a signal peptide

        for feature in data.get("features", []): 
            if feature.get("type") == "Signal peptide":
                start = int(feature["location"]["start"]["value"]) # position start
                end = int(feature["location"]["end"]["value"]) #position end  

                signal_peptide = seq[start - 1:end] 

                out.write(f">{acc}\n{signal_peptide}\n")

                found = True

        if not found:
            print(f"{acc}: no signal peptide")

print(f"Saved to {output_file}")
