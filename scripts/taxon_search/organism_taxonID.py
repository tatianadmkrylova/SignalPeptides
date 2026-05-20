#!/usr/bin/env python3
# To download a taxId for each protein via ProteinID
import sys
import csv
import requests
import time

input_file = sys.argv[1]
output_file = "organism_taxID.csv"

with open(input_file, newline="", encoding="utf-8") as f, \
     open(output_file, "w", newline="", encoding="utf-8") as out:

    reader = csv.DictReader(f)
    writer = csv.writer(out)
    writer.writerow(["ProteinId", "Length", "FullLength", "TaxonID_org"])

    for row in reader:
        acc = row["ProteinId"].strip()
        length_str = row["Length"].strip()

        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        r = requests.get(url, timeout=30)

        time.sleep(0.2)

        if r.status_code != 200:
            print(f"{acc}: not found")
            continue

        data = r.json()

        taxonID_org = data.get("organism", {}).get("taxonId", "")
        full_length = data.get("sequence", {}).get("length", "")

        if not taxonID_org:
            print(f"{acc}: no organism taxon ID")
            continue

        if not full_length:
            print(f"{acc}: no full protein length")
            continue

        writer.writerow([acc, length_str, full_length, taxonID_org])

print(f"Saved to {output_file}")