#!/usr/bin/env python3

import sys
import csv
import requests
import time

input_file = sys.argv[1]
output_file = "hosts_taxID.csv"

with open(input_file, newline="", encoding="utf-8") as f, \
     open(output_file, "w", newline="", encoding="utf-8") as out:

    reader = csv.DictReader(f)
    writer = csv.writer(out)
    writer.writerow(["ProteinId", "Length", "Host", "TaxonID"])

    for row in reader:
        acc = row["ProteinId"].strip()
        length_str = row["Length"].strip()
        host = row["Host"].strip()

        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        r = requests.get(url, timeout=30)

        time.sleep(0.2)

        if r.status_code != 200:
            print(f"{acc}: not found")
            continue

        data = r.json()

        taxonIDs = []
        for t in data["organismHosts"]:
            taxon = t.get("taxonId")
            if taxon:
                taxonIDs.append(str(taxon)) ## taxonID is numeric

        if not taxonIDs:
            print(f"{acc}: host field empty")
            continue

        taxonID = ", ".join(taxonIDs)

        writer.writerow([acc, length_str, host, taxonID])

print(f"Saved to {output_file}")