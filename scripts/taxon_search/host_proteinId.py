#!/usr/bin/env python3
# To download host common name via ProteinId in UniProt
import sys
import csv
import requests
import time

input_file = sys.argv[1]
output_file = "hosts.csv"

with open(input_file, newline="", encoding="utf-8") as f, \
     open(output_file, "w", newline="", encoding="utf-8") as out:

    reader = csv.DictReader(f)
    writer = csv.writer(out)
    writer.writerow(["ProteinId", "Length", "Host"])

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

        if "organismHosts" not in data:
            print(f"{acc}: no host in UniProt")
            continue

        hosts = []
        for h in data["organismHosts"]:
            name = h.get("scientificName") or h.get("commonName")
            if name:
                hosts.append(name)

        if not hosts:
            print(f"{acc}: host field empty")
            continue

        host = ", ".join(hosts)

        writer.writerow([acc, length_str, host])

print(f"Saved to {output_file}")