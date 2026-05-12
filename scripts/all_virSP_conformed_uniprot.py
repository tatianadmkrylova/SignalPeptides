#!/usr/bin/env python3

import csv
import requests
import time

output_file = "SP_allviruses_verified2.csv"

url = "https://rest.uniprot.org/uniprotkb/stream"

params = {
    "query": "taxonomy_id:10239 AND ft_signal:*", ## NCBI Taxonomy id:10239 for viruses
    "format": "json"
}

r = requests.get(url, params=params, timeout=300)
r.raise_for_status()
data = r.json()

total_entries = 0
entries_with_hosts = 0
signal_features = 0
signal_with_selected_eco = 0

with open(output_file, "w", newline="", encoding="utf-8") as out:
    writer = csv.writer(out)
    writer.writerow(["ProteinID", "Organism", "Hosts", "Length", "Feature"])

    for entry in data["results"]:
        total_entries += 1

        acc = entry["primaryAccession"]
        organism = entry["organism"]["scientificName"]
        hosts = entry.get("organism", {}).get("hosts", [])

        entries_with_hosts += 1

        host_names = "; ".join(
            host.get("scientificName", "")
            for host in hosts
        )

        for feature in entry.get("features", []):
            if feature.get("type") in ["Signal peptide", "Signal"]:
                signal_features += 1

                has_eco = any(
                    ev.get("evidenceCode") in ["ECO:0000269", "ECO:0000303", "ECO:0000305"] #https://www.uniprot.org/help/evidences
                    for ev in feature.get("evidences", [])
                )

                if has_eco:
                    signal_with_selected_eco += 1

                    length = feature["location"]["end"]["value"]

                    writer.writerow([acc, organism, host_names, length, feature])
                    print(acc, organism, host_names, length, feature)

        time.sleep(0.2)

print("Total entries:", total_entries)
print("Entries with hosts:", entries_with_hosts)
print("Signal features after host filter:", signal_features)
print("Signal features with selected ECO:", signal_with_selected_eco)
print(f"Saved to {output_file}")