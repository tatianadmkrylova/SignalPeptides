#!/usr/bin/env python3
# To download a full sequence of the protein via UniProt in fasta format
import sys
import csv
import requests
import time

input_file = sys.argv[1] ### Downloaded file from signalpeptide.de
output_file = sys.argv[2] ### Name and type of the output file

with open(input_file, newline="", encoding="utf-8") as f, open(output_file, "w") as out:
    reader = csv.DictReader(f)

    for row in reader:
        acc = row["Accession Number"].strip()

        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        r = requests.get(url)

        time.sleep(0.2)

        if r.status_code != 200: ### Error shows that the request did not successfully execute
            print(f"{acc}: not found")
            continue

        data = r.json() ### Data from JSON (request)

        # If accession was merged/demerged, get the new accession
        if "sequence" not in data and "inactiveReason" in data:
            if "mergeDemergeTo" in data["inactiveReason"]:
                new_acc = data["inactiveReason"]["mergeDemergeTo"][0]

                url = f"https://rest.uniprot.org/uniprotkb/{new_acc}.json"
                r = requests.get(url)

                time.sleep(0.2)

                if r.status_code != 200:
                    print(f"{acc}: redirected to {new_acc}, but not found")
                    continue

                data = r.json()
                acc = new_acc

        if "sequence" not in data: ### Error if there is not a sequence data for this accession number
            print(f"{acc}: no sequence in UniProt")
            continue

        seq = data["sequence"]["value"]

        out.write(f">{acc}\n{seq}\n") ### Fasta format for the output file

print(f"Saved to {output_file}")







