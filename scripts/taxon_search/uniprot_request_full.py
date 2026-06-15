#!/usr/bin/env python3

import sys
import csv
import time
import requests


if len(sys.argv) != 3: #### Error of usage 
    print("Usage: python script.py input.csv output.csv")
    sys.exit(1)


input_file = sys.argv[1] ### File which contains ProteinId 
output_file = sys.argv[2] ### Name and type of the output file


with open(input_file, newline="", encoding="utf-8") as f, \
     open(output_file, "w", newline="", encoding="utf-8") as out:

    reader = csv.DictReader(f) ### Read every line like a dictionnary
    writer = csv.writer(out)

    writer.writerow([ ### Headers of the output file
        "ProteinId",
        "Organism_class",
        "Pos_sp_start_uniprot",
        "Pos_sp_end_uniprot",
        "Organism",
        "TaxonID_org",
        "Organism_host",
        "TaxonID_org_host",
        "Evidence_signal_peptide"
    ])

    for row in reader: ### for line in the input file
        acc = row["ProteinId"].strip() #### Take an accession number

        if not acc: #### If there is not an accession number, continue
            continue

        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json" #### Connect to json Uniprot via accession number (ProteinId)

        try: 
            r = requests.get(url, timeout=30) ### requests to UniProt
            time.sleep(0.2)
            r.raise_for_status() ### is the request is successful ?

        except requests.RequestException as e:
            print(f"{acc}: request failed: {e}") ### Errors with the requests 

            writer.writerow([ #### Write the information into the output table
                acc,
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                ""
            ])

            continue

        data = r.json()

        organism = ""
        taxonID_org = ""
        organism_host = ""
        taxonID_org_host = ""
        organism_class = ""

        organism_data = data.get("organism", {})
        organism = organism_data.get("scientificName", "")
        taxonID_org = organism_data.get("taxonId", "")
        organism_class = data.get("organism", {}).get("lineage", [""])[0]

        organism_hosts = []
        taxonID_org_hosts = []

        for h in data.get("organismHosts", []):
            name = h.get("scientificName", "")
            taxon_id = h.get("taxonId", "")

            if name and name not in organism_hosts:
                organism_hosts.append(name)

            if taxon_id and str(taxon_id) not in taxonID_org_hosts:
                taxonID_org_hosts.append(str(taxon_id))

        organism_host = ", ".join(organism_hosts)
        taxonID_org_host = ", ".join(taxonID_org_hosts)

        if not taxonID_org:
            print(f"{acc}: no organism taxon ID")

        if not organism_host:
            print(f"{acc}: host field empty")

        if not taxonID_org_host:
            print(f"{acc}: no organism host taxon ID")

        signal_features = []

        for feature in data.get("features", []):
            if feature.get("type") == "Signal":
                location = feature.get("location", {})

                start = location.get("start", {})
                end = location.get("end", {})

                pos_sp_start = start.get("value", "")
                pos_sp_end = end.get("value", "")

                evidence_codes = [] ### if there is a several ECO (evidences)

                for ev in feature.get("evidences", []):
                    code = ev.get("evidenceCode", "")
                    if code and code not in evidence_codes:
                        evidence_codes.append(code)

                evidence_code = ";".join(evidence_codes)

                signal_features.append({ ### every evidence code has position start and end of signal peptide
                    "start": pos_sp_start,
                    "end": pos_sp_end,
                    "evidence": evidence_code
                })

        if not signal_features:
            print(f"{acc}: no signal peptide")

            writer.writerow([
                acc,
                organism_class,
                "",
                "",
                organism,
                taxonID_org,
                organism_host,
                taxonID_org_host,
                ""
            ])

            continue

        for signal in signal_features:
            if not signal["start"]:
                print(f"{acc}: signal peptide without start position")

            if not signal["end"]:
                print(f"{acc}: signal peptide without end position")

            if not signal["evidence"]:
                print(f"{acc}: no evidence for signal peptide")

            writer.writerow([
                acc,
                organism_class,
                signal["start"],
                signal["end"],
                organism,
                taxonID_org,
                organism_host,
                taxonID_org_host,
                signal["evidence"]
            ])


print(f"Saved to {output_file}")