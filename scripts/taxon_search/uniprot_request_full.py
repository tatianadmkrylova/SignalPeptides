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

    all_rows = [] ### List with dictionaries inside; each dictionary corresponds to one protein
    max_signal_count = 0 ### Max number of signal peptide annotations among all proteins
    max_trans_count = 0 ### Max number of transmembrane domain annotations among all proteins

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

            all_rows.append({ #### Save empty information for this protein
                "ProteinId": acc,
                "Length_UniProt": "",
                "PrimaryAccession": "",
                "SecondaryAccessions": "",
                "Organism_class": "",
                "Organism": "",
                "TaxonID_org": "",
                "Organism_host": "",
                "TaxonID_org_host": "",
                "signals": [],
                "trans": []
            })

            continue

        data = r.json()

        organism = ""
        length_UniProt = ""
        primaryAccession = ""
        secondaryAccessions = ""
        taxonID_org = ""
        organism_host = ""
        taxonID_org_host = ""
        organism_class = ""

        organism_data = data.get("organism", {})
        organism = organism_data.get("scientificName", "")
        length_UniProt_data = data.get("sequence", {})
        length_UniProt = length_UniProt_data.get("length", "")
        taxonID_org = organism_data.get("taxonId", "")
        organism_class = organism_data.get("lineage", [""])[0] if organism_data.get("lineage") else ""
        primaryAccession = data.get("primaryAccession", "")
        secondaryAccessions = data.get("secondaryAccessions", []) 
        secondaryAccessions_str = ", ".join(secondaryAccessions) ### to get all accession numbers with ","

        organism_hosts = []
        taxonID_org_hosts = []

        for h in data.get("organismHosts", []): ### it is possible to have many hosts 
            name = h.get("scientificName", "")
            taxon_id = h.get("taxonId", "")

            if name and name not in organism_hosts:
                organism_hosts.append(name)

            if taxon_id and str(taxon_id) not in taxonID_org_hosts:
                taxonID_org_hosts.append(str(taxon_id))

        organism_host = ", ".join(organism_hosts) ### to get all organism_host with ","
        taxonID_org_host = ", ".join(taxonID_org_hosts) ### to get all taxon ID organism_host with ","

        if not taxonID_org:
            print(f"{acc}: no organism taxon ID")

        if not organism_host:
            print(f"{acc}: host field empty")

        if not taxonID_org_host:
            print(f"{acc}: no organism host taxon ID")

        if not length_UniProt:
            print(f"{acc}: no length of protein")

        signal_features = [] #### in case if there are several signal peptide annotations
        trans_features = [] #### in case if there are several transmembrane domain annotations

        for feature in data.get("features", []):

            if feature.get("type") == "Signal":
                location = feature.get("location", {})

                start = location.get("start", {})
                end = location.get("end", {})

                pos_sp_start = start.get("value", "")
                pos_sp_end = end.get("value", "")

                evidence_codes = [] ### if there are several ECO evidences

                for ev in feature.get("evidences", []):
                    code = ev.get("evidenceCode", "")
                    if code and code not in evidence_codes:
                        evidence_codes.append(code)

                evidence_code = ";".join(evidence_codes)

                signal_features.append({ ### every signal peptide annotation has start, end and evidence code
                    "start": pos_sp_start,
                    "end": pos_sp_end,
                    "evidence": evidence_code
                })

            if feature.get("type") == "Transmembrane":
                location = feature.get("location", {})

                start_tr = location.get("start", {})
                end_tr = location.get("end", {})

                pos_tr_start = start_tr.get("value", "")
                pos_tr_end = end_tr.get("value", "")

                evidence_codes_tr = [] ### if there are several ECO evidences

                for ev in feature.get("evidences", []):
                    code = ev.get("evidenceCode", "")
                    if code and code not in evidence_codes_tr:
                        evidence_codes_tr.append(code)

                evidence_code_tr = ";".join(evidence_codes_tr)

                trans_features.append({ ### every transmembrane domain annotation has start, end and evidence code
                    "start": pos_tr_start,
                    "end": pos_tr_end,
                    "evidence": evidence_code_tr
                })

        if not signal_features:
            print(f"{acc}: no signal peptide")

        if not trans_features:
            print(f"{acc}: no transmembrane domain")

        for signal in signal_features:
            if not signal["start"]:
                print(f"{acc}: signal peptide without start position")

            if not signal["end"]:
                print(f"{acc}: signal peptide without end position")

            if not signal["evidence"]:
                print(f"{acc}: no evidence for signal peptide")

        for trans in trans_features:
            if not trans["start"]:
                print(f"{acc}: transmembrane domain without start position")

            if not trans["end"]:
                print(f"{acc}: transmembrane domain without end position")

            if not trans["evidence"]:
                print(f"{acc}: no evidence for transmembrane domain")

        max_signal_count = max(max_signal_count, len(signal_features)) ### to know the max number of signal peptide annotations
        max_trans_count = max(max_trans_count, len(trans_features)) ### to know the max number of transmembrane domain annotations

        all_rows.append({ #### Save all information into the list
            "ProteinId": acc,
            "Length_protein_UniProt":length_UniProt, 
            "PrimaryAccession": primaryAccession,
            "SecondaryAccessions": secondaryAccessions_str,
            "Organism_class": organism_class,
            "Organism": organism,
            "TaxonID_org": taxonID_org,
            "Organism_host": organism_host,
            "TaxonID_org_host": taxonID_org_host,
            "signals": signal_features,
            "trans": trans_features
        })

    header = [ ### Headers of the output file
        "ProteinId",
        "Length_protein_UniProt",
        "PrimaryAccession",
        "SecondaryAccessions",
        "Organism_class",
        "Organism",
        "TaxonID_org",
        "Organism_host",
        "TaxonID_org_host"
    ]

    for i in range(1, max_signal_count + 1): ### To write every signal peptide annotation in a separate column 
        header.extend([
            f"Pos_sp_start_uniprot_{i}",
            f"Pos_sp_end_uniprot_{i}",
            f"Evidence_signal_peptide_{i}"
        ])

    for i in range(1, max_trans_count + 1): ### To write every transmembrane domain annotation in a separate column 
        header.extend([
            f"Pos_tr_start_uniprot_{i}",
            f"Pos_tr_end_uniprot_{i}",
            f"Evidence_code_tr_{i}"
        ])

    writer.writerow(header)

    for item in all_rows: ### Write every protein into the output table
        row_out = [
            item["ProteinId"],
            item["Length_protein_UniProt"],
            item["PrimaryAccession"],
            item["SecondaryAccessions"],
            item["Organism_class"],
            item["Organism"],
            item["TaxonID_org"],
            item["Organism_host"],
            item["TaxonID_org_host"]
        ]

        signals = item["signals"]

        for signal in signals:
            row_out.extend([
                signal["start"],
                signal["end"],
                signal["evidence"]
            ])

        missing_signals = max_signal_count - len(signals)

        for _ in range(missing_signals):
            row_out.extend(["", "", ""])

        trans = item["trans"]

        for tran in trans:
            row_out.extend([
                tran["start"],
                tran["end"],
                tran["evidence"]
            ])

        missing_trans = max_trans_count - len(trans)

        for _ in range(missing_trans):
            row_out.extend(["", "", ""])

        writer.writerow(row_out)


print(f"Saved to {output_file}")