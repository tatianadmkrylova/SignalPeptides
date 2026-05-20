#!/usr/bin/env python3
# To make a fasta file which contains only confirmed sequences (check if proteinID in a column "confirmed")
import sys
import csv

csv_file = sys.argv[1]
fasta_file = sys.argv[2]
output_file = "sp_all_confirmed.fasta"

with open(csv_file, newline="", encoding="utf-8") as f:
    confirmed = {row["ProteinId"].strip() for row in csv.DictReader(f)}

with open(fasta_file, encoding="utf-8") as fasta, open(output_file, "w", encoding="utf-8") as out:
    acc = None
    seq = []

    for line in fasta:
        line = line.strip()

        if line.startswith(">"):
            if acc in confirmed:
                out.write(f">{acc}|SP_confirmed\n{''.join(seq)}\n")

            acc = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)

    if acc in confirmed:
        out.write(f">{acc}|SP_confirmed\n{''.join(seq)}\n")

print(f"Saved to {output_file}")