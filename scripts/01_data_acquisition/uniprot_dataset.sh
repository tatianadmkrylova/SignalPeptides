#!/usr/bin/env bash
set -euo pipefail

curl -L \
  "https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:10239%20AND%20ft_signal:*&format=tsv&fields=accession" \
  | tr '\t' ',' > ../../data/raw/uniprot/uniprotkb_taxonomy_viruses_AND_ft_sign_entry.csv

curl -L \
  "https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:40674%20AND%20ft_signal_exp:*&format=tsv&fields=accession" \
  | tr '\t' ',' > ../../data/raw/uniprot/uniprot_mammalia_signal_confirmed_exp.csv

curl -L \
  "https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:2%20AND%20ft_signal_exp:*&format=tsv&fields=accession" \
  | tr '\t' ',' > ../../data/raw/uniprot/uniprot_bacteria_signal_confirmed_exp.csv


### (taxonomy_id:10239) AND (reviewed:true) AND (existence:1) AND NOT (ft_signal:*) ## URL request to collect the proteins without signal peptide in viruses (control dataset) ### 3 346 results
### (taxonomy_id:10239) ### Viruses (taxonID)
### (reviewed:true) ### all reviewed entries
### (existence:1) ### existance at protein level (experimentally)
### In UniProtKB there are 5 types of evidence for the existence of a protein:
### 
### 1. Experimental evidence at protein level
### 2. Experimental evidence at transcript level
### 3. Protein inferred from homology
### 4. Protein predicted
### 5. Protein uncertain
### (ft_signal:*) ### feature: signal peptide 
