#!/usr/bin/env bash
set -euo pipefail

curl -L \
"https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid&format=tsv&query=%28%28%28%28keyword%3AKW-0167%29+OR+%28keyword%3AKW-0240%29+OR+%28keyword%3AKW-0696%29+OR+%28keyword%3AKW-0645%29%29+AND+virus+AND+%28virus_host_id%3A40674%29%29+NOT+%28ft_signal%3A*%29%29+AND+%28existence%3A1%29" \
| tr '\t' ',' > ../../data/raw/uniprot/uniprotkb_viruses_without_signalpeptide_newrequest.csv

