#!/usr/bin/env bash
set -euo pipefail

curl -L \
  -o uniprotkb_viruses_taxonomy_id_10239_AND_ft.csv \
  "https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:10239%20AND%20ft_signal:*&format=csv&fields=accession"

curl -L \
  -o uniprotkb_mammalia_taxonomy_id_40674_AND_ft_sign_exp.csv \
  "https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:40674%20AND%20ft_signal_exp:*&format=csv&fields=accession"

curl -L \
  -o uniprotkb_bacteria_taxonomy_id_2_AND_ft_sign_exp.csv \
  "https://rest.uniprot.org/uniprotkb/stream?query=taxonomy_id:2%20AND%20ft_signal_exp:*&format=csv&fields=accession"
