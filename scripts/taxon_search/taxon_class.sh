#!/bin/bash
# To download class of organism using id of taxon (via efetch ncbi)
echo "TaxID,Class"

while read taxid
do
  [ -z "$taxid" ] && continue

  class_name=$(efetch -db taxonomy -id "$taxid" -format xml | \
  xtract -pattern Taxon \
    -block "LineageEx/Taxon" \
      -if Rank -equals "class" \
      -element ScientificName)

  [ -z "$class_name" ] && class_name="NA"

  echo "${taxid},${class_name}"

  sleep 2

done < "$1"
