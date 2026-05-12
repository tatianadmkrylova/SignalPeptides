
#!/bin/bash

while read taxid
do
  [ -z "$taxid" ] && continue

 #echo "Processing TaxID: $taxid" >&2

  efetch -db taxonomy -id "$taxid" -format xml | \
  xtract -pattern Taxon \
    -first TaxId \
    -element TaxId \
    -block "LineageEx/Taxon" \
      -unless Rank -equals "no rank" \
      -tab "," \
      -sep "_" \
      -element Rank,ScientificName \
    -group Taxon \
      -unless Rank -equals "no rank" \
      -tab "," \
      -sep "_" \
      -element Rank,ScientificName

  sleep 0.5

done < "$1"
