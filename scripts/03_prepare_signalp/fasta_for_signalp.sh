#!/bin/bash
## Script to delete duplicates of sequences in fasta file

input_fasta=$1
output_all_csv=$2
output_unique_csv=$3
output_unique_fasta=$4 ### only unique sequence in the output file

### verification of the arguments
if [ -z "$input_fasta" ] || [ -z "$output_all_csv" ] || [ -z "$output_unique_csv" ] || [ -z "$output_unique_fasta" ]; then
    echo "Usage: $0 input.fasta output_all.csv output_unique.csv output_unique.fasta"
    exit 1
fi

if [ ! -f "$input_fasta" ]; then
    echo "Error: file '$input_fasta' not found"
    exit 1
fi

echo "ProteinId,Sequence" > "$output_all_csv"

id=""
seq=""

while read line
do
    if [[ "$line" == ">"* ]]; then
        if [ -n "$id" ]; then
            echo "${id},${seq}" >> "$output_all_csv"
        fi

        id="${line#>}"
        id="${id%%|*}"
        id="${id%% *}"
        seq=""
    else
        seq="${seq}${line}"
    fi
done < "$input_fasta"

if [ -n "$id" ]; then
    echo "${id},${seq}" >> "$output_all_csv"
fi

awk -F',' 'NR==1 || !seen[$2]++' "$output_all_csv" > "$output_unique_csv"

awk -F',' 'NR > 1 {print ">" $1 "\n" $2}' "$output_unique_csv" > "$output_unique_fasta"