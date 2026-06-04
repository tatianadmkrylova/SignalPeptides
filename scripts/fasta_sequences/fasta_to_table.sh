#!/bin/bash

input_fasta=$1 ### input fasta file (with duplicates)
output_all=$2 ### output csv file (with duplicates) in csv format
output_unique=$3 ### output csv file (without duplicates) 

if [ -z "$input_fasta" ] || [ -z "$output_all" ] || [ -z "$output_unique" ]; then ### control if all arguments are present
    echo "Usage: $0 input.fasta output_all.csv output_unique.csv"
    exit 1
fi

if [ ! -f "$input_fasta" ]; then
    echo "Error: file '$input_fasta' not found"  ### control if the first argument (input) is a file
    exit 1
fi

echo "ProteinId,Sequence" > "$output_all" ### the headers

id="" ### initiation of IDs
seq="" ### initiation of sequences

while read line 
do
    if [[ "$line" == ">"* ]]; then ### if line starts with ">" -> proteinId

        if [ -n "$id" ]; then ### if it's not empty (id)
            echo "${id},${seq}" >> "$output_all" ### save to 'output' file
        fi

        id="${line#>}" ### remove ">" in Id 
        id="${id%%|*}" ### remove "|" in Id 
        id="${id%% *}" ### remove " " in Id 
        seq=""

    else
        seq="${seq}${line}" ### if it's not a line (sequence), add to seq
    fi
done < "$input_fasta"

if [ -n "$id" ]; then
    echo "${id},${seq}" >> "$output_all"  ### if it's the last line (after there is no a ">")
fi

awk -F',' 'NR==1 || !seen[$2]++' "$output_all" > "$output_unique" ### to filter just unique sequences 