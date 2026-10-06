#!/usr/bin/env python3
# To download a dataset from signalpeptide.de

import sys
import re
import requests
import pandas as pd

name = sys.argv[1]        ### name of class for the output file
name_class = sys.argv[2]  ### for example, viruses, mammalia, bacteria

dfs = []  ### list of datasets for each page

base_url = (
    f"http://www.signalpeptide.de/index.php?"
    f"sess=&m=listspdb_{name_class}&orderby=id&sortdir=asc"
)

### Read the first page to find the total number of proteins
html = requests.get(base_url).text

### Find the total number of records from text like: 1 - 50 (of 13094)
match = re.search(r"\(of\s+(\d+)\)", html)

if not match:
    raise ValueError("Could not find the total number of records on the page")

total_records = int(match.group(1))
print(f"Total records found: {total_records}")

### Number of proteins per page on signalpeptide.de
page_size = 50

for start in range(0, total_records, page_size):
    url = f"{base_url}&start={start}"
    print("reading", url)

    tables = pd.read_html(url)
    df = max(tables, key=len)  ### find the longest table

    ### Delete empty rows
    df = df.dropna(how="all")

    ### Delete extra header-like rows that pandas can read as data
    df = df[
        ~df.astype(str)
        .apply(
            lambda row: row.str.contains(
                "Accession Number|Entry Name|Protein Name",
                regex=True
            ).any(),
            axis=1
        )
    ]

    ### Delete rows without accession number in the first column
    df = df[df.iloc[:, 0].notna()]

    print("rows:", len(df))

    if len(df) == 0:
        break

    dfs.append(df)

result = pd.concat(dfs, ignore_index=True).drop_duplicates()  ### merge dfs, ignore_index=True - enumerate again a dataset after merge, drop_duplicates() - delete the duplicated lines

### Final cleaning after merging all pages
result = result.dropna(how="all")

result = result[
    ~result.astype(str)
    .apply(
        lambda row: row.str.contains(
            "Accession Number|Entry Name|Protein Name",
            regex=True
        ).any(),
        axis=1
    )
]

result = result[result.iloc[:, 0].notna()]

result.to_csv(
    f"../../data/raw/signalpeptide_de/signalpeptide_{name}.csv",
    index=False
)  ### index = False to do not save an index of pandas like a separate column

print(f"Saved: signalpeptide_{name}.csv")
print(f"Total rows: {len(result)}")
