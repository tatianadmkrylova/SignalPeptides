#!/usr/bin/env python3
# To download a dataset from signalpeptide.de

import sys
import pandas as pd

name = sys.argv[1] ### name of class for the output file
name_class = sys.argv[2] ### for example, viruses, mammalia, bacteria 

dfs = [] ### list of datasets for each page (50 proteins par page by default)

for start in range(0, 100000, 25):
    url = f"http://www.signalpeptide.de/index.php?sess=&m=listspdb_{name_class}&start={start}&orderby=id&sortdir=asc"
    print("reading", url)

    tables = pd.read_html(url)
    df = max(tables, key=len) ### find the most long table (key = len, length of the table)

    print("rows:", len(df))

    if len(df) == 0:
        break

    dfs.append(df)

result = pd.concat(dfs, ignore_index=True).drop_duplicates() ### merge dfs, ignore_index=True - enumerate again a dataset after merge, drop_duplicates() - delete the duplicated lines 

result.to_csv(f"output/signalpeptide_{name}.csv", index=False) ### index = False to do not save an index of pandas like a separate column 

print(f"Saved: signalpeptide_{name}.csv")
print(f"Total rows: {len(result)}")
