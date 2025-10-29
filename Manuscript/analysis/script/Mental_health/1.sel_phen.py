import pandas as pd
import sys

in_file = sys.argv[1]
out_file = sys.argv[2]

column_names = ['eid']
df = pd.read_csv(in_file, sep="\t", usecols=["Field"], dtype=str)
column_names.extend([ f"{field}-0.0" for field in list(df.iloc[:, 0])])

print(column_names)
# Load CSV in chunks, only selecting necessary columns
df_chunk = pd.read_csv('/net/zootopia/disk1/chaon/data/UKB/pheno/ukb673999.csv',
                       chunksize=10000, usecols=column_names, dtype=str)

# Concatenate chunks into one DataFrame
df = pd.concat(df_chunk)

# Save to CSV
df.to_csv(out_file, index=False)
