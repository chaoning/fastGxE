import pandas as pd
import numpy as np
import sys, glob, os, re

wk_dir = sys.argv[1]
os.chdir(wk_dir)

# Change to "**/*.signal.txt" and recursive=True if files are in subfolders
files = glob.glob("*.signal.txt")
if not files:
    raise FileNotFoundError(f"No files matched in {wk_dir} with pattern *.signal.txt")

df_lst = []
bad, empty = [], []

for f in files:
    try:
        df = pd.read_csv(f, sep="\t", low_memory=False)
        # Treat all-NA rows as empty
        if df.empty or df.dropna(how="all").empty:
            empty.append(f)
            continue
        # (Optional) keep source filename
        df["__source"] = os.path.basename(f)
        df_lst.append(df)
    except Exception as e:
        bad.append((f, str(e)))

if not df_lst:
    raise ValueError(
        "No non-empty tables to merge. "
        f"Empty: {len(empty)} files, Failed: {len(bad)} files. "
        "Check input pattern/contents."
    )

# Concatenate (outer join on columns so slight schema mismatches won't crash)
df = pd.concat(df_lst, ignore_index=True, sort=False)

# Output (choose one)
df.to_csv("signal.merged.csv", index=False)
# df.to_csv("signal.merged.csv.gz", index=False, compression="gzip")

print(f"Merged {len(df_lst)} / {len(files)} files into shape {df.shape}.")
if empty:
    print(f"Skipped {len(empty)} empty/all-NA files (e.g., {empty[:3]})")
if bad:
    print(f"Failed to read {len(bad)} files (e.g., {bad[:3]})")
