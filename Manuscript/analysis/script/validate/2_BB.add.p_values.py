import pandas as pd
import os
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

# ------------------------- Paths -------------------------
input_file = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/validate/BB_known.txt"
output_file = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/validate/BB_known.addp.csv"
data_dir = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/testGxE"

# ------------------------- Read input -------------------------
df = pd.read_csv(input_file, sep=r"\s+")
# Assume the first column is fieldID, the second column is SNP
col_field = df.columns[0]
col_snp = df.columns[1]

# ------------------------- Worker function -------------------------
def process_field_group(fieldID: str, group: pd.DataFrame):
    """
    For each fieldID:
      - Read 10 shard files: {fieldID}.10_1.res ... {fieldID}.10_10.res
      - Concatenate them into one DataFrame (rows are SNPs)
      - Reindex by SNPs from the input group to preserve order
      - Return a DataFrame aligned with the original rows
    """
    frames = []
    for i in range(1, 11):
        file_path = f"{data_dir}/{fieldID}.10_{i}.res"
        if not os.path.exists(file_path):
            # Skip if shard file does not exist
            continue
        try:
            df2 = pd.read_csv(file_path, sep=r"\s+")
        except Exception as e:
            print(f"[WARN] Failed to read {file_path}: {e}")
            continue

        if "SNP" not in df2.columns:
            print(f"[WARN] No 'SNP' column in {file_path}, skipped.")
            continue

        df2 = df2.set_index("SNP")
        frames.append(df2)

    if not frames:
        # If no shards found, return an empty DataFrame aligned with group
        return group.index, pd.DataFrame(index=group.index)

    # Concatenate all shards (no duplicate SNPs expected)
    big = pd.concat(frames, axis=0, copy=False)

    # Align to the SNPs in the input group (preserve order)
    sub_df = big.reindex(group[col_snp].values).reset_index()  # restore SNP as column
    sub_df.index = group.index  # align row indices with the original input
    return group.index, sub_df

# ------------------------- Parallel execution -------------------------
field_groups = df.groupby(col_field, sort=False)  # keep original order
max_workers = min(32, (os.cpu_count() or 8))
results = []

with ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = {
        executor.submit(process_field_group, fg[0], fg[1]): fg[0]
        for fg in field_groups
    }
    for fut in as_completed(futures):
        try:
            results.append(fut.result())
        except Exception as e:
            fid = futures[fut]
            print(f"[ERROR] Field {fid} failed: {e}")

# ------------------------- Merge & Save -------------------------
if results:
    merged = pd.concat([r[1] for r in results], axis=0)
    df_out = pd.concat([df, merged], axis=1)
else:
    df_out = df.copy()

df_out.to_csv(output_file, index=False)
print(f"Updated file saved to {output_file}")
