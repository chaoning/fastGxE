import pandas as pd
import os
import numpy as np
from concurrent.futures import ThreadPoolExecutor

# Input and output file paths
input_file = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/validate/PT_known.txt"
output_file = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/validate/PT_known.addp.csv"
data_dir = "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/merge/"

# Read the input file
df = pd.read_csv(input_file, sep=r"\s+")
df["p"] = np.nan  # Add a column for p-values

# Group rows by fieldID
field_groups = df.groupby(df.columns[0])  # Assuming first column is `fieldID`

# Function to process each fieldID group
def process_field_group(fieldID, group):
    file_path = f"{data_dir}/{fieldID}.gz"
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return group.index, [np.nan] * len(group)

    try:
        # Read the file and set SNP as the index
        df2 = pd.read_csv(file_path, sep=r"\s+").set_index("SNP")
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return group.index, [np.nan] * len(group)

    # Map SNPs to their p-values
    p_values = group.iloc[:, 1].map(lambda snp: df2.loc[snp, df2.columns[-1]] if snp in df2.index else np.nan).tolist()
    return group.index, p_values

# Process groups in parallel and retain order
results = []
with ThreadPoolExecutor(max_workers=10) as executor:
    results = list(executor.map(lambda fg: process_field_group(fg[0], fg[1]), field_groups))

# Update the original dataframe with p-values
for group_index, p_values in results:
    df.loc[group_index, "p"] = p_values

# Save the updated dataframe to a new file
df.to_csv(output_file, index=False)

print(f"Updated file saved to {output_file}")
