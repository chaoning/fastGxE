import pandas as pd
import sys



def convert_to_plink_assoc(input_file_prefix, output_file):
    # Read the input file
    df_lst = []
    for i in range(10):
        input_file = f"{input_file_prefix}.10_{i+1}.res"
        df = pd.read_csv(input_file, sep=r"\s+", usecols=["SNP", "p_gxe"])
        df_lst.append(df)
    # Concatenate all DataFrames
    combined_df = pd.concat(df_lst, ignore_index=True)
    # Rename columns to match PLINK association format
    combined_df.rename(columns={"p_gxe": "P"}, inplace=True)
    combined_df = combined_df[["SNP", "P"]]
    # Save to output file
    combined_df.to_csv(output_file, sep='\t', index=False)

input_file_prefix = sys.argv[1]
output_file = sys.argv[2]

convert_to_plink_assoc(input_file_prefix, output_file)
