import sys
import pandas as pd
sys.path.append('/net/zootopia/disk1/chaon/WORK/structLMM/')
from cal_LD import cal_LD

bed_file = "/net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc"


# Load the uploaded file
file_path = '/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/FTO_leadingSNP/trait_snp.txt'
data = pd.read_csv(file_path, sep='\t', header=None, names=['SNP', 'Trait', 'Abbreviation'])

# Group by SNP and aggregate the third column by joining values with ';'
grouped_data = data.groupby('SNP').agg({
    'Trait': 'first',  # Keep the first Trait value for each SNP
    'Abbreviation': lambda x: ';'.join(x)  # Join the Abbreviation values with ';'
}).reset_index()

concatenated_data = grouped_data.apply(lambda x: f"{x['SNP']}:{x['Abbreviation']}", axis=1)

snp_lst = list(grouped_data.iloc[:, 0])

snp_info, corr = cal_LD(bed_file, snp_lst)
corr = corr * corr
corr.columns = list(concatenated_data)
corr.index = list(concatenated_data)
corr.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/FTO_leadingSNP/corr.csv")
