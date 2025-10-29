import os
import sys
import re

import numpy as np
import pandas as pd
from tqdm import tqdm

maxDist = float(sys.argv[1])
coding = int(sys.argv[2])

in_file = sys.argv[3]
out_file = sys.argv[4]

dfM = pd.read_csv(in_file, dtype={"chrom": "str"})
dfG = pd.read_csv("/net/zootopia/disk1/chaon/data/human_ref/gtf/gencode.v45.basic.annotation.gtf", sep="\t",
                  skiprows=5, header=None)
chr_lst = []
for chr in dfG.iloc[:, 0]:
    chr2 = chr.replace("chr", "")
    chr_lst.append(chr2)

dfG[0] = chr_lst


def extract_gene(data):
    # 提取 gene_type
    gene_type_match = re.search(r'gene_type "([^"]+)"', data)
    gene_type = gene_type_match.group(1) if gene_type_match else None
    # 提取 gene_name
    gene_name_match = re.search(r'gene_name "([^"]+)"', data)
    gene_name = gene_name_match.group(1) if gene_name_match else None
    dct_gene_type = {gene_name: gene_type}
    return gene_name, dct_gene_type

gene_num_lst = []
gene_info_lst = []
for i in tqdm(range(dfM.shape[0])):
    chr = dfM.iloc[i, 1]
    start = dfM.iloc[i, -3]
    end = dfM.iloc[i, -2]
    a = (dfG[0] == chr)
    b = (dfG[4] > start)
    c = (dfG[3] < end)
    dfGsub = dfG[a & b & c].copy()
    dct_gene_type = {}
    dct_gene_dis = {}
    for k in range(dfGsub.shape[0]):
        sinfo = dfGsub.iloc[k, 8]
        gene_name, dct = extract_gene(sinfo)
        dct_gene_type.update(dct)
        dct_gene_dis[gene_name] = [-1]

    dfGsub = dfG[a].copy()
    print(dfGsub.shape[0])
    dis_arr = np.zeros((dfGsub.shape[0], 2))
    dis_arr[:, 0] = np.absolute(dfGsub[3] - end)
    dis_arr[:, 1] = np.absolute(dfGsub[4] - start)
    dfGsub["dis"] = np.min(dis_arr, axis=1)
    dfGsub = dfGsub[dfGsub["dis"] < maxDist]
    for k in range(dfGsub.shape[0]):
        sinfo = dfGsub.iloc[k, 8]
        gene_name, dct = extract_gene(sinfo)
        dct_gene_type.update(dct)
        try:
            dct_gene_dis[gene_name].append(dfGsub.iloc[k, 9])
        except Exception as e:
            dct_gene_dis[gene_name] = [dfGsub.iloc[k, 9]]
    stri = ""
    for gene in dct_gene_type.keys():
        if coding == 1:
            if dct_gene_type[gene] == "protein_coding":
                stri = stri + f"{gene}({dct_gene_type[gene]}, {np.min(dct_gene_dis[gene])}); "
        else:
            stri = stri + f"{gene}({dct_gene_type[gene]}, {np.min(dct_gene_dis[gene])}); "
    gene_info_lst.append(stri)
    gene_num_lst.append(len(dct_gene_type))

dfM["nGenes"] = gene_num_lst
dfM["Genes"] = gene_info_lst

dfM.to_csv(out_file, index=False)
