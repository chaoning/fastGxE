import sys
import argparse
import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed
from tqdm import tqdm


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process GRL and SNP data.")
    parser.add_argument("GRL", help="Gene-related list file path")
    parser.add_argument("pcol", help="Column name for p-values in SNP file")
    parser.add_argument("out_file", help="output file")
    parser.add_argument("sig_file", help="sig file")
    parser.add_argument("num_traits", help="# of traits")
    return parser.parse_args()


def main():
    args = parse_arguments()

    snp_on_disk = Bed("/net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc", count_A1=False)
    dfBim = pd.read_csv("/net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc.bim", header=None, sep=r"\s+")
    dct = dict(zip(dfBim[1], range(dfBim.shape[0])))

    dfGRL = pd.read_csv(args.GRL)
    dfSig = pd.read_csv(args.sig_file, sep=r"\s+")
    print("The number of traits:", args.num_traits)
    dfSig = dfSig[dfSig[args.pcol] < 5.0e-8 / int(args.num_traits)].copy()

    dfGRLTrait_lst = []

    for index, row in dfGRL.iterrows():
        print(index)
        dfGRLi = dfGRL.iloc[[index], :].copy()
        chrom, start, end = row['chrom'], row['start'], row['end']
        dfSub = dfSig[(dfSig["chrom"] == chrom) & (dfSig["base"] >= start) & (dfSig["base"] <= end)].copy()

        for trait, group in dfSub.groupby('trait', sort=True):
            group_sorted = group.sort_values(args.pcol)
            snp_lst = group_sorted["SNP"].to_numpy()
            snp_index_lst = [dct[snp] for snp in snp_lst]
            snpdata = snp_on_disk[:, snp_index_lst].read().val
            r2 = np.square(pd.DataFrame(snpdata).corr())

            leading_snp_indices = []
            while not r2.empty:
                leading_index = r2.columns[0]
                leading_snp_indices.append(leading_index)
                corr_arr = r2.iloc[0, 1:].to_numpy()
                r2 = r2.iloc[1:, 1:]
                r2 = r2.loc[corr_arr < 0.1, corr_arr < 0.1]

            leading_snps = snp_lst[leading_snp_indices]
            dfGRLi["trait"] = trait
            dfGRLi["trait_leading_snp"] = ";".join(leading_snps)
            dfGRLTrait_lst.append(dfGRLi.copy())

    dfGRLTrait = pd.concat(dfGRLTrait_lst)
    dfGRLTrait.to_csv(args.out_file, index=False)


if __name__ == "__main__":
    main()
