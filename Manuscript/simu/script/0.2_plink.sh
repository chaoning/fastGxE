sbatch -p mulan -c 80 --wrap="plink --bfile /net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc --chr 22 --make-bed --out plink22"
sbatch -p mulan -c 80 --wrap="plink --bfile plink22 --freqx"
sbatch -p mulan -c 80 --wrap="plink --bfile plink22 --freq"

sbatch -p mulan -c 80 --wrap="plink --bfile plink22 --r2 --ld-window 999999 --ld-window-kb 1000 --ld-window-r2 0 --out plink_ld"
