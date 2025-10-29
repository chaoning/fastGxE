#!/bin/bash
#SBATCH --job-name=clump
#SBATCH --array=1-99%99
#SBATCH --cpus-per-task=15
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/covariantGxE/Aout/clump/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/covariantGxE/Aerr/clump/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


traits=()
while read -r line; do
    trait=$(echo "$line" | awk '{print $1}')
    traits+=("$trait")
done < <(tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt)

while read -r line; do
    trait=$(echo "$line" | awk '{print $1}')
    traits+=("$trait")
done < <(tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67_shortName.txt)

cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/plink/

for trait in "${traits[@]}"; do
    echo "$trait"
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		plink \
		  --bfile /net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc.100000 \
		  --clump $trait.gxe \
		  --clump-p1 1e-5 \
		  --clump-r2 0.01 \
		  --clump-kb 250 \
		  --out $trait.clumped.r20.01
	fi
done
