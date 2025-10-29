#!/bin/bash
#SBATCH --job-name=test_main
#SBATCH --array=1-108%108
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/AIE_split/Aout/test_main/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/AIE_split/Aerr/test_main/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


traits=()
snps=()

#
while IFS= read -r line; do
    IFS=',' read -ra parts <<< "$line"
    trait=${parts[13]}
    snp=${parts[14]}
    traits+=("$trait")
    snps+=("$snp")
done < <(tail -n +2 "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP.csv")

#
while IFS= read -r line; do
    IFS=',' read -ra parts <<< "$line"
    trait=${parts[13]}
    snp=${parts[14]}
    traits+=("$trait")
    snps+=("$snp")
done < <(tail -n +2 "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP.csv")

for i in {0..107}; do
	let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		echo "${traits[$i]}  ${snps[$i]}"
		python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/AIE_split/3.2.2_test_main_py.py ${traits[$i]} ${snps[$i]}
	fi
done
