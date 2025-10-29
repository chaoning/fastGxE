#!/bin/bash
#SBATCH --job-name=susie
#SBATCH --array=1-77%77
#SBATCH --cpus-per-task=5
#SBATCH --time=10-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/mmsusie/Aout/mmsusie/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/mmsusie/Aerr/mmsusie/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan

let k=0

tail -n +2 "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP_logP.csv" | while IFS= read -r line; do
    # Using space as a delimiter to split the line into an array
    IFS=',' read -ra parts <<< "$line"
    trait=${parts[13]}  # Adjusted to 0-based indexing
    snp=${parts[14]}    # Adjusted to 0-based indexing
    
	echo $trait $snp
    let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
        python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/mmsusie/mmSuSiE.py $trait $snp \
		/net/zootopia/disk1/chaon/WORK/GxE/pheno/Biological_samples/IINT/pheno.e42.$trait.txt \
		/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/testGxE/$trait.var 
    fi
done
