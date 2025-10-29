#!/bin/bash
#SBATCH --job-name=validate
#SBATCH --array=1-77%77
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/AIE_split/Aout/validate/BB.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/AIE_split/Aerr/validate/BB.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


k=0

tail -n +2 "/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP.csv" | while IFS= read -r line; do
    # Using space as a delimiter to split the line into an array
    IFS=',' read -ra parts <<< "$line"
    trait=${parts[13]}  # Adjusted to 0-based indexing
    snp=${parts[14]}    # Adjusted to 0-based indexing
    
    let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
    echo $trait $snp
        python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/AIE_split/1.1_mmsusie.py biomarker $trait $snp
    fi
	
done
