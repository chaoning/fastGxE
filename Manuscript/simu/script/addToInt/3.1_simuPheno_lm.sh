#!/bin/bash
#SBATCH --job-name=simu
#SBATCH --array=1-42%42
#SBATCH --cpus-per-task=10
#SBATCH --time=0-10:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/addToInt/Aout/lm/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/addToInt/Aerr/lm/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan

k=0
for unknown_h2SNP_h2E in 0.02 0.08 0.16 0.32 0.64 1; do
	
	for R2 in 0 0.01 0.1 0.3 0.5 0.7 0.9; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/addToInt/3.1_simuPheno_lm.py $unknown_h2SNP_h2E $R2
		fi
	done
	
done
