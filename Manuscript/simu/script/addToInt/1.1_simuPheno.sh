#!/bin/bash
#SBATCH --job-name=simu
#SBATCH --array=1-3300%300
#SBATCH --cpus-per-task=5
#SBATCH --time=0-10:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/addToInt/Aout/simu/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/addToInt/Aerr/simu/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan,main

k=0
for rep in {1..100}; do
	for unknown_h2SNP_h2E in 0.02 0.08 0.16; do
		
		for R2 in 0 0.01 0.1 0.3 0.5 0.7 0.9; do
			for R2GE in 0; do
				let k=${k}+1
				if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
					python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/addToInt/1.1_simuPheno.py \
					$unknown_h2SNP_h2E $R2 $R2GE 30 5 15 100000 $rep
				fi
			done
		done
		
		for R2 in 0.5; do
			for R2GE in 0 0.01 0.1 0.3; do
				let k=${k}+1
				if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
					python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/addToInt/1.1_simuPheno.py \
					$unknown_h2SNP_h2E $R2 $R2GE 30 5 15 100000 $rep
				fi
			done
		done
		
	done
done
