#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --array=1-33%33
#SBATCH --cpus-per-task=5
#SBATCH --time=1:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/addToInt/Aout/merge/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/addToInt/Aerr/merge/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan

k=0

for unknown_h2SNP_h2E in 0.02 0.08 0.16; do
	for R2 in 0 0.01 0.1 0.3 0.5 0.7 0.9; do
		for R2GE in 0; do
			let k=${k}+1
			if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
				python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/addToInt/1.2_merge_rep.py \
				h2_$unknown_h2SNP_h2E.R2_$R2.R2GE_$R2GE.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
			fi
		done
	done
	
	for R2 in 0.5; do
		for R2GE in 0.01 0.1 0.3; do
			let k=${k}+1
			if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
				python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/addToInt/1.2_merge_rep.py \
				h2_$unknown_h2SNP_h2E.R2_$R2.R2GE_$R2GE.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
			fi
		done
	done
done
