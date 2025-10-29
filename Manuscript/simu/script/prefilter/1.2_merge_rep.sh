#!/bin/bash
#SBATCH --job-name=simu
#SBATCH --array=1-30%30
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/simu_null.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/simu_null.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan,main

k=0
num_envi=2
for rGI in -0.9 -0.3 0 0.3 0.9; do
	for add_gxe_ratio in 0.006 0.02 0.25 1 2 4; do
			let k=${k}+1
			if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
				
				python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/prefilter/1.2_merge_rep.py \
				rGI_$rGI.add_gxe_ratio_$add_gxe_ratio.num_envi_$num_envi.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
				
			fi
	done
done
