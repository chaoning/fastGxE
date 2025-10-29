#!/bin/bash
#SBATCH --job-name=sum
#SBATCH --array=1-30%30
#SBATCH --cpus-per-task=20
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/prefilter/Aout/sum/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/prefilter/Aerr/sum/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan,main

num_envi=2
k=0
for rGI in -0.9 -0.3 0 0.3 0.9; do
	for add_gxe_ratio in 0.006 0.02 0.25 1 2 4; do
			let k=${k}+1
			if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
				
					python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/prefilter/3.1_power_sum.py $rGI $add_gxe_ratio $num_envi
			
			fi
	done
done
