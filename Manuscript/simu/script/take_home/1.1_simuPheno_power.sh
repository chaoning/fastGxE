#!/bin/bash
#SBATCH --job-name=simu_power
#SBATCH --array=5993-24500%160
#SBATCH --cpus-per-task=5
#SBATCH --time=05:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/take_home/Aout/simu_power/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/take_home/Aerr/simu_power/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


power_snp_index=-1
h2_add=30
h2_gxe=5
h2_nxe=15

k=0
for rep in {1..100}; do
	for num_envi_power in 1 2 5 10 20 30 40; do
		for power_snp_h2_gxe in 0.02 0.04 0.06 0.08 0.1 0.12 0.14; do
			for sample_size in 50000 100000 200000 300000 400000; do
				let k=${k}+1
				if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
					
					python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/take_home/1.1_simuPheno_power.py \
					$power_snp_index \
					$power_snp_h2_gxe \
					$num_envi_power \
					$h2_add \
					$h2_gxe \
					$h2_nxe \
					$sample_size \
					$rep
				
				fi
			done
		done
	done
done
