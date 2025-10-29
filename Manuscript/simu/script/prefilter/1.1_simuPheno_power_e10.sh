#!/bin/bash
#SBATCH --job-name=filter
#SBATCH --array=1-2500%2500
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/prefilter/Aout/simu_power/e2.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/prefilter/Aerr/simu_power/e2.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan,main


cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno_filter/

power_snp_index=-1
h2_add=30
h2_gxe=5
h2_nxe=15
sample_size=100000
num_envi_power=10

k=0
for rGI in -0.9 -0.3 0 0.3 0.9; do
	for add_gxe_ratio in 0.005 0.25 1 2 4; do
		for rep in {1..100}; do
			let k=${k}+1
			if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
				
				python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/prefilter/1.1_simuPheno_power.py \
				$power_snp_index \
				$add_gxe_ratio \
				$num_envi_power \
				$h2_add \
				$h2_gxe \
				$h2_nxe \
				$sample_size \
				$rep \
				$rGI
				
			fi
		done
	done
done
