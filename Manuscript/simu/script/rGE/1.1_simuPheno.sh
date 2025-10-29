#!/bin/bash
#SBATCH --job-name=simu
#SBATCH --array=1-500%500
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/rGE/Aout/simu/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/rGE/Aerr/simu/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_rGE/

power_snp_index=-1
h2_add=30
h2_gxe=5
h2_nxe=15
sample_size=100000
power_snp_h2_add=0.2

k=0
for rep in {1..100}; do
	
	for rGE in 0 0.001 0.01 0.1 0.5; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/rGE/1.1_simuPheno.py \
			$power_snp_index \
			$power_snp_h2_add \
			$rGE \
			$h2_add \
			$h2_gxe \
			$h2_nxe \
			$sample_size \
			$rep
		
		fi
	done
	
done
