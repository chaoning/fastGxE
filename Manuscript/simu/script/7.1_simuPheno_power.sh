#!/bin/bash
#SBATCH --job-name=simu_power
#SBATCH --array=1-1400%1400
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/simu_power/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/simu_power/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno/

power_snp_index=-1
h2_add=30
h2_gxe=5
h2_nxe=15
sample_size=100000


k=0
for rep in {1..100}; do
	num_envi_power=30
	
	for power_snp_h2_gxe in 0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.25; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/7.1_simuPheno_power.py \
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

for rep in {1..100}; do
	power_snp_h2_gxe=0.16
	
	for num_envi_power in 1 2 5 10 20 40; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/7.1_simuPheno_power.py \
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
