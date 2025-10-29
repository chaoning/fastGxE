#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --array=1-10%10
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/simu_power/merge.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/simu_power/merge.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


h2_add=30
h2_gxe=5
h2_nxe=15
sample_size=100000


k=0

	num_envi_power=30
	
	for power_snp_h2_gxe in 0.02 0.04 0.06 0.08; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			prefix=power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.h2_add_$h2_add.h2_gxe_$h2_gxe.h2_nxe_$h2_nxe.sample_$sample_size
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/7.2_merge_rep.py $prefix
		
		fi
	done
	

	power_snp_h2_gxe=0.08
	
	for num_envi_power in 1 2 5 10 20 40; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			prefix=power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.h2_add_$h2_add.h2_gxe_$h2_gxe.h2_nxe_$h2_nxe.sample_$sample_size
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/7.2_merge_rep.py $prefix
		
		fi
	done

