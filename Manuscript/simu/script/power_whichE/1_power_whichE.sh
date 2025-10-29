#!/bin/bash
#SBATCH --job-name=fastgxe
#SBATCH --array=1-2205%2205
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_whichE/Aout/test/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_whichE/Aerr/test/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


k=0

for power_snp_h2_gxe in 0.005 0.01 0.02 0.04 0.06 0.08 0.1 0.12 0.14; do
	for sample_size in 50000 100000 200000 300000 400000; do
		for num_envi_power in 1 2 5 10 20 30 40; do
			for topE in 1 2 5 10 20 30 40;do
				
					let k=${k}+1
					if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
						
						python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_whichE/1_power_whichE.py \
							power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.h2_add_30.h2_gxe_5.h2_nxe_15.sample_$sample_size $topE
					
					fi
					
			done
		done
	done
done
