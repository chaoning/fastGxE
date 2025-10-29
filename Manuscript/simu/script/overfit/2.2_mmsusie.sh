#!/bin/bash
#SBATCH --job-name=mmsusie
#SBATCH --array=1-4900%100
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/overfit/Aout/mmsusie/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/overfit/Aerr/mmsusie/%a.txt
#SBATCH --mem-per-cpu=6000MB
#SBATCH --partition=mulan

h2_add=30
h2_gxe=5
h2_nxe=15
sample_size=100000

k=0
for rep in {1..100}; do
	for num_envi_power in 1 2 5 10 20 30 40; do
		for power_snp_h2_gxe in 0.02 0.04 0.06 0.08 0.1 0.12 0.14; do
				let k=${k}+1
				if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
					prefix=power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.h2_add_$h2_add.h2_gxe_$h2_gxe.h2_nxe_$h2_nxe.sample_$sample_size
					python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/overfit/2.2_mmsusie.py $prefix $rep
				
				fi
		done
	done
done
