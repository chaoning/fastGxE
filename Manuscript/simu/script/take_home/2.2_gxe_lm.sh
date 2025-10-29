#!/bin/bash
#SBATCH --job-name=lm
#SBATCH --array=1-567%100
#SBATCH --cpus-per-task=5
#SBATCH --time=05:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/take_home/Aout/gxe_lm/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/take_home/Aerr/gxe_lm/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan,main
#SBATCH --exclude=mulan-mc[01-02]



k=0

for num_envi_power in 1 2 5 10 20 30 40; do
	for power_snp_h2_gxe in 0.005 0.01 0.02 0.04 0.06 0.08 0.1 0.12 0.14; do
		for sample_size in 50000 100000 200000 300000 400000 800000 1200000 1600000 2000000; do
			let k=${k}+1
			if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
				
				python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/take_home/2.2_gxe_lm.py \
				$num_envi_power \
				$power_snp_h2_gxe \
				$sample_size
			
			fi
		done
	done
done

