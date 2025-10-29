#!/bin/bash
#SBATCH --job-name=struct
#SBATCH --array=1-600%100
#SBATCH --cpus-per-task=5
#SBATCH --time=7-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_mask/Aout/struct/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_mask/Aerr/struct/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan,main
#SBATCH --exclude=mulan-mc[01-03]


bash

let k=0

sample_size=100000
cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_mask_structlmm/
for rmE in 1 2 5 10 20 30;do
	for rep in {1..100};do
		
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			prefix=power_snp_h2_gxe_0.08.num_envi_30.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_mask/2.2_structlmm.py \
				$prefix $sample_size $rmE $rep
		fi
	done
done
