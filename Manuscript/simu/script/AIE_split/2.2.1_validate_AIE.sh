#!/bin/bash
#SBATCH --job-name=validate
#SBATCH --array=1-13800%1000
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/AIE_split/Aout/validate_AIE/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/AIE_split/Aerr/validate_AIE/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


k=0
for rep in {1..100}; do
	for index in {1..138}; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/AIE_split/2.2.1_validate_AIE.py \
			/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/AIE_split_mmsusie/ \
			/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/AIE_split_group/ \
			$rep $index
			
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/AIE_split/2.2.1_validate_AIE.py \
			/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/AIE_split_mmsusie_random/ \
			/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/AIE_split_group_random/ \
			$rep $index
			
		fi
	done
done
