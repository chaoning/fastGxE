#!/bin/bash
#SBATCH --job-name=validate
#SBATCH --array=1-2760%1000
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/null/validate/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/null/validate/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


k=0
for rep in {1..10}; do
	for index in {1..138}; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/10.1_validate_mmsusie.py $rep $index
			
		fi
		
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/10.1_validate_mmsusie_random.py $rep $index
			
		fi
		
	done
done
