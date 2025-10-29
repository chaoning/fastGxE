#!/bin/bash
#SBATCH --job-name=topE
#SBATCH --array=1-700%700
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_whichE/Aout/topE/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_whichE/Aerr/topE/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


k=0


for rep in {1..100}; do
	for topE in 1 2 5 10 20 30 40;do
		
			let k=${k}+1
			if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
				
				python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_whichE/3_topE_null.py $rep $topE
			
			fi
			
	done
done

