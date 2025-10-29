#!/bin/bash
#SBATCH --job-name=mom
#SBATCH --array=1-36%36
#SBATCH --cpus-per-task=20
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Mental_health/Aout/data/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Mental_health/Aerr/data/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=main


k=0

tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/description/trait_Description.txt | while read -r line;
do
    trait=$(echo "$line" | awk '{print $1}')
    
   
   
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Mental_health/6.3_data.py $trait
		fi
	
done

