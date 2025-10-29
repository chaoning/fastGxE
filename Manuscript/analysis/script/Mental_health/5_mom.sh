#!/bin/bash
#SBATCH --job-name=mom
#SBATCH --array=1-36%36
#SBATCH --cpus-per-task=20
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Mental_health/Aout/mom/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Mental_health/Aerr/mom/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


k=0
cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/mom/

tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/description/trait_Description.txt | while read -r line;
do
    trait=$(echo "$line" | awk '{print $1}')
    
   
   
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		
		
		time fastgxe --mom \
			--bfile /net/zootopia/disk1/chaon/data/UKB/chip/ukb.keep.qc \
			--data /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/INT/$trait.txt --trait $trait \
			--threads 20 \
			--env-int Age:Confide --out $trait
		fi
	
done

