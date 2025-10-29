#!/bin/bash
#SBATCH --job-name=to_plink
#SBATCH --array=1-67%67
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/covariantGxE/Aout/to_plink/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/covariantGxE/Aerr/to_plink/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan



k=0


tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67_shortName.txt | while read -r line;
do
    trait=$(echo "$line" | awk '{print $1}')
    
   
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/covariantGxE/1.1_to_plink_assoc.py \
		/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/testGxE/$trait \
		/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/plink/$trait.gxe
	fi

	
done
