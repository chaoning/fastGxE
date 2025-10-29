#!/bin/bash
#SBATCH --job-name=PT
#SBATCH --array=1-32%32
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/testGxE_log/Aout/vs/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/testGxE_log/Aerr/vs/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


k=0


tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt | while IFS=$'\t' read -r col1 col2 col3 rest
do
    trait="$col1"
    shortname="$col3"
    echo "trait: $trait, shortname: $shortname"
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/testGxE_log/3_trait_vs.py "$trait" "$shortname"
	fi

done
