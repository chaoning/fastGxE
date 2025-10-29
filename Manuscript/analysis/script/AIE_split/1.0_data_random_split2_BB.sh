#!/bin/bash
#SBATCH --job-name=split
#SBATCH --array=1-67%67
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/AIE_split/Aout/split/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/AIE_split/Aerr/split/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan



k=0
npart=10
cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/AIE_split/pheno/

tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67.txt | while read -r line;
do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		trait=$(echo "$line" | awk '{print $1}')
		python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/AIE_split/1.0_data_random_split2.py \
		/net/zootopia/disk1/chaon/WORK/GxE/pheno/Biological_samples/IINT/pheno.e42.$trait.txt \
		$trait.part1.txt $trait.part2.txt
	fi
	
done
