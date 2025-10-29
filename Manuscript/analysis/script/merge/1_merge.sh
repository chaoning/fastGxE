#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --array=1-99%99
#SBATCH --cpus-per-task=5
#SBATCH --time=7-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/merge/Aout/merge/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/merge/Aerr/merge/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


traits_pt=()
while read -r line; do
    trait=$(echo "$line" | awk '{print $1}')
    traits_pt+=("$trait")
done < <(tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt)

traits_bb=()
while read -r line; do
    trait=$(echo "$line" | awk '{print $1}')
    traits_bb+=("$trait")
done < <(tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67.txt)



k=0

for trait in "${traits_pt[@]}"; do
    echo "$trait"
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/merge/1_merge.py \
			$trait \
			/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/testGxE/
	fi
done

for trait in "${traits_bb[@]}"; do
    echo "$trait"
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/merge/1_merge.py \
			$trait \
			/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/testGxE/
	fi
done
