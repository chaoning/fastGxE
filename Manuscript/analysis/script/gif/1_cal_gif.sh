#!/bin/bash
#SBATCH --job-name=gif
#SBATCH --array=1-99%99
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/gif/Aout/gif/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/gif/Aerr/gif/%a.txt
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
done < <(tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67_shortName.txt)

k=0

for trait in "${traits_pt[@]}"; do
    echo "$trait"
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/gif/1_cal_gif.py \
			/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/testGxE/$trait \
			/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/gif/traits/$trait.csv
	fi
done

for trait in "${traits_bb[@]}"; do
    echo "$trait"
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/gif/1_cal_gif.py \
			/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/testGxE/$trait \
			/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/gif/traits/$trait.csv
	fi
done

