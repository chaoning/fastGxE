#!/bin/bash
#SBATCH --job-name=r
#SBATCH --array=1-99%99
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/covariantGxE/Aout/r/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/covariantGxE/Aerr/r/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


traits1=()
while read -r line; do
    trait=$(echo "$line" | awk '{print $1}')
    traits1+=("$trait")
done < <(tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt)

traits2=()
while read -r line; do
    trait=$(echo "$line" | awk '{print $1}')
    traits2+=("$trait")
done < <(tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67_shortName.txt)


for trait in "${traits1[@]}"; do
   
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		echo "$trait"
		python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/covariantGxE/1.3_r_r2_0.01.py \
			$trait \
			/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/testGxE/ \
			/net/zootopia/disk1/chaon/WORK/GxE/pheno/Physical_measures/IINT/
	fi
	
done

for trait in "${traits2[@]}"; do
   
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		echo "$trait"
		python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/covariantGxE/1.3_r_r2_0.01.py \
			$trait \
			/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/testGxE/ \
			/net/zootopia/disk1/chaon/WORK/GxE/pheno/Biological_samples/IINT/
	fi
	
done
