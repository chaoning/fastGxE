#!/bin/bash
#SBATCH --job-name=PT
#SBATCH --array=1-99%99
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/testGxE_log/Aout/log/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/testGxE_log/Aerr/log/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


k=0

while read -r line; do
		trait=$(echo "$line" | awk '{print $1}')
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/testGxE_log/1_logTrait.py \
			/net/zootopia/disk1/chaon/WORK/GxE/pheno/Physical_measures/pheno.e42.$trait.txt \
			/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/pheno/$trait.txt
		fi
done < <(awk 'NR>1' /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt)


while read -r line; do
		trait=$(echo "$line" | awk '{print $1}')
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/testGxE_log/1_logTrait.py \
			/net/zootopia/disk1/chaon/WORK/GxE/pheno/Biological_samples/pheno.e42.$trait.txt \
			/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/pheno/$trait.txt
		fi
done < <(awk 'NR>1' /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67_shortName.txt)
