#!/bin/bash
#SBATCH --job-name=fastGxE
#SBATCH --array=1-360%120
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Mental_health/Aout/testgxe/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Mental_health/Aerr/testgxe/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan
#SBATCH --exclude=mulan-mc[01-03]


k=0
npart=10
cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/testGxE2/

tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/description/trait_Description.txt | while read -r line;
do
    trait=$(echo "$line" | awk '{print $1}')
    for ipart in $(seq 1 $npart); do
	
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		
		
			time fastgxe --test-gxe --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
				--bfile /net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc \
				--covar Genotype_batch:Sex_Age3 \
				--data /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/pheno/$trait.txt --trait $trait \
				--env-int Age:Confide --out $trait --split-task $npart $ipart
		
		fi
	done
	
done
