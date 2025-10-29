#!/bin/bash
#SBATCH --job-name=BB
#SBATCH --array=1-670%200
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/testGxE_log/Aout/testgxe_BB/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/testGxE_log/Aerr/testgxe_BB/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan
#SBATCH --exclude=mulan-mc[01-02]


k=0
npart=10
cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/pheno/

tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67_shortName.txt | while read -r line;
do
    trait=$(echo "$line" | awk '{print $1}')
    
   for ipart in $(seq 1 $npart); do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			fastgxe --test-gxe --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
				--bfile /net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc \
				--data $trait.txt --trait $trait \
				--covar Genotype_batch:Sex_Age3 \
				--env-int Age:Confide --out ../testGxE/$trait --split-task $npart $ipart
			fi
	done
	
done
