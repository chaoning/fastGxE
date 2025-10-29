#!/bin/bash
#SBATCH --job-name=fastGxE
#SBATCH --array=[205,206,207,250,421,456,505,533,607,613]
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Aout/testgxe_biomarker.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Aerr/testgxe_biomarker.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan



k=0
npart=10
cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/testGxE/

tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67_shortName.txt | while read -r line;
do
    trait=$(echo "$line" | awk '{print $1}')
    
   for ipart in $(seq 1 $npart); do
   
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		
		
		time fastgxe --test-gxe --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
			--bfile /net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc \
			--data /net/zootopia/disk1/chaon/WORK/GxE/pheno/Biological_samples/IINT/pheno.e42.$trait.txt --trait $trait \
			--env-int Age:Confide --out $trait --split-task $npart $ipart
		fi
	done
	
done
