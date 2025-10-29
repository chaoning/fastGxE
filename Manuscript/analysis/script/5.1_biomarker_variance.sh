#!/bin/bash
#SBATCH --job-name=fastGxE
#SBATCH --array=1-67%67
#SBATCH --cpus-per-task=5
#SBATCH --time=7-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Aout/fastgxe_biomarker.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Aerr/fastgxe_biomarker.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


k=0

cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/variance/

tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/trait67_shortName.txt | while read -r line;
do
    trait=$(echo "$line" | awk '{print $1}')
    
   
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
	
	time fastgxe --test-gxe --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
    	--data /net/zootopia/disk1/chaon/WORK/GxE/pheno/Biological_samples/IINT/pheno.e42.$trait.txt --trait $trait \
		--env-int Age:Confide --out $trait
	
	time fastgxe --test-gxe --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
    	--data /net/zootopia/disk1/chaon/WORK/GxE/pheno/Biological_samples/IINT/pheno.e42.$trait.txt --trait $trait \
		--env-int Age:Confide --out $trait\_noNxE --no-noisebye
	
	fi
	
done
