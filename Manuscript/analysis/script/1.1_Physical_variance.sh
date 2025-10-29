#!/bin/bash
#SBATCH --job-name=fastGxE
#SBATCH --array=1-32%32
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Aout/fastgxe_physical.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Aerr/fastgxe_physical.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan
#SBATCH --exclude=mulan-mc[01-06]


k=0

tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt | while read -r line;
do
    trait=$(echo "$line" | awk '{print $1}')
    
   
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
	cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/variance/
	time fastgxe --test-gxe --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
    	--data /net/zootopia/disk1/chaon/WORK/GxE/pheno/Physical_measures/IINT/pheno.e42.$trait.txt --trait $trait \
		--env-int Age:Confide --out $trait
	
	time fastgxe --test-gxe --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
    	--data /net/zootopia/disk1/chaon/WORK/GxE/pheno/Physical_measures/IINT/pheno.e42.$trait.txt --trait $trait \
		--env-int Age:Confide --out $trait\_noNxE --no-noisebye
	
	fi
	
done
