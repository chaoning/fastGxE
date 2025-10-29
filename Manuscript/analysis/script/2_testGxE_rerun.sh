#!/bin/bash
#SBATCH --job-name=fastGxE
#SBATCH --array=[241,242,243,248,249,253,254,259,260,261,265,266,267,271,272,274,275,277,278,279,280,282,283,284,285,286,288,289,294,296,297,299,301,305,306,307,308,309,310,311,312,313,318,319,320]
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Aout/testgxe_physical.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Aerr/testgxe_physical.%a.txt
#SBATCH --mem-per-cpu=6000MB
#SBATCH --partition=mulan



k=0
npart=10
cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/testGxE/

tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt | while read -r line;
do
    trait=$(echo "$line" | awk '{print $1}')
    
   for ipart in $(seq 1 $npart); do
   
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		
		
		time fastgxe --test-gxe --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
			--bfile /net/zootopia/disk1/chaon/data/UKB/imp/ukb_imp_info.qc \
			--data /net/zootopia/disk1/chaon/WORK/GxE/pheno/Physical_measures/IINT/pheno.e42.$trait.txt --trait $trait \
			--env-int Age:Confide --out $trait --split-task $npart $ipart
		fi
	done
	
done
