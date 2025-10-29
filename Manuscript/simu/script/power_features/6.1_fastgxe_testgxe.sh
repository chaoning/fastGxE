#!/bin/bash
#SBATCH --job-name=fastgxe
#SBATCH --array=1-100%100
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/Aout/fastgxe_testgxe/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/Aerr/fastgxe_testgxe/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_baselineLD_fastgxe/
k=0
for rep in {1..100}; do
    let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		
		prefix=h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
			time fastgxe --test-gxe  \
				--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
                --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
                --data $prefix.txt --trait trait$rep \
                --env-int age:alcohol_frequency \
                --threads 10 \
                --out ../res_baselineLD_fastgxe_testgxe/$prefix.rep_$rep
	
    fi
done
