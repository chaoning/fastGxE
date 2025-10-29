#!/bin/bash
#SBATCH --job-name=fastgxe
#SBATCH --array=1-3000%100
#SBATCH --cpus-per-task=5
#SBATCH --time=1:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/addToInt/Aout/fastgxe_testgxe2/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/addToInt/Aerr/fastgxe_testgxe2/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_addToint_fastgxe/
k=0
for rep in {1..100}; do
	for unknown_h2SNP_h2E in 0.02 0.08 0.16; do
		for R2 in 0 0.01 0.1 0.3 0.5 0.7 0.9; do
			for R2GE in 0; do
				let k=${k}+1
				if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
					
					prefix=h2_$unknown_h2SNP_h2E.R2_$R2.R2GE_$R2GE.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
					fastgxe --test-gxe  \
							--bfile ../pheno_addToint/$prefix.rep_$rep.power_snp \
							--grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
							--data $prefix.txt --trait trait$rep \
							--env-int age:alcohol_frequency \
							--threads 10 --num-random-snp 1 \
							--out ../res_addToint_fastgxe/$prefix.rep_$rep
				
				fi
			done
		done
		
		for R2 in 0.5; do
			for R2GE in 0.01 0.1 0.3; do
				let k=${k}+1
				if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
					
					prefix=h2_$unknown_h2SNP_h2E.R2_$R2.R2GE_$R2GE.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
					fastgxe --test-gxe  \
							--bfile ../pheno_addToint/$prefix.rep_$rep.power_snp \
							--grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
							--data $prefix.txt --trait trait$rep \
							--env-int age:alcohol_frequency \
							--threads 10 --num-random-snp 1 \
							--out ../res_addToint_fastgxe/$prefix.rep_$rep
				
				fi
			done
		done
		
	done
done
