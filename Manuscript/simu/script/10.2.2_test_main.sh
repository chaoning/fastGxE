#!/bin/bash
#SBATCH --job-name=validate
#SBATCH --array=1-6900%1000
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/null/test_main/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/null/test_main/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan



k=0
for rep in {1..10}; do
	for index in {1..138}; do
		for group in {1..5}; do
			let k=${k}+1
			if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			
			
				snp_name=$(cat /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/validate_mmsusie/$rep.$index.mmsusie_out.lead_snp)
				fastgxe --test-main  \
					--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
					--grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
					--covar age:alcohol_frequency \
					--env-int age:alcohol_frequency \
					--data /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/validate_AIE/$rep.$index.score_group_$group.txt \
					--trait trait \
					--threads 10 \
					--out /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/validate_AIE_test_main/$rep.$index.score_group_$group \
					--snp-range $snp_name $snp_name
				
				snp_name=$(cat /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/validate_mmsusie_random/$rep.$index.mmsusie_out.lead_snp)
				fastgxe --test-main  \
					--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
					--grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
					--covar age:alcohol_frequency \
					--env-int age:alcohol_frequency \
					--data /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/validate_AIE_random/$rep.$index.score_group_$group.txt \
					--trait trait \
					--threads 10 \
					--out /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/validate_AIE_test_main_random/$rep.$index.score_group_$group \
					--snp-range $snp_name $snp_name
				
				
			fi
		done
	done
done
