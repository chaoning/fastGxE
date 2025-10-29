#!/bin/bash
#SBATCH --job-name=fastgxe
#SBATCH --array=1-2500%2500
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/prefilter/Aout/fastgxe_testgxe/e10.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/prefilter/Aerr/fastgxe_testgxe/e10.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan
#SBATCH --exclude=mulan-mc[01-02]

num_envi=10
cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno_filter_fastgxe/
k=0
for rGI in -0.9 -0.3 0 0.3 0.9; do
	for add_gxe_ratio in 0.005 0.25 1 2 4; do
		for rep in {1..100}; do
			let k=${k}+1
			if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
				
					prefix=rGI_$rGI.add_gxe_ratio_$add_gxe_ratio.num_envi_$num_envi.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
					time fastgxe --test-gxe  \
						--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
						--grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
						--data $prefix.txt --trait trait$rep \
						--env-int age:alcohol_frequency \
						--threads 10 \
						--out ../res_filter_fastgxe/$prefix.rep_$rep --num-random-snp 1000
			
			fi
		done
	done
done
