#!/bin/bash
#SBATCH --job-name=fastgxe
#SBATCH --array=1-700%700
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/fastgxe_testgxe_baselineLD/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/fastgxe_testgxe_baselineLD/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan
#SBATCH --exclude=mulan-mc[01-03]


cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno_baselineLD_fastgxe/
h2_add=30
h2_gxe=5
h2_nxe=15
sample_size=100000


k=0
for rep in {1..100}; do
	
	num_envi_power=30
	
	for power_snp_h2_gxe in 0.02 0.04 0.06 0.08 0.1 0.12 0.14; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			prefix=power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.h2_add_$h2_add.h2_gxe_$h2_gxe.h2_nxe_$h2_nxe.sample_$sample_size
			bim_file=../pheno_baselineLD/$prefix.rep_$rep.power_snps.bim
			snp_start=$(awk 'NR==1 {print $2}' $bim_file)
			snp_end=$(awk 'END {print $2}' $bim_file)
			time fastgxe --test-gxe  \
				--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
                --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
                --data $prefix.txt --trait trait$rep \
                --env-int age:alcohol_frequency \
                --threads 10 \
                --out ../res_baselineLD_fastgxe/$prefix.rep_$rep \
				--snp-range $snp_start $snp_end
		
		fi
	done

done