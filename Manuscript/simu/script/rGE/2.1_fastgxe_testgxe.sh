#!/bin/bash
#SBATCH --job-name=fastgxe
#SBATCH --array=1-500%120
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/rGE/Aout/fastgxe_testgxe/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/rGE/Aerr/fastgxe_testgxe/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_rGE/
h2_add=30
h2_gxe=5
h2_nxe=15
sample_size=100000
power_snp_h2_add=0.2

k=0
for rep in {1..100}; do
	
	
	for rGE in 0 0.001 0.01 0.1 0.5; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			prefix=power_snp_h2_add_$power_snp_h2_add.rGE_$rGE.h2_add_$h2_add.h2_gxe_$h2_gxe.h2_nxe_$h2_nxe.sample_$sample_size.rep_$rep
			bim_file=$prefix.power_snp.bim
			snp_name=$(awk 'NR==1 {print $2}' $bim_file)
			time fastgxe --test-gxe  \
				--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
                --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
                --data $prefix.txt --trait pheno \
                --env-int age:alcohol_frequency \
                --threads 10 \
                --out ../res_rGE_fastgxe/$prefix \
				--snp-range $snp_name $snp_name
		
		fi
	done
	

done