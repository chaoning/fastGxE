#!/bin/bash
#SBATCH --job-name=fastgxe
#SBATCH --array=1-400%100
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_fastgxe_testgxe_nonactiE/Aout/power/fastgxe_testgxe/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_fastgxe_testgxe_nonactiE/Aerr/power/fastgxe_testgxe/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno/

prefix=h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
num_envi_power=10
power_snp_h2_gxe=0.08

k=0
for rep in {1..100}; do
	
	for nonE in 0 10 20 30;do
		arr=()
		while read -r col1; do
			arr+=("$col1")
		done < ../res_fastgxe_testgxe_nonactiE/power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.$prefix.rep_$rep.rndE_$nonE.txt

		
			let k=${k}+1
			if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
				
				bim_file=power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.$prefix.rep_$rep.power_snp.bim
				snp_start=$(awk 'NR==1 {print $2}' $bim_file)
				snp_end=$snp_start
				time fastgxe --test-gxe  \
					--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
					--grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
					--data ../pheno_fastgxe/power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.$prefix.txt --trait trait$rep \
					--env-int ${arr[@]} \
					--threads 10 \
					--out ../res_fastgxe_testgxe_nonactiE/nonE_$nonE.rep_$rep \
					--snp-range $snp_start $snp_end
			
			fi
			
	done
done
