#!/bin/bash
#SBATCH --job-name=fastgxe
#SBATCH --array=1-2000%120
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/Aout/power/fastgxe_testgxe/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/Aerr/power/fastgxe_testgxe/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno_feature_fastgxe/
h2_add=30
h2_gxe=5
h2_nxe=15
sample_size=100000
num_envi_power=30
power_snp_h2_gxe=0.08

k=0
for rep in {1..100}; do
	
	
	for nBin in {1..6}; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			prefix=maf_$nBin.power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.h2_add_$h2_add.h2_gxe_$h2_gxe.h2_nxe_$h2_nxe.sample_$sample_size
			bim_file=../pheno_feature/$prefix.rep_$rep.power_snps.bim
			snp_start=$(awk 'NR==1 {print $2}' $bim_file)
			snp_end=$(awk 'END {print $2}' $bim_file)
			time fastgxe --test-gxe  \
				--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
                --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
                --data $prefix.txt --trait trait$rep \
                --env-int age:alcohol_frequency \
                --threads 10 \
                --out ../res_fastgxe_testgxe_feature/$prefix.rep_$rep \
				--snp-range $snp_start $snp_end
			
			time fastgxe --test-gxe  \
			--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
                --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
                --data $prefix.txt --trait trait$rep \
                --env-int age:alcohol_frequency \
                --threads 10 \
                --out ../res_fastgxe_testgxe_feature/$prefix.rep_$rep.noNxE --no-noisebye \
				--snp-range $snp_start $snp_end
		
		fi
	done
	
	
	arr=()
	while read -r col1 col2; do
		arr+=("$col1")
	done < /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQplot_vf/variant_function.txt

	for val in "${arr[@]}"; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			
			prefix=variant_function.$val.power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.h2_add_$h2_add.h2_gxe_$h2_gxe.h2_nxe_$h2_nxe.sample_$sample_size
			bim_file=../pheno_feature/$prefix.rep_$rep.power_snps.bim
			snp_start=$(awk 'NR==1 {print $2}' $bim_file)
			snp_end=$(awk 'END {print $2}' $bim_file)
			time fastgxe --test-gxe  \
				--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
                --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
                --data $prefix.txt --trait trait$rep \
                --env-int age:alcohol_frequency \
                --threads 10 \
                --out ../res_fastgxe_testgxe_feature/$prefix.rep_$rep \
				--snp-range $snp_start $snp_end
			
			time fastgxe --test-gxe  \
			--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
                --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
                --data $prefix.txt --trait trait$rep \
                --env-int age:alcohol_frequency \
                --threads 10 \
                --out ../res_fastgxe_testgxe_feature/$prefix.rep_$rep.noNxE --no-noisebye \
				--snp-range $snp_start $snp_end
		
		fi
	done

	
	for nBin in {1..5}; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		
			prefix=ld_$nBin.power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.h2_add_$h2_add.h2_gxe_$h2_gxe.h2_nxe_$h2_nxe.sample_$sample_size
			bim_file=../pheno_feature/$prefix.rep_$rep.power_snps.bim
			snp_start=$(awk 'NR==1 {print $2}' $bim_file)
			snp_end=$(awk 'END {print $2}' $bim_file)
			time fastgxe --test-gxe  \
				--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
                --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
                --data $prefix.txt --trait trait$rep \
                --env-int age:alcohol_frequency \
                --threads 10 \
                --out ../res_fastgxe_testgxe_feature/$prefix.rep_$rep \
				--snp-range $snp_start $snp_end
			
			time fastgxe --test-gxe  \
			--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
                --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
                --data $prefix.txt --trait trait$rep \
                --env-int age:alcohol_frequency \
                --threads 10 \
                --out ../res_fastgxe_testgxe_feature/$prefix.rep_$rep.noNxE --no-noisebye \
				--snp-range $snp_start $snp_end

		fi
	done
	
done
