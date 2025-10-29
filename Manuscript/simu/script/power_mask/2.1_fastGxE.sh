#!/bin/bash
#SBATCH --job-name=fastGxE
#SBATCH --array=1-600%100
#SBATCH --cpus-per-task=5
#SBATCH --time=7-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_mask/Aout/fastGxE/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_mask/Aerr/fastGxE/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan
#SBATCH --exclude=mulan-mc[01-03]


bash

let k=0


cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_mask_fastGxE/
for rmE in 1 2 5 10 20 30;do
	for rep in {1..100};do
		
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			envi_arr=($(awk '{print $2}' envi.rm$rmE.$rep.txt))
			prefix=power_snp_h2_gxe_0.08.num_envi_30.h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
			bim_file=../pheno/$prefix.rep_$rep.power_snp.bim
			snp_name=$(awk 'NR==1 {print $2}' $bim_file)
			fastgxe --test-gxe  \
				--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
                --grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
                --data ../pheno_fastgxe/$prefix.txt --trait trait$rep \
                --env-int ${envi_arr[@]} \
                --threads 5 \
                --out rmE_$rmE.rep_$rep \
				--snp-range $snp_name $snp_name
		fi
	done
done
