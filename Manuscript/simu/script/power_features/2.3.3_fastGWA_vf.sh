#!/bin/bash
#SBATCH --job-name=fastGWA
#SBATCH --array=1-36000%300
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/Aout/power/fastGWA_vf/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/Aerr/power/fastGWA_vf/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan,main
#SBATCH --exclude=mulan-mc[01-02]

h2_add=30
h2_gxe=5
h2_nxe=15
sample_size=100000
num_envi_power=30
power_snp_h2_gxe=0.08


cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_fastGWA_feature/

k=0
for rep in {1..100}; do
	for envi_order in {1..40};do
		
		
		
		arr=()
		while read -r col1 col2; do
			arr+=("$col1")
		done < /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQplot_vf/variant_function.txt

		for val in "${arr[@]}"; do
			let k=${k}+1
			if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
				
				prefix=variant_function.$val.power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.h2_add_$h2_add.h2_gxe_$h2_gxe.h2_nxe_$h2_nxe.sample_$sample_size
				gcta --thread-num 5 --fastGWA-mlm \
									--bfile ../pheno_feature/$prefix.rep_$rep.power_snps \
									--pheno ../pheno_fastGWA_feature/$prefix.rep_$rep.txt \
									--grm-sparse /net/zootopia/disk1/chaon/data/UKB/GRM/gcta_imp_sp \
									--envir ../pheno_fastGWA_feature/E.$envi_order.txt \
									--geno 0.1 --maf 0.0001 \
									--qcovar ../pheno_fastGWA_feature/cov.txt \
									--out $prefix.$envi_order.$rep
				awk 'BEGIN{FS="[ \t]+";OFS="\t"} {print $1,$2,$3,$17}' $prefix.$envi_order.$rep.fastGWA \
								> $prefix.$envi_order.$rep.fastGWA.res
						rm $prefix.$envi_order.$rep.fastGWA
			
			fi
		done
		
		
	done
done

