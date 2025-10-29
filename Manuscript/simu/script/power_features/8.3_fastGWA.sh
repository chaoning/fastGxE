#!/bin/bash
#SBATCH --job-name=fastGWA
#SBATCH --array=16635-28000%1000
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/Aout/power/fastGWA_baselineLDpower/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/Aerr/power/fastGWA_baselineLDpower/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan,main
#SBATCH --exclude=mulan-mc[01-03]

h2_add=30
h2_gxe=5
h2_nxe=15
sample_size=100000

cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/res_baselineLD_fastGWA/

k=0
for rep in {1..100}; do
	for envi_order in {1..40};do
	
		num_envi_power=30
		
		for power_snp_h2_gxe in 0.02 0.04 0.06 0.08 0.1 0.12 0.14; do
			let k=${k}+1
			if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
				prefix=power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.h2_add_$h2_add.h2_gxe_$h2_gxe.h2_nxe_$h2_nxe.sample_$sample_size
				
				gcta --thread-num 5 --fastGWA-mlm \
									--bfile ../pheno_baselineLD/$prefix.rep_$rep.power_snps \
									--pheno ../pheno_baselineLD_fastGWA/$prefix.rep_$rep.txt \
									--grm-sparse /net/zootopia/disk1/chaon/data/UKB/GRM/gcta_imp_sp \
									--envir ../pheno_baselineLD_fastGWA/E.$envi_order.txt \
									--geno 0.1 --maf 0.0001 \
									--qcovar ../pheno_baselineLD_fastGWA/cov.txt \
									--out $prefix.$envi_order.$rep
				awk 'BEGIN{FS="[ \t]+";OFS="\t"} {print $1,$2,$3,$17}' $prefix.$envi_order.$rep.fastGWA \
								> $prefix.$envi_order.$rep.fastGWA.res
						rm $prefix.$envi_order.$rep.fastGWA
			
			fi
		done
		
	done
done
