#!/bin/bash
#SBATCH --job-name=fastGWA
#SBATCH --array=1-12000%500
#SBATCH --cpus-per-task=5
#SBATCH --time=2-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/fastGWA/h2_gxe.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/fastGWA/h2_gxe.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan,main
#SBATCH --exclude=mulan-mc[01-03]

bash
cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/res_fastGWA/
k=0
for rep in {1..100}; do
	for envi_order in {1..40};do
	
		for h2_gxe in 1 5 10; do
			let k=${k}+1
			if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
				prefix=h2_add_30.h2_gxe_$h2_gxe.h2_nxe_15.sample_100000
				gcta --thread-num 5 --fastGWA-mlm \
								--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
								--pheno ../pheno_fastGWA/$prefix.rep_$rep.txt \
								--grm-sparse /net/zootopia/disk1/chaon/data/UKB/GRM/gcta_imp_sp \
								--envir ../pheno_fastGWA/E.$envi_order.txt \
								--geno 0.1 --maf 0.0001 \
								--qcovar ../pheno_fastGWA/cov.txt \
								--out $prefix.$envi_order.$rep
			awk 'BEGIN{FS="[ \t]+";OFS="\t"} {print $1,$2,$3,$17}' $prefix.$envi_order.$rep.fastGWA \
                            > $prefix.$envi_order.$rep.fastGWA.res
                    rm $prefix.$envi_order.$rep.fastGWA
                    gzip -f $prefix.$envi_order.$rep.fastGWA.res
			fi
		done
		
	done
	
done
