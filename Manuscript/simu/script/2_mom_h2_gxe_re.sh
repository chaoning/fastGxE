#!/bin/bash
#SBATCH --job-name=mom
#SBATCH --array=[440,491]
#SBATCH --cpus-per-task=10
#SBATCH --time=2-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/mom_h2_gxe.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/mom_h2_gxe.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan

k=0
cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_fastgxe/

for rep in {1..100}; do
    
    for h2_gxe in 1 5 10; do
	
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			
            prefix=h2_add_30.h2_gxe_$h2_gxe.h2_nxe_15.sample_100000
			time fastgxe --mom  \
                --bfile /net/zootopia/disk1/chaon/data/UKB/chip/ukb.keep.qc \
                --data $prefix.txt --trait trait$rep \
                --env-int age:alcohol_frequency \
                --threads 10 \
                --out ../res_mom/$prefix.rep_$rep
				
        fi
		
		let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            prefix=h2_add_30.h2_gxe_$h2_gxe.h2_nxe_15.sample_100000
			
			time fastgxe --mom  \
                --bfile /net/zootopia/disk1/chaon/data/UKB/chip/ukb.keep.qc \
                --data $prefix.txt --trait trait$rep \
                --env-int age:alcohol_frequency \
                --threads 10 \
                --out ../res_mom/$prefix.rep_$rep.noNxE --no-noisebye
        fi
		
    done
	
done
