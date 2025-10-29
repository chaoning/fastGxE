#!/bin/bash
#SBATCH --job-name=mom
#SBATCH --array=[22,23,24,25,206,208,288,298,303,309]
#SBATCH --cpus-per-task=10
#SBATCH --time=2-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/mom_h2_nxe.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/mom_h2_nxe.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan

k=0
cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_fastgxe/

for rep in {1..100}; do
    
    for h2_nxe in 0 5 25; do
	
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			
            prefix=h2_add_30.h2_gxe_5.h2_nxe_$h2_nxe.sample_100000
			time fastgxe --mom  \
                --bfile /net/zootopia/disk1/chaon/data/UKB/chip/ukb.keep.qc \
                --data $prefix.txt --trait trait$rep \
                --env-int age:alcohol_frequency \
                --threads 10 \
                --out ../res_mom/$prefix.rep_$rep
				
        fi
		
		let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            prefix=h2_add_30.h2_gxe_5.h2_nxe_$h2_nxe.sample_100000
			
			time fastgxe --mom  \
                --bfile /net/zootopia/disk1/chaon/data/UKB/chip/ukb.keep.qc \
                --data $prefix.txt --trait trait$rep \
                --env-int age:alcohol_frequency \
                --threads 10 \
                --out ../res_mom/$prefix.rep_$rep.noNxE --no-noisebye
        fi
		
    done
	
done
