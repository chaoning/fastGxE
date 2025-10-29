#!/bin/bash
#SBATCH --job-name=structlmm
#SBATCH --array=1-400%200
#SBATCH --cpus-per-task=3
#SBATCH --time=2-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/structlmm_size/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/structlmm_size/%a.txt
#SBATCH --mem-per-cpu=8000MB
#SBATCH --partition=mulan
#SBATCH --exclude=mulan-mc[01-03]


cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_fastgxe/
k=0
for rep in {1..100}; do
	
    for size in 50000 200000 300000 400000; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            prefix=h2_add_30.h2_gxe_5.h2_nxe_15.sample_$size
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/5.2_structlmm.py $prefix $size $rep
        fi
    done
	
done
