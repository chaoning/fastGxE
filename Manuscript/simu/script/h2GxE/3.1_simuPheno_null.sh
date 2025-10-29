#!/bin/bash
#SBATCH --job-name=simu
#SBATCH --array=1-400%400
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/simu_null_indE.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/simu_null_indE.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan

k=0

for rep in {1..100}; do
    for h2_nxe in 0 5 15 25; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/h2GxE/3.1_simuPheno_null.py \
            30 5 $h2_nxe 100000 $rep
        fi
    done
done
