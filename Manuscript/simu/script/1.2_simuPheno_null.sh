#!/bin/bash
#SBATCH --job-name=simu
#SBATCH --array=1-1200%50
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/simu_null.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/simu_null.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan

k=0

for rep in {1..100}; do
    let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
        python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/1.2_simuPheno_null.py \
        30 5 15 100000 $rep
    fi

    for h2_add in 5 50; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/1.2_simuPheno_null.py \
            $h2_add 5 15 100000 $rep
        fi
    done

    for h2_gxe in 1 10; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/1.2_simuPheno_null.py \
            30 $h2_gxe 15 100000 $rep
        fi
    done

    for h2_nxe in 0 5 25; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/1.2_simuPheno_null.py \
            30 5 $h2_nxe 100000 $rep
        fi
    done

    for size in 50000 200000 300000 400000; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/1.2_simuPheno_null.py \
            30 5 15 $size $rep
        fi
    done
done
