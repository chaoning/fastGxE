#!/bin/bash
#SBATCH --job-name=simu
#SBATCH --array=1-12%12
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/simu_null_chip.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/simu_null_chip.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan

k=0


    let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
        python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/1.3_merge_rep_chip.py \
        h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
    fi

    for h2_add in 5 50; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/1.3_merge_rep_chip.py \
            h2_add_$h2_add.h2_gxe_5.h2_nxe_15.sample_100000
        fi
    done

    for h2_gxe in 1 10; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/1.3_merge_rep_chip.py \
            h2_add_30.h2_gxe_$h2_gxe.h2_nxe_15.sample_100000
        fi
    done

    for h2_nxe in 0 5 25; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/1.3_merge_rep_chip.py \
            h2_add_30.h2_gxe_5.h2_nxe_$h2_nxe.sample_100000
        fi
    done

    for size in 50000 200000 300000 400000; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/1.3_merge_rep_chip.py \
            h2_add_30.h2_gxe_5.h2_nxe_15.sample_$size
        fi
    done
