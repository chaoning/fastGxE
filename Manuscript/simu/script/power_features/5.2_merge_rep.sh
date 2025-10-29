#!/bin/bash
#SBATCH --job-name=simu
#SBATCH --array=1-1%1
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/Aout/merge/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/Aerr/merge/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan

k=0


    let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
        python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/5.2_merge_rep.py \
        h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
    fi