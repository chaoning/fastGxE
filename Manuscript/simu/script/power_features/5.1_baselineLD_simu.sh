#!/bin/bash
#SBATCH --job-name=simu
#SBATCH --array=1-100%100
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/Aout/baselineLD_null/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/Aerr/baselineLD_null/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


PYTHON="/net/zootopia/disk1/chaon/software/anaconda3/bin/python"

k=0

for rep in {1..100}; do
    let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
        "${PYTHON}" /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/5.1_baselineLD_simu.py \
        30 5 15 100000 $rep
    fi
done
