#!/bin/bash
#SBATCH --job-name=zoom
#SBATCH --array=1-300%300
#SBATCH --cpus-per-task=5
#SBATCH --time=7-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/locuszoom/Aout/BB.zoom.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/locuszoom/Aerr/BB.zoom.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan

let k=0

for val in {1..300};
do
    
    let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
        python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/locuszoom/2.locuszoom_data_BB.py $val
    fi
    
done
