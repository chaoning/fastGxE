#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --array=1-99%99
#SBATCH --cpus-per-task=5
#SBATCH --time=7-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/signal_vs/Aout/merge/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/signal_vs/Aerr/merge/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan



k=0

for index in {0..31};
do
    let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/signal_vs/1.1_PT.py $index
    fi
    
done

for index in {0..66};
do
    let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/signal_vs/1.2_BB.py $index
    fi
    
done
