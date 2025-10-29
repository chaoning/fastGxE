#!/bin/bash
#SBATCH --job-name=IINT
#SBATCH --array=1-36%36
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Mental_health/Aout/I/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Mental_health/Aerr/I/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan

k=0

tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/description/trait_Description.txt | while read -r line; do
    trait=$(echo "$line" | awk '{print $1}')
    
    let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
        echo $trait
        python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Mental_health/3.IINT.py \
        /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/pheno/$trait.txt \
        /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/mental_health/INT/$trait.txt
    fi
    
done
