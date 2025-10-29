#!/bin/bash
#SBATCH --job-name=MoM
#SBATCH --array=1-64%30
#SBATCH --cpus-per-task=20
#SBATCH --time=7-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Aout/MoM.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/Aerr/MoM.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan



k=0

cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/mom/

tail -n +2 /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/trait32.txt | while read -r line;
do
    trait=$(echo "$line" | awk '{print $1}')

    let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
        time fastgxe --mom  \
                --bfile /net/zootopia/disk1/chaon/data/UKB/chip/ukb.keep.qc \
                --data /net/zootopia/disk1/chaon/WORK/GxE/pheno/Physical_measures/IINT/pheno.e42.$trait.txt --trait $trait \
                --env-int Age:Confide \
                --threads 20 \
                --out $trait --num-random-vector 30
    fi
	
	let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
        time fastgxe --mom  \
                --bfile /net/zootopia/disk1/chaon/data/UKB/chip/ukb.keep.qc \
                --data /net/zootopia/disk1/chaon/WORK/GxE/pheno/Physical_measures/IINT/pheno.e42.$trait.txt --trait $trait \
                --env-int Age:Confide \
                --threads 20 \
                --out $trait\_noNxE --no-noisebye --num-random-vector 30
    fi
    
done
