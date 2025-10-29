#!/bin/bash
#SBATCH --job-name=GENIE
#SBATCH --array=1-8%8
#SBATCH --cpus-per-task=10
#SBATCH --time=2-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/GENIE_data/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/GENIE_data/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan

bash

let k=0


for add in 5 30 50; do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		prefix=h2_add_$add.h2_gxe_5.h2_nxe_15.sample_100000
		python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/2.1_data.py $prefix
	fi
done

for gxe in 1 10; do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		prefix=h2_add_30.h2_gxe_$gxe.h2_nxe_15.sample_100000
		python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/2.1_data.py $prefix
	fi
done

for nxe in 0 5 25; do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		prefix=h2_add_30.h2_gxe_5.h2_nxe_$nxe.sample_100000
		python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/2.1_data.py $prefix
	fi
done

