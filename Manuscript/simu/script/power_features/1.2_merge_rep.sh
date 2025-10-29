#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --array=1-20%20
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/Aout/simu_power/merge.%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/Aerr/simu_power/merge.%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


h2_add=30
h2_gxe=5
h2_nxe=15
sample_size=100000


k=0

num_envi_power=30

power_snp_h2_gxe=0.08

for nBin in {1..6}; do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		prefix=maf_$nBin.power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.h2_add_$h2_add.h2_gxe_$h2_gxe.h2_nxe_$h2_nxe.sample_$sample_size
		python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/1.2_merge_rep.py $prefix
	
	fi
done


arr=()
while read -r col1 col2; do
    arr+=("$col1")
done < /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQplot_vf/variant_function.txt

for val in "${arr[@]}"; do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		echo "$val"
		prefix=variant_function.$val.power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.h2_add_$h2_add.h2_gxe_$h2_gxe.h2_nxe_$h2_nxe.sample_$sample_size
		python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/1.2_merge_rep.py $prefix
	
	fi
done


for nBin in {1..5}; do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		prefix=ld_$nBin.power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.h2_add_$h2_add.h2_gxe_$h2_gxe.h2_nxe_$h2_nxe.sample_$sample_size
		python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_features/1.2_merge_rep.py $prefix

	fi
done
