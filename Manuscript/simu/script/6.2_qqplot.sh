#!/bin/bash
#SBATCH --job-name=QQ
#SBATCH --array=1-20%20
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/QQ_/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/QQ_/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_fastgxe/
k=0


for num in {1..6};do

	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		
		prefix=h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
		python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/6.2_qqplot.py $prefix \
		/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQplot_maf/mafbin_$num.frq

	fi
done



arr=()
while read -r col1 col2; do
    arr+=("$col1")
done < /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQplot_vf/variant_function.txt



for val in "${arr[@]}"; do
    echo "$val"
	
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		
		prefix=h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
		python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/6.2_qqplot.py $prefix \
		/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQplot_vf/variant_function.$val.txt

	fi
	
done

for num in {1..5};do

	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		
		prefix=h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
		python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/6.2_qqplot.py $prefix \
		/net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/QQplot_ldscore/ldscore_bin_$num.txt

	fi
done

