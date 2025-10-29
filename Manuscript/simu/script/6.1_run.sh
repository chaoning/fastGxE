#!/bin/bash
#SBATCH --job-name=QQ
#SBATCH --array=1-8%8
#SBATCH --cpus-per-task=20
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/QQ/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/QQ/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan


cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_fastgxe/
k=0

    let k=${k}+1
    if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		
		prefix=h2_add_30.h2_gxe_5.h2_nxe_15.sample_100000
		python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/6.1_qqplot.py $prefix
	
    fi

    for h2_add in 5 50; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            prefix=h2_add_$h2_add.h2_gxe_5.h2_nxe_15.sample_100000
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/6.1_qqplot.py $prefix
        fi
    done

    for h2_gxe in 1 10; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            prefix=h2_add_30.h2_gxe_$h2_gxe.h2_nxe_15.sample_100000
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/6.1_qqplot.py $prefix
        fi
    done

    for h2_nxe in 0 5 25; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
            prefix=h2_add_30.h2_gxe_5.h2_nxe_$h2_nxe.sample_100000
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/6.1_qqplot.py $prefix
        fi
    done
