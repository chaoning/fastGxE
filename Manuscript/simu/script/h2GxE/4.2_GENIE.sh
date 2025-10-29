#!/bin/bash
#SBATCH --job-name=GENIE
#SBATCH --array=1-400%200
#SBATCH --cpus-per-task=10
#SBATCH --time=2-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aout/GENIE/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/Aerr/GENIE/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan,main
#SBATCH --exclude=mulan-mc[01-02]

bash

let k=0

gen=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink.100K.chip
annot=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink.100K.chip.anno

cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/null/pheno_GENIE/
cov=0.cov.txt

for NxE in 0 5 15 25; do
	for rep in {1..100};do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			prefix=h2_add_30.h2_gxe_5.h2_nxe_$NxE.sample_100000
			phen=indE.$prefix.$rep.txt
			env=indE.$prefix.E.$rep.txt
			out=../res_GENIE/indE.$prefix.$rep.G_GxE_NxE.out2
			GENIE -g $gen -p $phen -c $cov -e $env -m G+GxE+NxE  -k 10 -jn 2  -o $out --annot $annot -t 10 -s 1
		fi
	done
done
