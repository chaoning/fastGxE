#!/bin/bash
#SBATCH --job-name=fastgxe
#SBATCH --array=2128-24500%300
#SBATCH --cpus-per-task=3
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/take_home/Aout/fastgxe_testgxe/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/take_home/Aerr/fastgxe_testgxe/%a.txt
#SBATCH --mem-per-cpu=6000MB
#SBATCH --partition=mulan
#SBATCH --exclude=mulan-mc[01-03]

cd /net/zootopia/disk1/chaon/WORK/GxE/simuR1/power/pheno_take_home_fastgxe/
h2_add=30
h2_gxe=5
h2_nxe=15

k=0
for rep in {1..100}; do
	for num_envi_power in 1 2 5 10 20 30 40; do
		for power_snp_h2_gxe in 0.02 0.04 0.06 0.08 0.1 0.12 0.14; do
			for sample_size in 50000 100000 200000 300000 400000; do
				let k=${k}+1
				if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
					prefix=power_snp_h2_gxe_$power_snp_h2_gxe.num_envi_$num_envi_power.h2_add_$h2_add.h2_gxe_$h2_gxe.h2_nxe_$h2_nxe.sample_$sample_size
					bim_file=../pheno_take_home/$prefix.rep_$rep.power_snps.bim
					snp_start=$(awk 'NR==1 {print $2}' $bim_file)
					snp_end=$(awk 'END {print $2}' $bim_file)
					time fastgxe --test-gxe  \
						--bfile /net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/plink22 \
						--grm /net/zootopia/disk1/chaon/data/UKB/GRM/ukb_imp \
						--data $prefix.txt --trait trait$rep \
						--env-int age:alcohol_frequency \
						--threads 10 --num-random-snp 1000 --maxiter 100 \
						--out ../res_fastgxe_testgxe_take_home/$prefix.rep_$rep \
						--snp-range $snp_start $snp_end
				
				fi
			done
		done
	done
done
