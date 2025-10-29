#!/bin/bash
#SBATCH --job-name=mmSuSiE
#SBATCH --array=1-1300%1300
#SBATCH --cpus-per-task=5
#SBATCH --time=1-00:00:00
#SBATCH --output=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_mmSuSiE/Aout/mmSuSiE/%a.txt
#SBATCH --error=/net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_mmSuSiE/Aerr/mmSuSiE/%a.txt
#SBATCH --mem-per-cpu=5000MB
#SBATCH --partition=mulan



k=0
for rep in {1..100}; do
	
	num_envi_power=30
	
	for power_snp_h2_gxe in 0.02 0.04 0.06 0.08 0.1 0.12 0.14 ; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_mmSuSiE/mmSuSiE.py $power_snp_h2_gxe $num_envi_power $rep
		
		fi
	done
	

	power_snp_h2_gxe=0.08
	
	for num_envi_power in 1 2 5 10 20 40; do
		let k=${k}+1
		if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		
			python /net/zootopia/disk1/chaon/WORK/GxE/simuR1/script/power_mmSuSiE/mmSuSiE.py $power_snp_h2_gxe $num_envi_power $rep
			
		fi
	done

done
