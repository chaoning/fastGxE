
sbatch -p mulan -c 20 --wrap="python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/FUMA/7.1_add_logP.py \
/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP \
/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/testGxE/"

sbatch -p mulan -c 80 -w mulan-mc12 --wrap="python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/FUMA/7.1_add_logP.py \
/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP \
/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/log_trans/testGxE/"
