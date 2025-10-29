cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/FUMA/
sbatch -p mulan -c 40 --wrap="python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/FUMA/4.trait_specific.py \
GenomicRiskLoci.hg38.keygenes.csv \
p_gxe \
GenomicRiskLoci.hg38.keygenes.TraitSpecific.csv \
/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/signal/signal.txt 32 "


cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/
sbatch -p mulan -c 40 -w mulan-mc12 --wrap="python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/FUMA/4.trait_specific.py \
GenomicRiskLoci.hg38.keygenes.csv \
p_gxe \
GenomicRiskLoci.hg38.keygenes.TraitSpecific.csv \
/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/signal/signal.txt 67"

cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/FUMA_log/
sbatch -p mulan -c 40 -w mulan-mc12 --wrap="python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/FUMA/4.trait_specific.py \
GenomicRiskLoci.hg38.csv \
p_gxe_log \
GenomicRiskLoci.hg38.TraitSpecific.csv \
/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/testGxE_log/5_sum_signal.txt 32"
