
cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/FUMA/

python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/FUMA/2.bed_file.py \
GenomicRiskLoci.csv \
GenomicRiskLoci.bed


liftOver /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/FUMA/GenomicRiskLoci.bed \
    /net/zootopia/disk1/chaon/data/human_ref/hg19ToHg38.over.chain.gz \
    GenomicRiskLoci.hg38.bed \
    GenomicRiskLoci.hg38.unmapped.bed


cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/

python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/FUMA/2.bed_file.py \
GenomicRiskLoci.csv \
GenomicRiskLoci.bed


liftOver GenomicRiskLoci.bed \
    /net/zootopia/disk1/chaon/data/human_ref/hg19ToHg38.over.chain.gz \
    GenomicRiskLoci.hg38.bed \
    GenomicRiskLoci.hg38.unmapped.bed


cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/FUMA_log/

python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/FUMA/2.bed_file.py \
GenomicRiskLoci.csv \
GenomicRiskLoci.bed


liftOver GenomicRiskLoci.bed \
    /net/zootopia/disk1/chaon/data/human_ref/hg19ToHg38.over.chain.gz \
    GenomicRiskLoci.hg38.bed \
    GenomicRiskLoci.hg38.unmapped.bed

