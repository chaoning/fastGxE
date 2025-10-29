
cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/FUMA/

for maxDist in 0 300000;
do
    python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/FUMA/3.addGenes.hg38.py \
    $maxDist 1 \
    GenomicRiskLoci.hg38.csv \
    GenomicRiskLoci.hg38.addGene.$maxDist.protein1.csv &
done

for maxDist in 0 300000;
do
    python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/FUMA/3.addGenes.hg38.py \
    $maxDist 0 \
    GenomicRiskLoci.hg38.csv \
    GenomicRiskLoci.hg38.addGene.$maxDist.protein0.csv &
done


cd /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/biomarker/results/FUMA/

for maxDist in 0 300000;
do
    python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/FUMA/3.addGenes.hg38.py \
    $maxDist 1 \
    GenomicRiskLoci.hg38.csv \
    GenomicRiskLoci.hg38.addGene.$maxDist.protein1.csv &
done

for maxDist in 0 300000;
do
    python /net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/FUMA/3.addGenes.hg38.py \
    $maxDist 0 \
    GenomicRiskLoci.hg38.csv \
    GenomicRiskLoci.hg38.addGene.$maxDist.protein0.csv &
done
