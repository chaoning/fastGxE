library(locuszoomr)
library(EnsDb.Hsapiens.v75)
library(dplyr)
library(ensembldb)
library(AnnotationForge)
library(AnnotationHub)
ah <- AnnotationHub()




setwd("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/results/locuszoom/data/")

df <- read.csv("/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/physical/results/FUMA/GenomicRiskLoci.hg38.keygenes.TraitSpecific.addP.csv", head=T)
ensDb_v111 <- ah[["AH116291"]]
print(ensembldb::listGenebiotypes(ensDb_v111))


for (i in 1:dim(df)[1]) {
    loci <- df$GenomicLocus_hg38[i]
    start1 <- df$start_hg38[i]
    end1 <- df$end_hg38[i]
    trait <- df$trait[i]
    snp <- df$trait_leading_snp[i]
    gene <- df$Genes[i]
    TraitName <- df$TraitName[i]
    up <- df$up[i]
    down <- df$down[i]
    
    in_file = paste0("PT", ".", loci, ".", trait, ".", snp)
    print(in_file)
    file = paste0(in_file, ".GWAS.txt")
    data <- read.table(file, head = TRUE)
    file <- paste0(in_file, ".r2.txt")
    data2 <- read.table(file, head = TRUE)
    data$r2 <- data2$r2

    data <- data %>%
          mutate(p_gxe = ifelse(p_gxe < 1.0e-30, 1.0e-30, p_gxe))

    
    loc <- locus(data = data, chrom="chrom", pos="base", p="p_gxe", gene=gene,
             ens_db = ensDb_v111, flank = c(up, down), LD = "r2")
    
    file = paste("../fig/", in_file, ".protein_coding.pdf", sep="")
    log_min_value <- -log10(min(data$p_gxe, na.rm = TRUE))
    
    upper_limit <- min(log_min_value, 50)
    
    pdf(file = file, width = 7/2.54, height = 6.5/2.54)
    pf <- quote({
      abline(v = start1, col = "darkgrey", lwd = 1.5, lty = 2)
      abline(v = end1, col = "darkgrey", lwd = 1.5, lty = 2)
    })
    locus_plot(loc, pcutoff=5.0e-8/32, legend_pos = 'topright',
    filter_gene_biotype = c("protein_coding"), 
    labels = c(snp),
          gene_col = '#376795', 
          exon_col = '#ef8a47', exon_border = '#ef8a47',
          cex = 0.6,
       cex.axis = 0.7,
       cex.lab = 0.8,
       cex.text = 0.6, ylim = c(0, upper_limit), panel.last = pf)
    dev.off()
}
