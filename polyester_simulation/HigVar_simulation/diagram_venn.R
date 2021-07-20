setwd("~/Documents/16S_rrna/HigVar_simulation")
edgeR<- read.csv("edgeR_high_var.csv")
Deseq2<- read.csv("Deseq2_high_vars.csv")
stats_16 <- read.delim("comparisons_compound_or_gene.txt")
aldex2<-read.csv("Aldo_result_high_var.csv")


degs16S<-stats_16[stats_16$pvalue<0.05,c(1,4)]
degs_deseq2<-Deseq2[Deseq2$padj<0.05,c(1,7)]
degs_edgeR<-edgeR[edgeR$PValue<0.05,c(1,4)]
degs_aldex2<-aldex2[aldex2$wi.eBH<0.05,c(1,13)]

degsList<-list('16Srdna_staats'=degs16S$Data,
               'Deseq2'=degs_deseq2$X,
               'edgeR'=degs_edgeR$X,
               'Aldex2'=degs_aldex2$X)

library(ggvenn)

png("venn_diagram_high_var.png", units = "cm",width=25, height = 25, res=300)
ggvenn(degsList,c('16Srdna_staats','Deseq2','edgeR','Aldex2'))
dev.off()
