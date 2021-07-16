library(DESeq2)

start_time <- Sys.time()
#load("~/Documents/16S_rrna/LowVar_simulation/sim_counts_matrix.rda")
#sim.low.var<-as.data.frame(counts_matrix)
sim.alta.var <- read.csv("~/Documents/16S_rrna/HigVar_simulation/sim-alta-var.csv", row.names = 1)
mapping_file <- read.delim("~/Documents/16S_rrna/HigVar_simulation/mapping_file.txt")
mapping_file$description<-as.factor(mapping_file$description)

deseq_data<-DESeqDataSetFromMatrix(sim.alta.var, colData = mapping_file, design = ~description)
deseq_analysis <- DESeq(deseq_data)
#write.csv(data.frame(results(deseq_analysis)), file="DEGS_low_vars.csv")
end_time <- Sys.time()
end_time - start_time
