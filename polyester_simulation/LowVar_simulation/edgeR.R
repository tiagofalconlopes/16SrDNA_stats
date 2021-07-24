library(edgeR)

start_time <- Sys.time()

load("~/Documents/16S_rrna/LowVar_simulation/sim_counts_matrix.rda")
sim.low.var<-as.data.frame(counts_matrix)
mapping_file <- read.delim("~/Documents/16S_rrna/LowVar_simulation/mapping_file.txt")
mapping_file$description<-as.factor(mapping_file$description)

geneList<-DGEList(counts = sim.low.var, group = mapping_file$description)
normFactors <- calcNormFactors(geneList)
design <- model.matrix(sample~description, data = mapping_file)

geneList_norm <- estimateDisp(normFactors,design)
et <- exactTest(geneList_norm)
end_time <- Sys.time()
end_time - start_time

#write.csv(et$table,"edgeR_low_var.csv")
