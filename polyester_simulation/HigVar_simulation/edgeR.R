library(edgeR)

sim.alta.var <- read.csv("~/Documents/16S_rrna/HigVar_simulation/sim-alta-var.csv", row.names = 1)
mapping_file <- read.delim("~/Documents/16S_rrna/HigVar_simulation/mapping_file.txt")
mapping_file$description<-as.factor(mapping_file$description)

geneList<-DGEList(counts = sim.alta.var, group = mapping_file$description)
normFactors <- calcNormFactors(geneList)
design <- model.matrix(sample~description, data = mapping_file)

geneList_norm <- estimateDisp(normFactors,design)
et <- exactTest(geneList_norm)
write.csv(et$table,"edgeR_high_var.csv")
