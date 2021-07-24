library(ALDEx2)

sim.alta.var <- read.csv("~/Documents/16S_rrna/HigVar_simulation/sim-alta-var.csv", row.names = 1)
mapping_file <- read.delim("~/Documents/16S_rrna/HigVar_simulation/mapping_file.txt")
mapping_file$description<-as.factor(mapping_file$description)

mtx_CX <- as.matrix(sim.alta.var[,c(1:ncol(sim.alta.var))])
Xald<-aldex(mtx_CX,mapping_file$description)
df <- data.frame(row.names(Xald),Xald)
write.csv(df,"Aldo_result_high_var.csv")
