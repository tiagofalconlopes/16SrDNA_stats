library(seqc)

#Gather the dataset for gene expression, using refseq annotation, sequenced on Illumina platform
#on BGI center
rna_seq_exp_bank=ILM_refseq_gene_BGI
qPCR_exp=taqman

#Pick two samples of each replicate from the rna-seq analysis, random lane and flowcell
#colname= (sample letter A-D)_(replicate number 1-5)_L0(lane number 1-8)_FlowCell(fc letter A-B)
sample_lets<-c('A','B','C','D')
set.seed(10)
rna_col_names<-c('EntrezID','Symbol','IsERCC')
for(s in sample_lets){
  for(rep in c(1:4)){
    lane=sample(1:8,2)
    for(l in lane){
      fc=sample(c('A','B'),1)
      s_name=paste(s,'_',rep,'_','L0',l,'_FlowCell',fc, sep='')
      rna_col_names<-c(rna_col_names, s_name)    
    }
  }
}

expr_matrix<-rna_seq_exp_bank[,rna_col_names]

####TAQMAN analysis####

#Remove detection flags
qPCR_exp<-qPCR_exp[,grep("detection",colnames(qPCR_exp), invert = TRUE)]
row.names(qPCR_exp)<-paste(row.names(qPCR_exp),qPCR_exp$Symbol,sep = '_')


# housekeepers_symbols<-c("UBC","ACTB", "ACTN2","TOP1","CYC1","GAPDH","18S")
# #Add um gene aleatorio pra criar parametro de variavao
# housekeepers_symbols<-c(housekeepers_symbols,"MDM4")
# housekeepers_expr<-na.omit(qPCR_exp[match(housekeepers_symbols, qPCR_exp$Symbol, nomatch = NULL),])
# heatmap(as.matrix(housekeepers_expr[,3:18]),Colv = NA)

#https://bitesizebio.com/24894/4-easy-steps-to-analyze-your-qpcr-data-using-double-delta-ct-analysis/
DDqCT<-function(data, housekeeper, groups){
  hk_row<-grep(housekeeper,data$Symbol) #Retrieve the row(s) for the housekeeper
  for(gr in groups){
    group_samples<-grep(gr, colnames(data)) #Retrieve the column for the samples
    group_samples = group_samples[!(group_samples %in% c(1:2))] #Remove 2 first cols Entrez/Symbol
    hk_average<-colMeans(data[hk_row,group_samples])
    hk_average<-mean(hk_average)
    #Calculate mean of each gene - mean of the HK
    data[,paste("deltaCt",gr)]<-apply(data[,group_samples], MARGIN=1,hk_average=hk_average,FUN = function(x, hk_average){mean(x)-hk_average})
  }
  return(data)
}

norm_data<-DDqCT(qPCR_exp, "UBC", sample_lets)

