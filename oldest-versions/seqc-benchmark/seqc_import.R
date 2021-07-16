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



DqCT<-function(data, housekeeper){
  hk_row<-grep(housekeeper,data$Symbol) #Retrieve the row(s) for the housekeeper
  norm_data<-data.frame(Symbol=data[,"Symbol"])
  for(column in (3:ncol(data))){
    hk_average<-mean(data[hk_row,column])
    new_column<-c()
    for(gene in (1:nrow(data))){
      new_column<-c(new_column,(data[gene,column]-hk_average))
    }
    norm_data[,colnames(data)[column]]<-new_column
}
  return(norm_data)
}


norm_data<-DqCT(qPCR_exp, "GAPDHS")

DE_analysis<-function(data, groups){
  
  #Create a df with group name | Sample idxs
  groups_mapping<-data.frame(group=character(),idx=numeric())
  for(gr in groups){
    group_samples<-grep(gr, colnames(data)) #Retrieve the column for the samples
    group_samples = group_samples[!(group_samples %in% c(1:2))] #Remove 2 first cols Entrez/Symbol
    groups_mapping<-rbind(groups_mapping,data.frame(gr,group_samples))  
  }
  
  #Perform the ANOVA analysis for each gene
  n_col<-(factorial(length(groups))/(factorial(length(groups)-2)*2))+3
  comparisonsDF<-data.frame(matrix(NA,  ncol = n_col))
  comparisonsDF<-comparisonsDF[-1,]
  
  for(gene in c(1:nrow(data))){
    temp_df<-data.frame(group=character(),idx=numeric(),expr=numeric())
    temp_vec<-data[gene,groups_mapping$group_samples]
    temp_df<-cbind(groups_mapping,t(temp_vec))
    
    aov_results<-aov(temp_df[,3] ~ temp_df$gr)
    aov_p<-summary(aov_results)[[1]][["Pr(>F)"]][1]
    
    tukres<-TukeyHSD(aov_results)
    tukey_fdr<-p.adjust(tukres$`temp_df$gr`[,4],method='fdr')
    newrow<-data[gene,c(1,2)]
    newrow<-c(newrow,aov_p)
    newrow<-c(newrow,tukey_fdr)
    comparisonsDF<-rbind(setNames(comparisonsDF,names(newrow)),newrow)
  }
  comparisonsDF$p_adj<-p.adjust(comparisonsDF[,3],method='fdr')
  return(comparisonsDF)
}

DE<-DE_analysis(norm_data, groups = sample_lets)

