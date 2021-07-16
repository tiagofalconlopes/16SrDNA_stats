##################################################################################################################################
#############################################16S rRNA profilling downstream analysis##############################################
##################################################################################################################################
#
start_time <- Sys.time()
#
# This pipeline assumes that your data has two or more groups of samples!!!
#
# Adapted for one dependent variable column (compound/gene list of variables)
#
#
#
# Necessary packages
library("psych") # For the geometric mean calculation
library ("pvclust") # For the hierarchical clustering
library("gplots") # For the heat maps
library("PMCMRplus") # For non parametric pairwised tests
library("ggbiplot") # For PCA plot
library("vegan") # For PERMANOVA analysis and plot, and for diversity indexes
# This is all you need to change:
work_directory <- "/Users/tiagofalcon/Desktop/socorro" # Path to your work directory. This is 
# the directory where you may find the two necessary files to begin the analysis and where you gonna save the outputs.
file_name <- "GSE95224_counts.spf" # Change the file name. It is a tab delimited OTU table.
map_file_data <- "sim_rep_info.txt" # Mapping file with groups description names.
#
##Import the data from .rdata
start_time <- Sys.time()
load("~/Documents/16S_rrna/LowVar_simulation/sim_counts_matrix.rda")
input<-as.data.frame(counts_matrix)
input$Level<-row.names(input)
input<-input[,c(25,c(1:24))]


#
#
#create a transformed_df backup to modify the colnames by color names to use as side colors ...
#
#
#### Aitchison's log ratio transformation...
#
# Go to the work directory.
#setwd(work_directory)
# Read the input file... In this case we used the OTU table with raw 
# reads count per sample from a 16S profilling analysis. The standard 
# file is the 16S BMP pipeline output OTU table. Could be either the TXT
# or SPF file.
#input <- read.delim(file_name)
#
map_groups <- read.delim(map_file_data)
map_groups$lib_sizes<-NULL
colnames(map_groups)<-c("Sample.ID","Description")
map_groups$Description<-as.factor(map_groups$Description)
#
#
# Now let's evaluate the different compounds. For this purpose, we gonna 
# exclude unclassified taxa. Aitchison's transformation should be done again.
# To change to others taxonomic leves, modify the "input$level_N" column and
# the name of the taxon level.
#
# Replace all zeros by a really small value (0.0000001). This step is 
# necessary once Aitchinson's log ratio transformation do not accept zeros.
# This value is kind of arbritary. As there are many cells with "zeros"
# and they are not fixed to specific compositions (variables), we opted
# to use a very low value in place of the zero.
input[input == 0] <- 0.0000001
# Retrieve only the values
numbers_to_norm <- input[,!lapply(input, class) == "character"]
# Transform each column to proportion.
prop_otus <- t(t(numbers_to_norm)/rowSums(t(numbers_to_norm)))
# Check if each column's sum is 1.
colSums(prop_otus)
# Convert prop_otus to data frame format.
df_prop_otus <- as.data.frame(prop_otus)
# Now, divide each cell in a column by the geometric mean of this 
# column and get the LN of this division...
# ... and get the interest values. The log transformation is ln.
gambiarra_methods <- function(x){
  for (column in 1:ncol(df_prop_otus)) {
    df_prop_otus[,column] <- log(df_prop_otus[,column]/geometric.mean(df_prop_otus[,column]))
  }
  transformed_data <- df_prop_otus
  return(transformed_data)
}
# Send the transformed data frame to a variable.
transformed_df <- gambiarra_methods(df_prop_otus)
# Put the taxonomic levels back ...
end_file <- cbind (input[,lapply(input, class) == "character"], transformed_df)#
# Make sure that values in data frame are numeric
for (i in 1:ncol(transposed_taxa)){
  transposed_taxa[,i] <- as.numeric(as.character(transposed_taxa[,i]))
}

transposed_taxa <- t(end_file)
colnames(transposed_taxa) <- transposed_taxa[1,]
transposed_taxa <- transposed_taxa[-1,]
transposed_taxa <- as.data.frame(transposed_taxa)
# Make sure that values in data frame are numeric
 for (i in 1:ncol(transposed_taxa)){
   transposed_taxa[,i] <- as.numeric(as.character(transposed_taxa[,i]))
 }
# Change samples per groups names
merged_df <- merge(map_groups, transposed_taxa, by.x='Sample.ID', by.y= 'row.names')
transposed_taxa <- merged_df[,c(2:ncol(merged_df))]
# Perform normality test and the adequate  comparisons
data_info <- c()
data_t <- c()
data_fdr <- c()
data_pvalue <-c()
data_info_wil <- c()
data_wil_w <- c()
data_wil_fdr <- c()
data_wil_pvalue<-c()


for (i in 2:ncol(transposed_taxa)){
  z<-colnames(transposed_taxa[i])
  x <- shapiro.test(transposed_taxa[,i])
  if (length(levels(map_groups$Description)) <= 2){
    warning("You don't have more than two groups")} else if (x$p.value >= 0.05 && length(levels(transposed_taxa$Description)) >= 3){
      tuk_res <- TukeyHSD(aov(transposed_taxa[,i] ~ transposed_taxa$Description))
      tuk_res$`transposed_taxa$Description`[,4] <- p.adjust(tuk_res$`transposed_taxa$Description`[,4],method = 'fdr')
      write.table(z, file = "comparisons_compound_or_gene_tukey.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
      write.table(tuk_res$`transposed_taxa$Description`, file = "comparisons_compound_or_gene_tukey.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    } else if (x$p.value < 0.05 && length(levels(map_groups$Description)) >= 3){
      krusk_res <- kwAllPairsDunnTest(x=transposed_taxa[,i], g=transposed_taxa$Description, p.adjust.method="fdr")
      write.table(z, file = "comparisons_compound_or_gene_kruskal.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
      write.table(krusk_res$method,file = "comparisons_compound_or_gene_kruskal.txt", append = TRUE, row.names = TRUE, col.names = TRUE, sep='\t')
      write.table(krusk_res$p.value,file = "comparisons_compound_or_gene_kruskal.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
      write.table(krusk_res$p.adjust.method,file = "comparisons_compound_or_gene_kruskal.txt", append = TRUE, row.names = TRUE, col.names = TRUE, sep='\t')
    }
  if (length(levels(map_groups$Description)) != 2){
    warning("You don't have two groups")} else if (x$p.value >= 0.05 && length(levels(transposed_taxa$Description)) == 2){
      ttest_res <- t.test(transposed_taxa[,i] ~ transposed_taxa$Description)
      data_info <- c(data_info, z)
      data_t <- c(data_t, ttest_res[1])
      data_fdr <- c(data_fdr, ttest_res[3])
      data_pvalue <- c(data_pvalue, ttest_res[["p.value"]] )
      data_ttest_sum <- cbind(data_info, data_t, data_fdr,data_pvalue)
      colnames(data_ttest_sum)[1:4] <- c('Data', 't_value', 'FDR','pvalue')
    } else if (x$p.value < 0.05 && length(levels(transposed_taxa$Description)) == 2){
      wilcox_res <-wilcox.test(transposed_taxa[,i] ~ transposed_taxa$Description)

      data_info_wil <- c(data_info_wil, z)
      data_wil_w <- c(data_wil_w, wilcox_res[1])
      data_wil_fdr <- c(data_wil_fdr, wilcox_res[3])
      data_wil_pvalue <- c(data_wil_pvalue, wilcox_res[["p.value"]] )
      data_wil_sum <- cbind(data_info_wil, data_wil_w, data_wil_fdr,data_wil_pvalue)
      colnames(data_wil_sum)[1:4] <- c('Data', 'w_value', 'FDR','pvalue')
    }
}
#
# if (length(levels(transposed_taxa$Description)) != 2){
#     warning("You don't have two groups")} else if (exists('data_ttest_sum') == TRUE){
#         data_ttest_sum[,3] <- p.adjust(data_ttest_sum[,3], method = 'fdr')
#         write.table(data_ttest_sum, file = "comparisons_compound_or_gene.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
#     }
# if (length(levels(transposed_taxa$Description)) != 2){
#     warning("You don't have two groups")} else if (exists('data_wil_sum') == TRUE){
#         data_wil_sum[,3] <- p.adjust(data_wil_sum[,3], method = 'fdr')
#         write.table(data_wil_sum, file = "comparisons_compound_or_gene.txt", append = TRUE, row.names = FALSE, col.names = TRUE, sep='\t')
#     }
#
#
#
# Create the boxplot for each taxa per comparison group
#
end_time <- Sys.time()
end_time - start_time
