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
# This is all you need to change:
work_directory <- "/Users/tiagofalcon/Desktop/socorro" # Path to your work directory. This is 
# the directory where you may find the two necessary files to begin the analysis and where you gonna save the outputs.
file_name <- "sim-alta-var.csv" # Change the file name. It is a tab delimited OTU table.
map_file_data <- "mapping_file.txt" # Mapping file with groups description names.
#
##Import the data from .rdata
input<-read.csv(file_name)


#
#
# Necessary packages
library("psych") # For the geometric mean calculation
library ("pvclust") # For the hierarchical clustering
library("gplots") # For the heat maps
library("PMCMRplus") # For non parametric pairwised tests
library("ggbiplot") # For PCA plot
library("vegan") # For PERMANOVA analysis and plot, and for diversity indexes
#
#
#create a transformed_df backup to modify the colnames by color names to use as side colors ...
#
#
#### Aitchison's log ratio transformation...
#
# Go to the work directory.
setwd(work_directory)
# Read the input file... In this case we used the OTU table with raw 
# reads count per sample from a 16S profilling analysis. The standard 
# file is the 16S BMP pipeline output OTU table. Could be either the TXT
# or SPF file.
input <- read.delim(file_name)
#
map_groups <- read.delim(map_file_data)
map_groups<-map_groups[,c("X.SampleID","Description")]
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


transposed_taxa <- t(end_file)
colnames(transposed_taxa) <- transposed_taxa[1,]
transposed_taxa <- transposed_taxa[-1,]
transposed_taxa <- as.data.frame(transposed_taxa)
# Make sure that values in data frame are numeric
for (i in 1:ncol(transposed_taxa)){
  transposed_taxa[,i] <- as.numeric(as.character(transposed_taxa[,i]))
}
# Change samples per groups names
merged_df <- merge(map_groups, transposed_taxa, by.x='X.SampleID', by.y= 'row.names')
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
if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (exists('data_ttest_sum') == TRUE){
        data_ttest_sum[,3] <- p.adjust(data_ttest_sum[,3], method = 'fdr')
        write.table(data_ttest_sum, file = "comparisons_compound_or_gene.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
    }
if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (exists('data_wil_sum') == TRUE){
        data_wil_sum[,3] <- p.adjust(data_wil_sum[,3], method = 'fdr')
        write.table(data_wil_sum, file = "comparisons_compound_or_gene.txt", append = TRUE, row.names = FALSE, col.names = TRUE, sep='\t')
    }
#
#
#
# Create the boxplot for each taxa per comparison group
for_boxplot<-cbind (combine__rows_compound_or_gene_m[,lapply(combine__rows_compound_or_gene_m, class) == "factor"], df_prop_otus)
transposed_for_boxplot <- t(for_boxplot)
colnames(transposed_for_boxplot) <- transposed_for_boxplot[1,]
transposed_for_boxplot <- transposed_for_boxplot[-1,]
transposed_for_boxplot <- as.data.frame(transposed_for_boxplot)
transposed_for_boxplot <- merge (map_groups, transposed_for_boxplot,by.x = "X.SampleID", by.y="row.names")
transposed_for_boxplot <- transposed_for_boxplot[,c(2:ncol(transposed_for_boxplot))]
for (i in 2:ncol(transposed_for_boxplot)){
  name = paste("boxplot_compound_or_gene",i-1,".pdf",sep="")
  pdf(name, width=10, height=4)
  boxplot(as.numeric(as.character(transposed_for_boxplot[,i])) ~ Description, data = transposed_for_boxplot, main = colnames(transposed_for_boxplot[i]), cex.axis=.5, ylab = "Proportion", xlab = "groups", col="gray")
  dev.off()
}
#
#
#
# Plot the heat map of taxa proportion per groups
#color pattern
transformed_df2 <- transformed_df
transformed_df2[nrow(transformed_df2)+1,] <- NA
for (i in levels(map_groups$Description)){
  data_i <- as.vector(subset(map_groups, Description == i)[,1])
  for (name_col in data_i) {
    transformed_df2[nrow(transformed_df2), which( colnames(transformed_df2)== name_col )] <- color_data[match(i,levels(map_groups$Description))]
  }
}
colors_in <-as.vector(as.character(transformed_df2[nrow(transformed_df2),]))
#
for (i in 2:ncol(transposed_for_boxplot)){
  transposed_for_boxplot[,i] <- as.numeric(as.character(transposed_for_boxplot[,i]))
}
for_samp_prop <- aggregate(transposed_for_boxplot[,-1], by = list(transposed_for_boxplot[,1]), FUN = sum)
row.names(for_samp_prop) <- for_samp_prop[,1]
for_samp_prop <- for_samp_prop[,-1]
for_samp_prop <- for_samp_prop/rowSums(for_samp_prop)
mat_samp <- as.matrix(for_samp_prop)
pdf('plot_heat_map_proportion_compound_or_gene_groups.pdf')
heatmap.2(mat_samp,dendrogram = c("both"),key.title = "Proportion",key.xlab = "Proportion",
          cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, main = "Proportion 16S Compound/gene between groups",
          notecol="black", density.info="none", trace="none")
dev.off()
# log 10 proportion plot
pdf('plot_heat_map_proportion_compound_or_gene_groups_log_scale.pdf')
heatmap.2(log10(mat_samp +1),dendrogram = c("both"),key.title = "Proportion log scale +1",
          key.xlab = "Proportion log scale",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, 
          main = "Proportion log scale 16S Compound/gene between groups",notecol="black", density.info="none", trace="none")
dev.off()
#
# Z-score heat map
transposed_for_boxplot_zs <- transposed_for_boxplot
for (column in 2:ncol(transposed_for_boxplot_zs)) {
  transposed_for_boxplot_zs[,column] <- (transposed_for_boxplot_zs[,column]-mean(transposed_for_boxplot_zs[,column]))/sd(transposed_for_boxplot_zs[,column])
}
for_samp_zscore <- aggregate(transposed_for_boxplot_zs[,-1], by = list(transposed_for_boxplot_zs[,1]), FUN = sum)
row.names(for_samp_zscore) <- for_samp_zscore[,1]
for_samp_zscore <- for_samp_zscore[,-1]
mat_samp_zscore <- as.matrix(for_samp_zscore)
pdf('plot_heat_map_zscore_compound_or_gene_groups.pdf')
heatmap.2(mat_samp_zscore,dendrogram = c("both"),key.title = "Z-score",
          key.xlab = "Z-score",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, 
          main = "Z-score compound/gene between groups",notecol="black", density.info="none", trace="none")
dev.off()
#
# Sample's correlation per Compound/gene (heat map correlation)
cor_compound_or_gene_between_samples <- cor(transformed_df)
pdf('plot_heat_map_correlation_between_samples_compound_or_gene.pdf')
heatmap.2(cor_compound_or_gene_between_samples,dendrogram = c("both"),key.title = "Correlation",
          key.xlab = "Correlation",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, 
          main = "Correlation 16S Compound/gene between samples",notecol="black", density.info="none", 
          trace="none", ColSideColors = colors_in, 
          RowSideColors = colors_in)
dev.off()
# Samples distances per Compound/gene (pvclust of euclidean distances)
result_compound_or_gene_pvclust <- pvclust(end_file_compound_or_gene[,-1], method.dist = "euclidean", method.hclust="complete", nboot=1000)
pdf('rplot_pvclust_compound_or_gene.pdf')
plot(result_compound_or_gene_pvclust, main="Compound/gene taxa plot", cex = .5, cex.pv=.3)
dev.off()
# Compound/gene's association (Perform all comparisons Compound/gene x Compound/gene and save association tables' results)
# Rename columns in sense to make file understandable by lmGC function
colnames(transposed_taxa) <- sub("p__", "", colnames(transposed_taxa))
colnames(transposed_taxa) <- sub("\\[", "", colnames(transposed_taxa))
colnames(transposed_taxa) <- sub("\\]", "", colnames(transposed_taxa))
#for (i in 2:ncol(transposed_taxa)){
#  for (j in 2:ncol(transposed_taxa)){
#    cname <- c(colnames(transposed_taxa[i]),colnames(transposed_taxa[j]))
#    write.table(cname, file = "regression_compound_or_gene.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
#    reg_value <- summary(lm(transposed_taxa[,i]~ transposed_taxa[,j] ,data=transposed_taxa))
#    write.table(reg_value$coefficients,file = "regression_compound_or_gene.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
#    write.table(reg_value[9],file = "regression_compound_or_gene.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
#    write.table(reg_value[10],file = "regression_compound_or_gene.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
#  }
#}
#
#
#
# PCA plot and contribution analysis
# Create a data frame with the variavles (taxa names) as columns and the groups (samples) in lines
mat_pca<-merged_df[,c(3:ncol(merged_df))]
# PCA calculation
comp_Data<-prcomp(mat_pca)
# The PCA plot
pca_plot <- ggbiplot( comp_Data, groups = merged_df$Description,  obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = TRUE)
pca_plot <- pca_plot + scale_color_discrete(name = '')
pca_plot <- pca_plot + theme(legend.direction = 'horizontal', legend.position = 'top')
pdf('pca_plot_compound_or_gene.pdf')
pca_plot
dev.off()
# Now, calculate each taxa contribution to samples' dispersion
# Use the compa_Data variable created earlier
pca_data_scores<-comp_Data$rotation
summary_pca_data<-summary(comp_Data)
loadings <- pca_data_scores^2
# The line 2 of the summary_pca_data is the exactly value of each component
summary_pca_data
# Pass the table of the summary_pca_data to another variable:
pca_descriptive<- summary_pca_data$importance
# Use the line 2 withe the valiues of the components
pca_comp_value<- pca_descriptive [2,]
# Calculate taxa contributions
taxa_contribution<-loadings %*% pca_comp_value # use the symbol "%" to separate the multiplication signal of from the names of matices and vectors
# See each taxon contribution
taxa_contribution
# Check of the sum is 1 (or really close to it since the numbers approximation)
sum (taxa_contribution)
# Compare the value and the visual of the PCA plot. The largest contribution values should be represented as the longer arraows in the plot.
# Save the taxa contributions to a file
write.table(taxa_contribution, "taxa_contribution_compound_or_gene.txt",sep="\t", col.names = NA)
#
#
#
# Check groups heterogeneity based on samples dispersion through a Permutational Multivariate Analysis of Variance (PERMANOVA)
# Create a matrix and use groups names as the line names
mat_pca <- as.matrix(mat_pca)
row.names (mat_pca) <- merged_df$Description
# Calculate samples' distances and ANOVA
samp_dis<-vegdist(mat_pca,method="euclidean")
mod_disper<-betadisper(samp_dis, rownames(mat_pca))
anova(mod_disper)
mult_tukey_res <- TukeyHSD(mod_disper)
write.table(mult_tukey_res$group, "PERMANOVA_mult_tukey_res_compound_or_gene.txt", sep="\t", row.names = TRUE, col.names = NA)
# Now the permutation comparison to check groups homogeneity
per_compound_or_gene<-permutest(mod_disper,pairwise=TRUE, control=permControl(nperm=999))
per_compound_or_gene
write.table(per_compound_or_gene$pairwise$observed, "PERMANOVA_observed_permuted_pvalue_compound_or_gene.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
write.table(per_compound_or_gene$pairwise$permuted, "PERMANOVA_observed_permuted_pvalue_compound_or_gene.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
write.table(per_compound_or_gene$statistic, "PERMANOVA_observed_permuted_pvalue_compound_or_gene.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
# Calculate the general p-value
res_adonis<-adonis(mat_pca~merged_df$Description, permutations=999, method="euclidean")
write.table(res_adonis$aov.tab, "PERMANOVA_adonis_res_compound_or_gene.txt", sep="\t", row.names = TRUE, col.names = NA)
pdf("permanova_homogeneity_compound_or_gene.pdf")
plot(mod_disper, xlab = paste( "PCoA1 explains ", round((mod_disper$eig[1]/sum(mod_disper$eig)*100), digits = 5), "% of sample's dispersion", sep=""), ylab = paste( "PCoA2 explains ", round((mod_disper$eig[2]/sum(mod_disper$eig)*100), digits = 5), "% of sample's dispersion", sep=""), main = "Compound or gene dispersion plot")
dev.off()
save.image("data_results.RData")
#
end_time <- Sys.time()
end_time - start_time
