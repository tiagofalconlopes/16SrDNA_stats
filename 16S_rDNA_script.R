##################################################################################################################################
#############################################16S rRNA profilling downstream analysis##############################################
##################################################################################################################################
#
start_time <- Sys.time()
#
# This pipeline assumes that your data has two or more groups of samples!!!
#
# This is all you need to change:
work_directory <- "/PATH/TO/WORK/DIRECTORY/" # Path to your work directory. This is 
# the directory where you may find the two necessary files to begin the analysis and where you gonna save the outputs.
file_name <- "otu_table_tax.spf" # Change the file name. It is a tab delimited OTU table.
map_file_data <- "mapping_file.txt" # Mapping file with groups description names.
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
#
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
# File for the diversity steps...
for_diversity <- input[,!lapply(input, class) == "factor"]
# Aitchison...
# Replace all zeros by a really small value (0.0000001). This step is 
# necessary once Aitchinson's log ratio transformation do not accept zeros.
# This value is kind of arbritary. As there are many cells with "zeros"
# and they are not fixed to specific compositions (variables), we opted
# to use a very low value in place of the zero.
input[input == 0] <- 0.0000001
# Retrieve only the values
numbers_to_norm <- input[,!lapply(input, class) == "factor"]
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
end_file <- cbind (input[,lapply(input, class) == "factor"], transformed_df)
# Save the transformed output file in case you need it later. "log_ratio_transformed_data.txt"
# contains the OTU table with the Aitchinson's log ratio transformation.
write.table (end_file, 'log_ratio_transformed_data.txt', sep = '\t', row.names = FALSE)
#
#
#
#
### Diversity indexes
# Shannon - results in the output file "shannon_diversity.txt".
t_for_diversity <- t(for_diversity)
shannon_div <- diversity(t_for_diversity, index = "shannon")
shannon_div <- data.frame(shannon_div)
shannon_div_fix_col <- cbind(row.names(shannon_div), shannon_div$shannon_div)
colnames(shannon_div_fix_col) <- c("Samples", "Shannon_diversity_index")
write.table(shannon_div_fix_col, "shannon_diversity.txt", sep="\t", row.names = FALSE)
#
# Simpson - results in the output file "simpson_diversity.txt"
simpson_div <- diversity(t_for_diversity, index = "simpson")
simpson_div <- data.frame(simpson_div)
simpson_div_fix_col <- cbind(row.names(simpson_div), simpson_div$simpson_div)
colnames(simpson_div_fix_col) <- c("Samples", "Simpson_diversity_index")
write.table(simpson_div_fix_col, "simpson_diversity.txt", sep="\t", row.names = FALSE)
#
div_merg <- merge (shannon_div, simpson_div, by = "row.names") 
map_groups <- read.delim(map_file_data)
div_map_merg <- merge (div_merg, map_groups, by.x = "Row.names", by.y = "X.SampleID")
# Shannon comparisons
sh_sp_test <- shapiro.test(div_map_merg$shannon_div)
if (length(levels(div_map_merg$Description)) <= 2){
  warning("You don't have more than two groups")} else if (sh_sp_test$p.value >= 0.05 && length(levels(div_map_merg$Description)) >= 3){
    sh_tuk_res <- TukeyHSD(aov(div_map_merg$shannon_div ~ div_map_merg$Description))
    sh_tuk_res$`div_map_merg$Description`[,4] <- p.adjust(sh_tuk_res$`div_map_merg$Description`[,4], method = 'fdr')
    write.table(sh_tuk_res$`div_map_merg$Description`, file = "comparisons_species_shannon_tukey.txt", row.names = TRUE, col.names = NA, sep='\t')
} else if (sh_sp_test$p.value < 0.05 && length(levels(div_map_merg$Description)) >= 3){
    sh_krusk_res <- kwAllPairsDunnTest(x=div_map_merg$shannon_div, g=div_map_merg$Description, p.adjust.method="fdr")
    write.table(sh_krusk_res$method,file = "comparisons_species_shannon_kruskal.txt", append = TRUE, row.names = TRUE, col.names = TRUE, sep='\t')
    write.table(sh_krusk_res$p.value,file = "comparisons_species_shannon_kruskal.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
}
if (length(levels(div_map_merg$Description)) != 2){
  warning("You don't have two groups")} else if (sh_sp_test$p.value >= 0.05 && length(levels(div_map_merg$Description)) == 2){
    ttest_res <- t.test(div_map_merg$shannon_div ~ div_map_merg$Description)
    sh_t <- c(ttest_res[1])
    sh_fdr <- c(ttest_res[3])
    sh_ttest_sum <- cbind(sh_t, sh_fdr)
    colnames(sh_ttest_sum)[1:2] <- c('t_value', 'FDR')
      } else if (sh_sp_test$p.value < 0.05 && length(levels(div_map_merg$Description)) == 2){
    wilcox_res <-wilcox.test(div_map_merg$shannon_div ~ div_map_merg$Description)
    sh_w <- c(wilcox_res[1])
    sh_w_fdr <- c(wilcox_res[3])
    sh_w_sum <- cbind(sh_w, sh_w_fdr)
    colnames(sh_w_sum)[1:2] <- c('w_value', 'FDR')
}
#
if (length(levels(div_map_merg$Description)) != 2){
  warning("You don't have two groups")} else if (exists('sh_ttest_sum') == TRUE){
    sh_ttest_sum[,2] <- p.adjust(sh_ttest_sum[,2], method = 'fdr')
    write.table(sh_ttest_sum, file = "comparisons_shannon_t_test.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
      } else if (exists('sh_w_sum') == TRUE){
    sh_w_sum[,2] <- p.adjust(sh_w_sum[,2], method = 'fdr')
    write.table(sh_w_sum, file = "comparisons_shannon_wilcox.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
}
#
# Simpson comparisons
si_sp_test <- shapiro.test(div_map_merg$simpson_div)
if (length(levels(div_map_merg$Description)) <= 2){
  warning("You don't have more than two groups")} else if (si_sp_test$p.value >= 0.05 && length(levels(div_map_merg$Description)) >= 3){
    si_tuk_res <- TukeyHSD(aov(div_map_merg$simpson_div ~ div_map_merg$Description))
    si_tuk_res$`div_map_merg$Description`[,4] <- p.adjust(si_tuk_res$`div_map_merg$Description`[,4], method = 'fdr')
    write.table(si_tuk_res$`div_map_merg$Description`, file = "comparisons_species_simpson_tukey.txt", row.names = TRUE, col.names = NA, sep='\t')
} else if (si_sp_test$p.value < 0.05 && length(levels(div_map_merg$Description)) >= 3){
    si_krusk_res <- kwAllPairsDunnTest(x=div_map_merg$simpson_div, g=div_map_merg$Description, p.adjust.method="fdr")
    write.table(si_krusk_res$method,file = "comparisons_species_simpson_kruskal.txt", append = TRUE, row.names = TRUE, col.names = TRUE, sep='\t')
    write.table(si_krusk_res$p.value,file = "comparisons_species_simpson_kruskal.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
}
if (length(levels(div_map_merg$Description)) != 2){
  warning("You don't have two groups")} else if (si_sp_test$p.value >= 0.05 && length(levels(div_map_merg$Description)) == 2){
    ttest_res <- t.test(div_map_merg$simpson_div ~ div_map_merg$Description)
    si_t <- c(ttest_res[1])
    si_fdr <- c(ttest_res[3])
    si_ttest_sum <- cbind(si_t, si_fdr)
    colnames(si_ttest_sum)[1:2] <- c('t_value', 'FDR')
      } else if (si_sp_test$p.value < 0.05 && length(levels(div_map_merg$Description)) == 2){
    wilcox_res <-wilcox.test(div_map_merg$simpson_div ~ div_map_merg$Description)
    si_w <- c(wilcox_res[1])
    si_w_fdr <- c(wilcox_res[3])
    si_w_sum <- cbind(si_w, si_w_fdr)
    colnames(si_w_sum)[1:2] <- c('w_value', 'FDR')
}
#
if (length(levels(div_map_merg$Description)) != 2){
  warning("You don't have two groups")} else if (exists('si_ttest_sum') == TRUE){
    si_ttest_sum[,2] <- p.adjust(si_ttest_sum[,2], method = 'fdr')
    write.table(si_ttest_sum, file = "comparisons_simpson_t_test.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
      } else if (exists('si_w_sum') == TRUE){
    si_w_sum[,2] <- p.adjust(si_w_sum[,2], method = 'fdr')
    write.table(si_w_sum, file = "comparisons_simpson_wilcox.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
}
#
#
# Species richness using the logarithm base in Shannon - output in "species_richness_groups.txt"
merge_div <- merge(t_for_diversity, map_groups, by.x = "row.names", by.y = "X.SampleID")
merge_div_clean <- merge_div[,2:(ncol(t_for_diversity)+1)]
richness_div_groups <- specnumber(merge_div_clean, groups = merge_div$Description)
df_richness_div_groups <- data.frame(richness_div_groups)
df_richness_div_groups <- cbind(row.names(df_richness_div_groups), df_richness_div_groups$richness_div_groups)
colnames(df_richness_div_groups) <- c("Samples", "richness")
write.table(df_richness_div_groups, "species_richness_groups.txt", sep="\t", row.names = FALSE)
# Richness comparisons
richness_div_samp <- specnumber(merge_div_clean, groups = merge_div$Row.names)
df_div_dif <- data.frame(richness_div_samp, merge_div$Row.names, merge_div$Description)
sp_test <- shapiro.test(df_div_dif$richness_div_samp)
if (length(levels(df_div_dif$merge_div.Description)) <= 2){
  warning("You don't have three or more groups")} else if (sp_test$p.value >= 0.05 && length(levels(df_div_dif$merge_div.Description)) >= 3){
    tuk_res <- TukeyHSD(aov(df_div_dif$richness_div_samp ~ df_div_dif$merge_div.Description))
    tuk_res$`df_div_dif$merge_div.Description`[,4] <- p.adjust(tuk_res$`df_div_dif$merge_div.Description`[,4], method = 'fdr')
    write.table(tuk_res$`df_div_dif$merge_div.Description`, file = "comparisons_species_richness_tukey.txt", row.names = TRUE, col.names = NA, sep='\t')
  } else if (sp_test$p.value < 0.05 && length(levels(df_div_dif$merge_div.Description)) >= 3){
    krusk_res <- kwAllPairsDunnTest(x=df_div_dif$richness_div_samp, g=df_div_dif$merge_div.Description, p.adjust.method="fdr")
    write.table(krusk_res$method,file = "comparisons_species_richness_kruskal.txt", append = TRUE, row.names = TRUE, col.names = TRUE, sep='\t')
    write.table(krusk_res$p.value,file = "comparisons_species_richness_kruskal.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
}
if (length(levels(df_div_dif$merge_div.Description)) != 2){
  warning("You don't have two groups")} else if (sp_test$p.value >= 0.05 && length(levels(df_div_dif$merge_div.Description)) == 2){
    ttest_res <- t.test(df_div_dif$richness_div_samp ~ df_div_dif$merge_div.Description)
    sr_t <- c(ttest_res[1])
    sr_fdr <- c(ttest_res[3])
    sr_ttest_sum <- cbind(sr_t, sr_fdr)
    colnames(sr_ttest_sum)[1:2] <- c('t_value', 'FDR')
  } else if (sp_test$p.value < 0.05 && length(levels(df_div_dif$merge_div.Description)) == 2){
    wilcox_res <-wilcox.test(df_div_dif$richness_div_samp ~ df_div_dif$merge_div.Description)
    sr_w <- c(wilcox_res[1])
    sr_w_fdr <- c(wilcox_res[3])
    sr_w_sum <- cbind(sr_w, sr_w_fdr)
    colnames(sr_w_sum)[1:2] <- c('w_value', 'FDR')
}
#
if (length(levels(df_div_dif$merge_div.Description)) != 2){
  warning("You don't have two groups")} else if (exists('sr_ttest_sum') == TRUE){
    sr_ttest_sum[,2] <- p.adjust(sr_ttest_sum[,2], method = 'fdr')
    write.table(sr_ttest_sum, file = "comparisons_species_richness_t_test.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
  } else if (exists('sr_w_sum') == TRUE){
    sr_w_sum[,2] <- p.adjust(sr_w_sum[,2], method = 'fdr')
    write.table(sr_w_sum, file = "comparisons_species_richness_wilcox.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
}
#
#
# Plot richness
pdf("richness_plot.pdf", width=10, height=4)
boxplot(as.numeric(as.character(df_div_dif$richness_div_samp)) ~ merge_div.Description, data = df_div_dif, main = "Species richness Plot", cex.axis=.5, ylab = "Richness", xlab = "groups", col="gray")
dev.off()
# Rarefaction values and plot
sp_raref <- rarefy(t_for_diversity, min(rowSums(t_for_diversity)))
write.table(sp_raref,file = "rarefaction.txt", row.names = TRUE, col.names = TRUE, sep = '\t')
pdf("rarefaction_plot_per_sample.pdf")
rarecurve(t_for_diversity, xlab = "Reads counts per sample", ylab = "Number of species", label = TRUE)
dev.off()
#
# Beta diversity
new_dat <- t_for_diversity
#Calculate betadiversity
bdiv_res<-bdiv_res <- betadiver(new_dat, "z") # see "betadiver(help=TRUE)" for options
# Calculate the dispersion based on betadiversity
bdiv_disper <- with(merge_div, betadisper(bdiv_res, Description))
pdf("betadiversity.pdf")
plot(bdiv_disper, xlab = paste( "PCoA1 explains ", round((bdiv_disper$eig[1]/sum(bdiv_disper$eig)*100), digits = 5), "% of sample's dispersion", sep=""), ylab = paste( "PCoA2 explains ", round((bdiv_disper$eig[2]/sum(bdiv_disper$eig)*100), digits = 5), "% of sample's dispersion", sep=""), main = "Beta diversity dispersion plot")
dev.off()
#
#
#
### Statistical analysis over transformed data ...
#
# 
#
# Get the hierarchical clusterization based on all taxa based on the Euclidean distances ...
result_all_taxa_pvclust <- pvclust(transformed_df, method.dist = "euclidean", method.hclust="complete", nboot=1000)
pdf('rplot_pvclust_all_taxa.pdf')
plot(result_all_taxa_pvclust, main="All taxa plot", cex = .5, cex.pv=.3)
dev.off()
# ... and now the correlation between samples heat map
mypalette<-colorRampPalette(c("red","yellow", "green"))(n=50)
color_data <- c("blue", "cyan", "lightgreen",
                "yellow", "gray", "gold",
                "green", "lightblue", "pink",
                "orange", "darkblue", "violet",
                "red", "brown", "darkgreen") # for 15 groups maximum - change the order and/or add colors in your preference.
#create a transformed_df backup to modify the colnames by color names to use as side colors ...
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
cor_all_taxa_between_samples <- cor(transformed_df)
pdf('plot_heat_map_correlation_between_samples_all_taxa.pdf')
heatmap.2(cor_all_taxa_between_samples, dendrogram = c("both"), key.title = "Correlation", 
          key.xlab = "Correlation", cexRow=0.3, Colv="Rowv", cexCol=0.3, keysize = 1, col= mypalette, 
          main = "Correlation 16S all taxa between samples", notecol="black", density.info="none", 
          trace="none", ColSideColors = colors_in, 
          RowSideColors = colors_in)
dev.off()
#
#
#
# Now let's evaluate the different taxonomic levels. For this purpose, we gonna 
# exclude unclassified taxa. Aitchison's transformation should be done again.
# To change to others taxonomic leves, modify the "input$level_N" column and
# the name of the taxon level.
#
# First, for the ###kingdom### level
if (length(levels(input$Level_1)) == 2){
  kingdom_no_unclassified <- subset(input,!(input$Level_1=="Unclassified"))
  kingdom_no_unclassified_edt <-cbind (kingdom_no_unclassified$Level_1, kingdom_no_unclassified[,!lapply(kingdom_no_unclassified, class) == "factor"])
  # One row per taxon...
  combine__rows_kingdom_m <- aggregate(kingdom_no_unclassified_edt[,-1], by = list(kingdom_no_unclassified_edt[,1]), FUN = sum)
  # Prepare to and transform data again... We are using the same variable names (overwriting) except for the final transformed file
  # Retrieve only the values
  numbers_to_norm <- combine__rows_kingdom_m[,!lapply(combine__rows_kingdom_m, class) == "factor"]
  # Transform each column to proportion.
  prop_otus<-t(t(numbers_to_norm)/rowSums(t(numbers_to_norm)))
  # Check if each column's sum is 1.
  colSums(prop_otus)
  # Convert prop_otus to data frame format.
  df_prop_otus<-as.data.frame(prop_otus)
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
  end_file_kingdom <- cbind (combine__rows_kingdom_m[,lapply(combine__rows_kingdom_m, class) == "factor"], transformed_df)
  # Save the transformed output file in case you need it later.
  write.table (end_file_kingdom, 'kingdom_log_ratio_transformed_data.txt', sep = '\t', row.names = FALSE)
  #
  #
  #
  # Create a barplot with top five taxa based on the first group ...
  combine__rows_kingdom_m2 <- t(combine__rows_kingdom_m)
  colnames(combine__rows_kingdom_m2)<- combine__rows_kingdom_m2[1,]
  combine__rows_kingdom_m2 <- data.frame(combine__rows_kingdom_m2)
  combine__rows_kingdom_m2 <- combine__rows_kingdom_m2[-1,]
  combine__rows_kingdom_m2[, ncol(combine__rows_kingdom_m2)+1] <- NA
  for (i in levels(map_groups$Description)){
    data_i <- as.vector(subset(map_groups, Description == i)[,1])
    for (name_col in data_i) {
      combine__rows_kingdom_m2[which(row.names(combine__rows_kingdom_m2) == name_col), ncol(combine__rows_kingdom_m2)] <- i
    }
  }
  colnames(combine__rows_kingdom_m2)[ncol(combine__rows_kingdom_m2)] <- c("group_ID")
  for (i in 1:(ncol(combine__rows_kingdom_m2)-1)){
    combine__rows_kingdom_m2[,i] <- as.numeric(as.character(combine__rows_kingdom_m2[,i]))
  }
  combine__rows_kingdom_m2 <- aggregate(combine__rows_kingdom_m2[,-ncol(combine__rows_kingdom_m2)], by = list(combine__rows_kingdom_m2[,ncol(combine__rows_kingdom_m2)]), FUN = mean)
  t_combine__rows_kingdom_m2 <- t(combine__rows_kingdom_m2)
  colnames(t_combine__rows_kingdom_m2) <- t_combine__rows_kingdom_m2[1,]
  t_combine__rows_kingdom_m2 <- t_combine__rows_kingdom_m2[-1,]
  t_combine__rows_kingdom_m2 <- data.frame(t_combine__rows_kingdom_m2)
  for (i in 1:ncol(t_combine__rows_kingdom_m2)){
    t_combine__rows_kingdom_m2[,i] <- as.numeric(as.character(t_combine__rows_kingdom_m2[,i]))
  }
  prop_kingdom_data <- t(t(t_combine__rows_kingdom_m2)/rowSums(t(t_combine__rows_kingdom_m2)))
  colSums(prop_kingdom_data)
  prop_kingdom_data <- prop_kingdom_data[order(prop_kingdom_data[,1], decreasing = TRUE),]
  pdf('barplot_top2_kingdom.pdf')
  par(xpd = NA)
  barplot(prop_kingdom_data[1:2,],main = "Top 5 most representative kingdom",
          xlab = "Groups", cex.names = .3,  ylab = "Proportion", yaxt="n",
          col = color_data[1:5], ylim=c(0,1.2))
  ticks<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
  axis(2,at=ticks,labels=ticks)
  legend("topright", row.names(prop_kingdom_data[1:2,]), cex=0.5, fill = color_data[1:2])
  dev.off()
  #
  #
  #
  # Normality tests and comparisons
  # Transpose the end_file_taxa Level table and create the transposed data frame
  transposed_taxa <- t(end_file_kingdom)
  colnames(transposed_taxa) <- transposed_taxa[1,]
  transposed_taxa <- transposed_taxa[-1,]
  transposed_taxa <- as.data.frame(transposed_taxa)
  # Make sure that values in data frame are numeric
  for (i in 1:ncol(transposed_taxa)){
    transposed_taxa[,i] <- as.numeric(as.character(transposed_taxa[,i]))
  }
  # Change samples per groups names
  merged_df <- merge(map_groups, transposed_taxa, by.x='X.SampleID', by.y= 'row.names')
  transposed_taxa <- merged_df[,c(4:ncol(merged_df))]
  # Perform normality test and the adequate comparisons
  kingdom_taxa <- c()
  kingdom_t <- c()
  kingdom_fdr <- c()
  kingdom_taxa_wil <- c()
  kingdom_wil_w <- c()
  kingdom_wil_fdr <- c()
  for (i in 2:ncol(transposed_taxa)){
    z<-colnames(transposed_taxa[i])
    x <- shapiro.test(transposed_taxa[,i])
    if (length(levels(transposed_taxa$Description)) <= 2){
      warning("You don't have more than two groups")} else if (x$p.value >= 0.05 && length(levels(transposed_taxa$Description)) >= 3){
        tuk_res <- TukeyHSD(aov(transposed_taxa[,i] ~ transposed_taxa$Description))
        write.table(z, file = "comparisons_kingdom_tukey.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
        tuk_res$`transposed_taxa$Description`[,4] <- p.adjust(tuk_res$`transposed_taxa$Description`[,4],method = 'fdr')
        write.table(tuk_res$`transposed_taxa$Description`, file = "comparisons_kingdom_tukey.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
      } else if (x$p.value < 0.05 && length(levels(transposed_taxa$Description)) >= 3){
        krusk_res <- kwAllPairsDunnTest(x=transposed_taxa[,i], g=transposed_taxa$Description, p.adjust.method="fdr")
        write.table(z, file = "comparisons_kingdom_kruskal.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
        write.table(krusk_res$method,file = "comparisons_kingdom_kruskal.txt", append = TRUE, row.names = TRUE, col.names = TRUE, sep='\t')
        write.table(krusk_res$p.value,file = "comparisons_kingdom_kruskal.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
      }
    if (length(levels(transposed_taxa$Description)) != 2){
      warning("You don't have two groups")} else if (x$p.value >= 0.05 && length(levels(transposed_taxa$Description)) == 2){
        ttest_res <- t.test(transposed_taxa[,i] ~ transposed_taxa$Description)
        kingdom_taxa <- c(kingdom_taxa, z)
        kingdom_t <- c(kingdom_t, ttest_res[1])
        kingdom_fdr <- c(kingdom_fdr, ttest_res[3])
        kingdom_ttest_sum <- cbind(kingdom_taxa, kingdom_t, kingdom_fdr)
        colnames(kingdom_ttest_sum)[1:3] <- c('kingdom', 't_value', 'FDR')
      } else if (x$p.value < 0.05 && length(levels(transposed_taxa$Description)) == 2){
        wilcox_res <-wilcox.test(transposed_taxa[,i] ~ transposed_taxa$Description)
        kingdom_taxa_wil <- c(kingdom_taxa_wil, z)
        kingdom_wil_w <- c(kingdom_wil_w, wilcox_res[1])
        kingdom_wil_fdr <- c(kingdom_wil_fdr, wilcox_res[3])
        kingdom_wil_sum <- cbind(kingdom_taxa_wil, kingdom_wil_w, kingdom_wil_fdr)
        colnames(kingdom_wil_sum)[1:3] <- c('kingdom', 'w_value', 'FDR')
      }
  }
} else if (length(levels(input$Level_1)) == 1) { print ("You have only one dependent variable at Level_1")} else { print ("You have some problems!!!")}
#
if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (exists('kingdom_ttest_sum') == TRUE){
        kingdom_ttest_sum[,3] <- p.adjust(kingdom_ttest_sum[,3], method = 'fdr')
        write.table(kingdom_ttest_sum, file = "comparisons_kingdom_t_test.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
    }
if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (exists('kingdom_wil_sum') == TRUE){
        kingdom_wil_sum[,3] <- p.adjust(kingdom_wil_sum[,3], method = 'fdr')
        write.table(kingdom_wil_sum, file = "comparisons_kingdom_wilcox.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
    }
#
#
####Phylum### level
phylum_no_unclassified <- subset(input,!(input$Level_2=="Unclassified"))
phylum_no_unclassified_edt <-cbind (phylum_no_unclassified$Level_2, phylum_no_unclassified[,!lapply(phylum_no_unclassified, class) == "factor"])
# One row per taxon...
combine__rows_phylum_m <- aggregate(phylum_no_unclassified_edt[,-1], by = list(phylum_no_unclassified_edt[,1]), FUN = sum)
# Prepare to and transform data again... We are using the same variable names (overwriting) except for the final transformed file
# Retrieve only the values
numbers_to_norm <- combine__rows_phylum_m[,!lapply(combine__rows_phylum_m, class) == "factor"]
# Transform each column to proportion.
prop_otus<-t(t(numbers_to_norm)/rowSums(t(numbers_to_norm)))
# Check if each column's sum is 1.
colSums(prop_otus)
# Convert prop_otus to data frame format.
df_prop_otus<-as.data.frame(prop_otus)
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
end_file_phylum <- cbind (combine__rows_phylum_m[,lapply(combine__rows_phylum_m, class) == "factor"], transformed_df)
# Save the transformed output file in case you need it later.
write.table (end_file, 'phylum_log_ratio_transformed_data.txt', sep = '\t', row.names = FALSE)
#
#
#
# Create a barplot with top five taxa based on the first group ...
combine__rows_phylum_m2 <- t(combine__rows_phylum_m)
colnames(combine__rows_phylum_m2)<- combine__rows_phylum_m2[1,]
combine__rows_phylum_m2 <- data.frame(combine__rows_phylum_m2)
combine__rows_phylum_m2 <- combine__rows_phylum_m2[-1,]
combine__rows_phylum_m2[, ncol(combine__rows_phylum_m2)+1] <- NA
for (i in levels(map_groups$Description)){
  data_i <- as.vector(subset(map_groups, Description == i)[,1])
  for (name_col in data_i) {
    combine__rows_phylum_m2[which(row.names(combine__rows_phylum_m2) == name_col), ncol(combine__rows_phylum_m2)] <- i
  }
}
colnames(combine__rows_phylum_m2)[ncol(combine__rows_phylum_m2)] <- c("group_ID")
for (i in 1:(ncol(combine__rows_phylum_m2)-1)){
  combine__rows_phylum_m2[,i] <- as.numeric(as.character(combine__rows_phylum_m2[,i]))
}
combine__rows_phylum_m2 <- aggregate(combine__rows_phylum_m2[,-ncol(combine__rows_phylum_m2)], by = list(combine__rows_phylum_m2[,ncol(combine__rows_phylum_m2)]), FUN = mean)
t_combine__rows_phylum_m2 <- t(combine__rows_phylum_m2)
colnames(t_combine__rows_phylum_m2) <- t_combine__rows_phylum_m2[1,]
t_combine__rows_phylum_m2 <- t_combine__rows_phylum_m2[-1,]
t_combine__rows_phylum_m2 <- data.frame(t_combine__rows_phylum_m2)
for (i in 1:ncol(t_combine__rows_phylum_m2)){
  t_combine__rows_phylum_m2[,i] <- as.numeric(as.character(t_combine__rows_phylum_m2[,i]))
}
prop_phylum_data <- t(t(t_combine__rows_phylum_m2)/rowSums(t(t_combine__rows_phylum_m2)))
colSums(prop_phylum_data)
prop_phylum_data <- prop_phylum_data[order(prop_phylum_data[,1], decreasing = TRUE),]
pdf('barplot_top5_phylum.pdf')
par(xpd = NA)
barplot(prop_phylum_data[1:5,],main = "Top 5 most representative Phylum",
        xlab = "Groups", cex.names = .3,  ylab = "Proportion", yaxt="n",
        col = color_data[1:5], ylim=c(0,1.2))
ticks<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
axis(2,at=ticks,labels=ticks)
legend("topright", row.names(prop_phylum_data[1:5,]), cex=0.5, fill = color_data[1:5])
dev.off()
#
#
#
# Normality tests and comparisons
# Transpose the end_file_taxa Level table and create the transposed data frame
transposed_taxa <- t(end_file_phylum)
colnames(transposed_taxa) <- transposed_taxa[1,]
transposed_taxa <- transposed_taxa[-1,]
transposed_taxa <- as.data.frame(transposed_taxa)
# Make sure that values in data frame are numeric
for (i in 1:ncol(transposed_taxa)){
  transposed_taxa[,i] <- as.numeric(as.character(transposed_taxa[,i]))
}
# Change samples per groups names
merged_df <- merge(map_groups, transposed_taxa, by.x='X.SampleID', by.y= 'row.names')
transposed_taxa <- merged_df[,c(4:ncol(merged_df))]
# Perform normality test and the adequate comparisons
phylum_taxa <- c()
phylum_t <- c()
phylum_fdr <- c()
phylum_taxa_wil <- c()
phylum_wil_w <- c()
phylum_wil_fdr <- c()
for (i in 2:ncol(transposed_taxa)){
  z<-colnames(transposed_taxa[i])
  x <- shapiro.test(transposed_taxa[,i])
  if (length(levels(transposed_taxa$Description)) <= 2){
    warning("You don't have more than two groups")} else if (x$p.value >= 0.05 && length(levels(transposed_taxa$Description)) >= 3){
      tuk_res <- TukeyHSD(aov(transposed_taxa[,i] ~ transposed_taxa$Description))
      write.table(z, file = "comparisons_phylum_tukey.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
      tuk_res$`transposed_taxa$Description`[,4] <- p.adjust(tuk_res$`transposed_taxa$Description`[,4],method = 'fdr')
      write.table(tuk_res$`transposed_taxa$Description`, file = "comparisons_phylum_tukey.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    } else if (x$p.value < 0.05 && length(levels(transposed_taxa$Description)) >= 3){
      krusk_res <- kwAllPairsDunnTest(x=transposed_taxa[,i], g=transposed_taxa$Description, p.adjust.method="fdr")
      write.table(z, file = "comparisons_phylum_kruskal.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
      write.table(krusk_res$method,file = "comparisons_phylum_kruskal.txt", append = TRUE, row.names = TRUE, col.names = TRUE, sep='\t')
      write.table(krusk_res$p.value,file = "comparisons_phylum_kruskal.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    }
  if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (x$p.value >= 0.05 && length(levels(transposed_taxa$Description)) == 2){
      ttest_res <- t.test(transposed_taxa[,i] ~ transposed_taxa$Description)
      phylum_taxa <- c(phylum_taxa, z)
      phylum_t <- c(phylum_t, ttest_res[1])
      phylum_fdr <- c(phylum_fdr, ttest_res[3])
      phylum_ttest_sum <- cbind(phylum_taxa, phylum_t, phylum_fdr)
      colnames(phylum_ttest_sum)[1:3] <- c('Phylum', 't_value', 'FDR')
    } else if (x$p.value < 0.05 && length(levels(transposed_taxa$Description)) == 2){
      wilcox_res <-wilcox.test(transposed_taxa[,i] ~ transposed_taxa$Description)
      phylum_taxa_wil <- c(phylum_taxa_wil, z)
      phylum_wil_w <- c(phylum_wil_w, wilcox_res[1])
      phylum_wil_fdr <- c(phylum_wil_fdr, wilcox_res[3])
      phylum_wil_sum <- cbind(phylum_taxa_wil, phylum_wil_w, phylum_wil_fdr)
      colnames(phylum_wil_sum)[1:3] <- c('Phylum', 'w_value', 'FDR')
    }
}
#
#
if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (exists('phylum_ttest_sum') == TRUE){
        phylum_ttest_sum[,3] <- p.adjust(phylum_ttest_sum[,3], method = 'fdr')
        write.table(phylum_ttest_sum, file = "comparisons_phylum_t_test.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
    }
if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (exists('phylum_wil_sum') == TRUE){
        phylum_wil_sum[,3] <- p.adjust(phylum_wil_sum[,3], method = 'fdr')
        write.table(phylum_wil_sum, file = "comparisons_phylum_wilcox.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
    }
#
# Create the boxplot for each taxa per comparison group
for_boxplot<-cbind (combine__rows_phylum_m[,lapply(combine__rows_phylum_m, class) == "factor"], df_prop_otus)
transposed_for_boxplot <- t(for_boxplot)
colnames(transposed_for_boxplot) <- transposed_for_boxplot[1,]
transposed_for_boxplot <- transposed_for_boxplot[-1,]
transposed_for_boxplot <- as.data.frame(transposed_for_boxplot)
transposed_for_boxplot <- merge (map_groups, transposed_for_boxplot,by.x = "X.SampleID", by.y="row.names")
transposed_for_boxplot <- transposed_for_boxplot[,c(4,5:ncol(transposed_for_boxplot))]
for (i in 2:ncol(transposed_for_boxplot)){
  name = paste("boxplot_phylum_",i-1,".pdf",sep="")
  pdf(name, width=10, height=4)
  boxplot(as.numeric(as.character(transposed_for_boxplot[,i])) ~ Description, data = transposed_for_boxplot, main = colnames(transposed_for_boxplot[i]), cex.axis=.5, ylab = "Proportion", xlab = "groups", col="gray")
  dev.off()
}
#
#
#
# Plot the heat map of taxa proportion per groups
for (i in 2:ncol(transposed_for_boxplot)){
  transposed_for_boxplot[,i] <- as.numeric(as.character(transposed_for_boxplot[,i]))
}
for_samp_prop <- aggregate(transposed_for_boxplot[,-1], by = list(transposed_for_boxplot[,1]), FUN = sum)
row.names(for_samp_prop) <- for_samp_prop[,1]
for_samp_prop <- for_samp_prop[,-1]
for_samp_prop <- for_samp_prop/rowSums(for_samp_prop)
mat_samp <- as.matrix(for_samp_prop)
pdf('plot_heat_map_proportion_phylum_groups.pdf')
heatmap.2(mat_samp,dendrogram = c("both"),key.title = "Proportion",key.xlab = "Proportion",
          cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, main = "Proportion 16S Phylum between groups",
          notecol="black", density.info="none", trace="none")
dev.off()
# log 10 proportion plot
pdf('plot_heat_map_proportion_phylum_groups_log_scale.pdf')
heatmap.2(log10(mat_samp +1),dendrogram = c("both"),key.title = "Proportion log scale +1",
          key.xlab = "Proportion log scale",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, 
          main = "Proportion log scale 16S Phylum between groups",notecol="black", density.info="none", trace="none")
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
pdf('plot_heat_map_zscore_phylum_groups.pdf')
heatmap.2(mat_samp_zscore,dendrogram = c("both"),key.title = "Z-score",
          key.xlab = "Z-score",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, 
          main = "Z-score Phylum between groups",notecol="black", density.info="none", trace="none")
dev.off()
# Sample's correlation per Phylum (heat map correlation)
cor_phylum_between_samples <- cor(transformed_df)
pdf('plot_heat_map_correlation_between_samples_phylum.pdf')
heatmap.2(cor_phylum_between_samples,dendrogram = c("both"),key.title = "Correlation",
          key.xlab = "Correlation",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, 
          main = "Correlation 16S Phylum between samples",notecol="black", density.info="none", 
          trace="none", ColSideColors = colors_in, 
          RowSideColors = colors_in)
dev.off()
# Samples distances per Phylum (pvclust of euclidean distances)
result_phylum_pvclust <- pvclust(end_file_phylum[,-1], method.dist = "euclidean", method.hclust="complete", nboot=1000)
pdf('rplot_pvclust_phylum.pdf')
plot(result_phylum_pvclust, main="Phylum taxa plot", cex = .5, cex.pv=.3)
dev.off()
# Phylum's association (Perform all comparisons Phylum x Phylum and save association tables' results)
# Rename columns in sense to make file understandable by lmGC function
colnames(transposed_taxa) <- sub("p__", "", colnames(transposed_taxa))
colnames(transposed_taxa) <- sub("\\[", "", colnames(transposed_taxa))
colnames(transposed_taxa) <- sub("\\]", "", colnames(transposed_taxa))
for (i in 2:ncol(transposed_taxa)){
  for (j in 2:ncol(transposed_taxa)){
    cname <- c(colnames(transposed_taxa[i]),colnames(transposed_taxa[j]))
    write.table(cname, file = "regression_phylum.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
    reg_value <- summary(lm(transposed_taxa[,i]~ transposed_taxa[,j] ,data=transposed_taxa))
    write.table(reg_value$coefficients,file = "regression_phylum.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    write.table(reg_value[9],file = "regression_phylum.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    write.table(reg_value[10],file = "regression_phylum.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
  }
}
#
#
#
# PCA plot and contribution analysis
# Create a data frame with the variavles (taxa names) as columns and the groups (samples) in lines
mat_pca<-merged_df[,c(5:ncol(merged_df))]
# PCA calculation
comp_Data<-prcomp(mat_pca)
# The PCA plot
pca_plot <- ggbiplot( comp_Data, groups = merged_df$Description,  obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = TRUE)
pca_plot <- pca_plot + scale_color_discrete(name = '')
pca_plot <- pca_plot + theme(legend.direction = 'horizontal', legend.position = 'top')
pdf('pca_plot_phylum.pdf')
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
taxa_contribution<-loadings %*% pca_comp_value # use the symbol "%" to separate the multiplication signal from the names of matices and vectors
# See each taxon contribution
taxa_contribution
# Check of the sum is 1 (or really close to it since the numbers approximation)
sum (taxa_contribution)
# Compare the value and the visual of the PCA plot. The largest contribution values should be represented as the longer arraows in the plot.
# Save the taxa contributions to a file
write.table(taxa_contribution, "taxa_contribution_phylum.txt",sep="\t", col.names = NA)
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
write.table(mult_tukey_res$group, "PERMANOVA_mult_tukey_res_phylum.txt", sep="\t", row.names = TRUE, col.names = NA)
# Now the permutation comparison to check groups homogeneity
per_phylum<-permutest(mod_disper,pairwise=TRUE, control=permControl(nperm=999))
per_phylum
write.table(per_phylum$pairwise$observed, "PERMANOVA_observed_permuted_pvalue_phylum.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
write.table(per_phylum$pairwise$permuted, "PERMANOVA_observed_permuted_pvalue_phylum.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
write.table(per_phylum$statistic, "PERMANOVA_observed_permuted_pvalue_phylum.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
# Calculate the general p-value
res_adonis<-adonis(mat_pca~merged_df$Description, permutations=999, method="euclidean")
write.table(res_adonis$aov.tab, "PERMANOVA_adonis_res_phylum.txt", sep="\t", row.names = TRUE, col.names = NA)
pdf("permanova_homogeneity_phylum.pdf")
plot(mod_disper, xlab = paste( "PCoA1 explains ", round((mod_disper$eig[1]/sum(mod_disper$eig)*100), digits = 5), "% of sample's dispersion", sep=""), ylab = paste( "PCoA2 explains ", round((mod_disper$eig[2]/sum(mod_disper$eig)*100), digits = 5), "% of sample's dispersion", sep=""), main = "Phylum dispersion plot")
dev.off()
#
#
#
####Class###
class_no_unclassified <- subset(input,!(input$Level_3=="Unclassified"))
class_no_unclassified_edt <-cbind (class_no_unclassified$Level_3, class_no_unclassified[,!lapply(class_no_unclassified, class) == "factor"])
# One row per taxon...
combine__rows_class_m <- aggregate(class_no_unclassified_edt[,-1], by = list(class_no_unclassified_edt[,1]), FUN = sum)
# Prepare to and transform data again... We are using the same variable names (overwriting) except for the final transformed file
# Retrieve only the values
numbers_to_norm <- combine__rows_class_m[,!lapply(combine__rows_class_m, class) == "factor"]
# Transform each column to proportion.
prop_otus<-t(t(numbers_to_norm)/rowSums(t(numbers_to_norm)))
# Check if each column's sum is 1.
colSums(prop_otus)
# Convert prop_otus to data frame format.
df_prop_otus<-as.data.frame(prop_otus)
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
end_file_class <- cbind (combine__rows_class_m[,lapply(combine__rows_class_m, class) == "factor"], transformed_df)
# Save the transformed output file in case you need it later.
write.table (end_file, 'class_log_ratio_transformed_data.txt', sep = '\t', row.names = FALSE)
#
#
#
# Create a barplot with top five taxa based on the first group ...
combine__rows_class_m2 <- t(combine__rows_class_m)
colnames(combine__rows_class_m2)<- combine__rows_class_m2[1,]
combine__rows_class_m2 <- data.frame(combine__rows_class_m2)
combine__rows_class_m2 <- combine__rows_class_m2[-1,]
combine__rows_class_m2[, ncol(combine__rows_class_m2)+1] <- NA
for (i in levels(map_groups$Description)){
  data_i <- as.vector(subset(map_groups, Description == i)[,1])
  for (name_col in data_i) {
    combine__rows_class_m2[which(row.names(combine__rows_class_m2) == name_col), ncol(combine__rows_class_m2)] <- i
  }
}
colnames(combine__rows_class_m2)[ncol(combine__rows_class_m2)] <- c("group_ID")
for (i in 1:(ncol(combine__rows_class_m2)-1)){
  combine__rows_class_m2[,i] <- as.numeric(as.character(combine__rows_class_m2[,i]))
}
combine__rows_class_m2 <- aggregate(combine__rows_class_m2[,-ncol(combine__rows_class_m2)], by = list(combine__rows_class_m2[,ncol(combine__rows_class_m2)]), FUN = mean)
t_combine__rows_class_m2 <- t(combine__rows_class_m2)
colnames(t_combine__rows_class_m2) <- t_combine__rows_class_m2[1,]
t_combine__rows_class_m2 <- t_combine__rows_class_m2[-1,]
t_combine__rows_class_m2 <- data.frame(t_combine__rows_class_m2)
for (i in 1:ncol(t_combine__rows_class_m2)){
  t_combine__rows_class_m2[,i] <- as.numeric(as.character(t_combine__rows_class_m2[,i]))
}
prop_class_data <- t(t(t_combine__rows_class_m2)/rowSums(t(t_combine__rows_class_m2)))
colSums(prop_class_data)
prop_class_data <- prop_class_data[order(prop_class_data[,1], decreasing = TRUE),]
pdf('barplot_top5_class.pdf')
par(xpd = NA)
barplot(prop_class_data[1:5,],main = "Top 5 most representative Class",
        xlab = "Groups", cex.names = .3,  ylab = "Proportion", yaxt="n",
        col = color_data[1:5], ylim=c(0,1.2))
ticks<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
axis(2,at=ticks,labels=ticks)
legend("topright", row.names(prop_class_data[1:5,]), cex=0.5, fill = color_data[1:5])
dev.off()
#
#
#
# Normality tests and comparisons
# Transpose the end_file_taxa Level table and create the transposed data frame
transposed_taxa <- t(end_file_class)
colnames(transposed_taxa) <- transposed_taxa[1,]
transposed_taxa <- transposed_taxa[-1,]
transposed_taxa <- as.data.frame(transposed_taxa)
# Make sure that values in data frame are numeric
for (i in 1:ncol(transposed_taxa)){
  transposed_taxa[,i] <- as.numeric(as.character(transposed_taxa[,i]))
}
# Change samples per groups names
merged_df <- merge(map_groups, transposed_taxa, by.x='X.SampleID', by.y= 'row.names')
transposed_taxa <- merged_df[,c(4:ncol(merged_df))]
# Perform normality test and the adequate comparisons
class_taxa <- c()
class_t <- c()
class_fdr <- c()
class_taxa_wil <- c()
class_wil_w <- c()
class_wil_fdr <- c()
for (i in 2:ncol(transposed_taxa)){
  z<-colnames(transposed_taxa[i])
  x <- shapiro.test(transposed_taxa[,i])
  if (length(levels(transposed_taxa$Description)) <= 2){
    warning("You don't have more than two groups")} else if (x$p.value >= 0.05 && length(levels(transposed_taxa$Description)) >= 3){
      tuk_res <- TukeyHSD(aov(transposed_taxa[,i] ~ transposed_taxa$Description))
      write.table(z, file = "comparisons_class_tukey.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
      tuk_res$`transposed_taxa$Description`[,4] <- p.adjust(tuk_res$`transposed_taxa$Description`[,4],method = 'fdr')
      write.table(tuk_res$`transposed_taxa$Description`, file = "comparisons_class_tukey.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    } else if (x$p.value < 0.05 && length(levels(transposed_taxa$Description)) >= 3){
      krusk_res <- kwAllPairsDunnTest(x=transposed_taxa[,i], g=transposed_taxa$Description, p.adjust.method="fdr")
      write.table(z, file = "comparisons_class_kruskal.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
      write.table(krusk_res$method,file = "comparisons_class_kruskal.txt", append = TRUE, row.names = TRUE, col.names = TRUE, sep='\t')
      write.table(krusk_res$p.value,file = "comparisons_class_kruskal.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    }
  if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (x$p.value >= 0.05 && length(levels(transposed_taxa$Description)) == 2){
      ttest_res <- t.test(transposed_taxa[,i] ~ transposed_taxa$Description)
      class_taxa <- c(class_taxa, z)
      class_t <- c(class_t, ttest_res[1])
      class_fdr <- c(class_fdr, ttest_res[3])
      class_ttest_sum <- cbind(class_taxa, class_t, class_fdr)
      colnames(class_ttest_sum)[1:3] <- c('class', 't_value', 'FDR')
    } else if (x$p.value < 0.05 && length(levels(transposed_taxa$Description)) == 2){
      wilcox_res <-wilcox.test(transposed_taxa[,i] ~ transposed_taxa$Description)
      class_taxa_wil <- c(class_taxa_wil, z)
      class_wil_w <- c(class_wil_w, wilcox_res[1])
      class_wil_fdr <- c(class_wil_fdr, wilcox_res[3])
      class_wil_sum <- cbind(class_taxa_wil, class_wil_w, class_wil_fdr)
      colnames(class_wil_sum)[1:3] <- c('class', 'w_value', 'FDR')
    }
}
#
if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (exists('class_ttest_sum') == TRUE){
        class_ttest_sum[,3] <- p.adjust(class_ttest_sum[,3], method = 'fdr')
        write.table(class_ttest_sum, file = "comparisons_class_t_test.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
    }
if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (exists('class_wil_sum') == TRUE){
        class_wil_sum[,3] <- p.adjust(class_wil_sum[,3], method = 'fdr')
        write.table(class_wil_sum, file = "comparisons_class_wilcox.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
    }
#
#
# Create the boxplot for each taxa per comparison group
for_boxplot<-cbind (combine__rows_class_m[,lapply(combine__rows_class_m, class) == "factor"], df_prop_otus)
transposed_for_boxplot <- t(for_boxplot)
colnames(transposed_for_boxplot) <- transposed_for_boxplot[1,]
transposed_for_boxplot <- transposed_for_boxplot[-1,]
transposed_for_boxplot <- as.data.frame(transposed_for_boxplot)
transposed_for_boxplot <- merge (map_groups, transposed_for_boxplot,by.x = "X.SampleID", by.y="row.names")
transposed_for_boxplot <- transposed_for_boxplot[,c(4,5:ncol(transposed_for_boxplot))]
for (i in 2:ncol(transposed_for_boxplot)){
  name = paste("boxplot_class_",i-1,".pdf",sep="")
  pdf(name, width=10, height=4)
  boxplot(as.numeric(as.character(transposed_for_boxplot[,i])) ~ Description, data = transposed_for_boxplot, main = colnames(transposed_for_boxplot[i]), cex.axis=.5, ylab = "Proportion", xlab = "groups", col="gray")
  dev.off()
}
#
#
#
# Plot the heat map of taxa proportion per groups
for (i in 2:ncol(transposed_for_boxplot)){
  transposed_for_boxplot[,i] <- as.numeric(as.character(transposed_for_boxplot[,i]))
}
for_samp_prop <- aggregate(transposed_for_boxplot[,-1], by = list(transposed_for_boxplot[,1]), FUN = sum)
row.names(for_samp_prop) <- for_samp_prop[,1]
for_samp_prop <- for_samp_prop[,-1]
for_samp_prop <- for_samp_prop/rowSums(for_samp_prop)
mat_samp <- as.matrix(for_samp_prop)
pdf('plot_heat_map_proportion_class_groups.pdf')
heatmap.2(mat_samp,dendrogram = c("both"),key.title = "Proportion",key.xlab = "Proportion",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, main = "Proportion 16S class between groups",notecol="black", density.info="none", trace="none")
dev.off()
# log 10 proportion plot
pdf('plot_heat_map_proportion_class_groups_log_scale.pdf')
heatmap.2(log10(mat_samp +1),dendrogram = c("both"),key.title = "Z-score",
          key.xlab = "Proportion log scale",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,
          col= mypalette, main = "Proportion log scale 16S class between groups",notecol="black", 
          density.info="none", trace="none")
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
pdf('plot_heat_map_zscore_Class_groups.pdf')
heatmap.2(mat_samp_zscore,dendrogram = c("both"),key.title = "Proportion log scale +1",
          key.xlab = "Z-score",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, 
          main = "Z-score Class between groups",notecol="black", density.info="none", trace="none")
dev.off()
#
# Sample's correlation per class (heat map correlation)
cor_class_between_samples <- cor(transformed_df)
pdf('plot_heat_map_correlation_between_samples_class.pdf')
heatmap.2(cor_class_between_samples,dendrogram = c("both"),key.title = "Correlation",key.xlab = "Correlation",
          cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, main = "Correlation 16S class between samples",
          notecol="black", density.info="none", trace="none",ColSideColors = colors_in, RowSideColors = colors_in)
dev.off()
# Samples distances per class (pvclust of euclidean distances)
result_class_pvclust <- pvclust(end_file_class[,-1], method.dist = "euclidean", method.hclust="complete", nboot=1000)
pdf('rplot_pvclust_class.pdf')
plot(result_class_pvclust, main="class taxa plot", cex = .5, cex.pv=.3)
dev.off()
# class's association (Perform all comparisons class x class and save association tables' results)
# Rename columns in sense to make file understandable by lmGC function
colnames(transposed_taxa) <- sub("c__", "", colnames(transposed_taxa))
colnames(transposed_taxa) <- sub("\\[", "", colnames(transposed_taxa))
colnames(transposed_taxa) <- sub("\\]", "", colnames(transposed_taxa))
for (i in 2:ncol(transposed_taxa)){
  for (j in 2:ncol(transposed_taxa)){
    cname <- c(colnames(transposed_taxa[i]),colnames(transposed_taxa[j]))
    write.table(cname, file = "regression_class.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
    reg_value <- summary(lm(transposed_taxa[,i]~ transposed_taxa[,j] ,data=transposed_taxa))
    write.table(reg_value$coefficients,file = "regression_class.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    write.table(reg_value[9],file = "regression_class.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    write.table(reg_value[10],file = "regression_class.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
  }
}
#
#
#
# PCA plot and contribution analysis
# Create a data frame with the variavles (taxa names) as columns and the groups (samples) in lines
mat_pca<-merged_df[,c(5:ncol(merged_df))]
# PCA calculation
comp_Data<-prcomp(mat_pca)
# The PCA plot
pca_plot <- ggbiplot( comp_Data, groups = merged_df$Description,  obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = TRUE)
pca_plot <- pca_plot + scale_color_discrete(name = '')
pca_plot <- pca_plot + theme(legend.direction = 'horizontal', legend.position = 'top')
pdf('pca_plot_class.pdf')
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
taxa_contribution<-loadings %*% pca_comp_value # use the symbol "%" to separate the multiplication signal from the names of matices and vectors
# See each taxon contribution
taxa_contribution
# Check of the sum is 1 (or really close to it since the numbers approximation)
sum (taxa_contribution)
# Compare the value and the visual of the PCA plot. The largest contribution values should be represented as the longer arraows in the plot.
# Save the taxa contributions to a file
write.table(taxa_contribution, "taxa_contribution_class.txt",sep="\t", col.names = NA)
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
write.table(mult_tukey_res$group, "PERMANOVA_mult_tukey_res_class.txt", sep="\t", row.names = TRUE, col.names = NA)
# Now the permutation comparison to check groups homogeneity
per_class<-permutest(mod_disper,pairwise=TRUE, control=permControl(nperm=999))
per_class
write.table(per_class$pairwise$observed, "PERMANOVA_observed_permuted_pvalue_class.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
write.table(per_class$pairwise$permuted, "PERMANOVA_observed_permuted_pvalue_class.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
write.table(per_class$statistic, "PERMANOVA_observed_permuted_pvalue_class.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
# Calculate the general p-value
res_adonis<-adonis(mat_pca~merged_df$Description, permutations=999, method="euclidean")
write.table(res_adonis$aov.tab, "PERMANOVA_adonis_res_class.txt", sep="\t", row.names = TRUE, col.names = NA)
pdf("permanova_homogeneity_class.pdf")
plot(mod_disper, xlab = paste( "PCoA1 explains ", round((mod_disper$eig[1]/sum(mod_disper$eig)*100), digits = 5), "% of sample's dispersion", sep=""), ylab = paste( "PCoA2 explains ", round((mod_disper$eig[2]/sum(mod_disper$eig)*100), digits = 5), "% of sample's dispersion", sep=""), main = "class dispersion plot")
dev.off()
#
#
#
####Order###
order_no_Unclassified <- subset(input,!(input$Level_4=="Unclassified"))
order_no_Unclassified_edt <-cbind (order_no_Unclassified$Level_4, order_no_Unclassified[,!lapply(order_no_Unclassified, class) == "factor"])
# One row per taxon...
combine__rows_order_m <- aggregate(order_no_Unclassified_edt[,-1], by = list(order_no_Unclassified_edt[,1]), FUN = sum)
# Prepare to and transform data again... We are using the same variable names (overwriting) except for the final transformed file
# Retrieve only the values
numbers_to_norm <- combine__rows_order_m[,!lapply(combine__rows_order_m, class) == "factor"]
# Transform each column to proportion.
prop_otus<-t(t(numbers_to_norm)/rowSums(t(numbers_to_norm)))
# Check if each column's sum is 1.
colSums(prop_otus)
# Convert prop_otus to data frame format.
df_prop_otus<-as.data.frame(prop_otus)
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
end_file_order <- cbind (combine__rows_order_m[,lapply(combine__rows_order_m, class) == "factor"], transformed_df)
# Save the transformed output file in case you need it later.
write.table (end_file, 'order_log_ratio_transformed_data.txt', sep = '\t', row.names = FALSE)
#
#
#
# Create a barplot with top five taxa based on the first group ...
combine__rows_order_m2 <- t(combine__rows_order_m)
colnames(combine__rows_order_m2)<- combine__rows_order_m2[1,]
combine__rows_order_m2 <- data.frame(combine__rows_order_m2)
combine__rows_order_m2 <- combine__rows_order_m2[-1,]
combine__rows_order_m2[, ncol(combine__rows_order_m2)+1] <- NA
for (i in levels(map_groups$Description)){
  data_i <- as.vector(subset(map_groups, Description == i)[,1])
  for (name_col in data_i) {
    combine__rows_order_m2[which(row.names(combine__rows_order_m2) == name_col), ncol(combine__rows_order_m2)] <- i
  }
}
colnames(combine__rows_order_m2)[ncol(combine__rows_order_m2)] <- c("group_ID")
for (i in 1:(ncol(combine__rows_order_m2)-1)){
  combine__rows_order_m2[,i] <- as.numeric(as.character(combine__rows_order_m2[,i]))
}
combine__rows_order_m2 <- aggregate(combine__rows_order_m2[,-ncol(combine__rows_order_m2)], by = list(combine__rows_order_m2[,ncol(combine__rows_order_m2)]), FUN = mean)
t_combine__rows_order_m2 <- t(combine__rows_order_m2)
colnames(t_combine__rows_order_m2) <- t_combine__rows_order_m2[1,]
t_combine__rows_order_m2 <- t_combine__rows_order_m2[-1,]
t_combine__rows_order_m2 <- data.frame(t_combine__rows_order_m2)
for (i in 1:ncol(t_combine__rows_order_m2)){
  t_combine__rows_order_m2[,i] <- as.numeric(as.character(t_combine__rows_order_m2[,i]))
}
prop_order_data <- t(t(t_combine__rows_order_m2)/rowSums(t(t_combine__rows_order_m2)))
colSums(prop_order_data)
prop_order_data <- prop_order_data[order(prop_order_data[,1], decreasing = TRUE),]
pdf('barplot_top5_order.pdf')
par(xpd = NA)
barplot(prop_order_data[1:5,],main = "Top 5 most representative Order",
        xlab = "Groups", cex.names = .3,  ylab = "Proportion", yaxt="n",
        col = color_data[1:5], ylim=c(0,1.2))
ticks<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
axis(2,at=ticks,labels=ticks)
legend("topright", row.names(prop_order_data[1:5,]), cex=0.5, fill = color_data[1:5])
dev.off()
#
#
#
# Normality tests and comparisons
# Transpose the end_file_taxa Level table and create the transposed data frame
transposed_taxa <- t(end_file_order)
colnames(transposed_taxa) <- transposed_taxa[1,]
transposed_taxa <- transposed_taxa[-1,]
transposed_taxa <- as.data.frame(transposed_taxa)
# Make sure that values in data frame are numeric
for (i in 1:ncol(transposed_taxa)){
  transposed_taxa[,i] <- as.numeric(as.character(transposed_taxa[,i]))
}
# Change samples per groups names
merged_df <- merge(map_groups, transposed_taxa, by.x='X.SampleID', by.y= 'row.names')
transposed_taxa <- merged_df[,c(4:ncol(merged_df))]
# Perform normality test and the adequate comparisons
order_taxa <- c()
order_t <- c()
order_fdr <- c()
order_taxa_wil <- c()
order_wil_w <- c()
order_wil_fdr <- c()
for (i in 2:ncol(transposed_taxa)){
  z<-colnames(transposed_taxa[i])
  x <- shapiro.test(transposed_taxa[,i])
  if (length(levels(transposed_taxa$Description)) <= 2){
    warning("You don't have more than two groups")} else if (x$p.value >= 0.05 && length(levels(transposed_taxa$Description)) >= 3){
      tuk_res <- TukeyHSD(aov(transposed_taxa[,i] ~ transposed_taxa$Description))
      write.table(z, file = "comparisons_order_tukey.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
      tuk_res$`transposed_taxa$Description`[,4] <- p.adjust(tuk_res$`transposed_taxa$Description`[,4],method = 'fdr')
      write.table(tuk_res$`transposed_taxa$Description`, file = "comparisons_order_tukey.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    } else if (x$p.value < 0.05 && length(levels(transposed_taxa$Description)) >= 3){
      krusk_res <- kwAllPairsDunnTest(x=transposed_taxa[,i], g=transposed_taxa$Description, p.adjust.method="fdr")
      write.table(z, file = "comparisons_order_kruskal.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
      write.table(krusk_res$method,file = "comparisons_order_kruskal.txt", append = TRUE, row.names = TRUE, col.names = TRUE, sep='\t')
      write.table(krusk_res$p.value,file = "comparisons_order_kruskal.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    }
  if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (x$p.value >= 0.05 && length(levels(transposed_taxa$Description)) == 2){
      ttest_res <- t.test(transposed_taxa[,i] ~ transposed_taxa$Description)
      order_taxa <- c(order_taxa, z)
      order_t <- c(order_t, ttest_res[1])
      order_fdr <- c(order_fdr, ttest_res[3])
      order_ttest_sum <- cbind(order_taxa, order_t, order_fdr)
      colnames(order_ttest_sum)[1:3] <- c('order', 't_value', 'FDR')
    } else if (x$p.value < 0.05 && length(levels(transposed_taxa$Description)) == 2){
      wilcox_res <-wilcox.test(transposed_taxa[,i] ~ transposed_taxa$Description)
      order_taxa_wil <- c(order_taxa_wil, z)
      order_wil_w <- c(order_wil_w, wilcox_res[1])
      order_wil_fdr <- c(order_wil_fdr, wilcox_res[3])
      order_wil_sum <- cbind(order_taxa_wil, order_wil_w, order_wil_fdr)
      colnames(order_wil_sum)[1:3] <- c('order', 'w_value', 'FDR')
    }
}
#
if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (exists('order_ttest_sum') == TRUE){
        order_ttest_sum[,3] <- p.adjust(order_ttest_sum[,3], method = 'fdr')
        write.table(order_ttest_sum, file = "comparisons_order_t_test.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
    }
if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (exists('order_wil_sum') == TRUE){
        order_wil_sum[,3] <- p.adjust(order_wil_sum[,3], method = 'fdr')
        write.table(order_wil_sum, file = "comparisons_order_wilcox.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
    }
#
#
# Create the boxplot for each taxa per comparison group
for_boxplot<-cbind (combine__rows_order_m[,lapply(combine__rows_order_m, class) == "factor"], df_prop_otus)
transposed_for_boxplot <- t(for_boxplot)
colnames(transposed_for_boxplot) <- transposed_for_boxplot[1,]
transposed_for_boxplot <- transposed_for_boxplot[-1,]
transposed_for_boxplot <- as.data.frame(transposed_for_boxplot)
transposed_for_boxplot <- merge (map_groups, transposed_for_boxplot,by.x = "X.SampleID", by.y="row.names")
transposed_for_boxplot <- transposed_for_boxplot[,c(4,5:ncol(transposed_for_boxplot))]
for (i in 2:ncol(transposed_for_boxplot)){
  name = paste("boxplot_order_",i-1,".pdf",sep="")
  pdf(name, width=10, height=4)
  boxplot(as.numeric(as.character(transposed_for_boxplot[,i])) ~ Description, data = transposed_for_boxplot, main = colnames(transposed_for_boxplot[i]), cex.axis=.5, ylab = "Proportion", xlab = "groups", col="gray")
  dev.off()
}
#
#
#
# Plot the heat map of taxa proportion per groups
for (i in 2:ncol(transposed_for_boxplot)){
  transposed_for_boxplot[,i] <- as.numeric(as.character(transposed_for_boxplot[,i]))
}
for_samp_prop <- aggregate(transposed_for_boxplot[,-1], by = list(transposed_for_boxplot[,1]), FUN = sum)
row.names(for_samp_prop) <- for_samp_prop[,1]
for_samp_prop <- for_samp_prop[,-1]
for_samp_prop <- for_samp_prop/rowSums(for_samp_prop)
mat_samp <- as.matrix(for_samp_prop)
pdf('plot_heat_map_proportion_order_groups.pdf')
heatmap.2(mat_samp,dendrogram = c("both"),key.title = "Proportion",key.xlab = "Proportion",cexRow=0.3,Colv="Rowv",
          cexCol=0.3,keysize = 1,col= mypalette, main = "Proportion 16S order between groups",notecol="black", 
          density.info="none", trace="none")
dev.off()
# log 10 proportion plot
pdf('plot_heat_map_proportion_order_groups_log_scale.pdf')
heatmap.2(log10(mat_samp +1),dendrogram = c("both"),key.title = "Proportion log scale +1",
          key.xlab = "Proportion log scale",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, 
          main = "Proportion log scale 16S order between groups",notecol="black", density.info="none", trace="none")
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
pdf('plot_heat_map_zscore_Order_groups.pdf')
heatmap.2(mat_samp_zscore,dendrogram = c("both"),key.title = "Z-score",
          key.xlab = "Z-score",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, 
          main = "Z-score Order between groups",notecol="black", density.info="none", trace="none")
dev.off()
#
# Sample's correlation per order (heat map correlation)
cor_order_between_samples <- cor(transformed_df)
pdf('plot_heat_map_correlation_between_samples_order.pdf')
heatmap.2(cor_order_between_samples,dendrogram = c("both"),key.title = "Correlation",key.xlab = "Correlation",
          cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, main = "Correlation 16S Order between samples",
          notecol="black", density.info="none", trace="none",ColSideColors = colors_in, 
          RowSideColors = colors_in)
dev.off()
# Samples distances per order (pvclust of euclidean distances)
result_order_pvclust <- pvclust(end_file_order[,-1], method.dist = "euclidean", method.hclust="complete", nboot=1000)
pdf('rplot_pvclust_order.pdf')
plot(result_order_pvclust, main="order taxa plot", cex = .5, cex.pv=.3)
dev.off()
# order's association (Perform all comparisons order x order and save association tables' results)
# Rename columns in sense to make file understandable by lmGC function
colnames(transposed_taxa) <- sub("o__", "", colnames(transposed_taxa))
colnames(transposed_taxa) <- sub("\\[", "", colnames(transposed_taxa))
colnames(transposed_taxa) <- sub("\\]", "", colnames(transposed_taxa))
for (i in 2:ncol(transposed_taxa)){
  for (j in 2:ncol(transposed_taxa)){
    cname <- c(colnames(transposed_taxa[i]),colnames(transposed_taxa[j]))
    write.table(cname, file = "regression_order.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
    reg_value <- summary(lm(transposed_taxa[,i]~ transposed_taxa[,j] ,data=transposed_taxa))
    write.table(reg_value$coefficients,file = "regression_order.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    write.table(reg_value[9],file = "regression_order.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    write.table(reg_value[10],file = "regression_order.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
  }
}
#
#
#
# PCA plot and contribution analysis
# Create a data frame with the variavles (taxa names) as columns and the groups (samples) in lines
mat_pca<-merged_df[,c(5:ncol(merged_df))]
# PCA calculation
comp_Data<-prcomp(mat_pca)
# The PCA plot
pca_plot <- ggbiplot( comp_Data, groups = merged_df$Description,  obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = TRUE)
pca_plot <- pca_plot + scale_color_discrete(name = '')
pca_plot <- pca_plot + theme(legend.direction = 'horizontal', legend.position = 'top')
pdf('pca_plot_order.pdf')
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
taxa_contribution<-loadings %*% pca_comp_value # use the symbol "%" to separate the multiplication signal from the names of matices and vectors
# See each taxon contribution
taxa_contribution
# Check of the sum is 1 (or really close to it since the numbers approximation)
sum (taxa_contribution)
# Compare the value and the visual of the PCA plot. The largest contribution values should be represented as the longer arraows in the plot.
# Save the taxa contributions to a file
write.table(taxa_contribution, "taxa_contribution_order.txt",sep="\t", col.names = NA)
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
write.table(mult_tukey_res$group, "PERMANOVA_mult_tukey_res_order.txt", sep="\t", row.names = TRUE, col.names = NA)
# Now the permutation comparison to check groups homogeneity
per_order<-permutest(mod_disper,pairwise=TRUE, control=permControl(nperm=999))
per_order
write.table(per_order$pairwise$observed, "PERMANOVA_observed_permuted_pvalue_order.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
write.table(per_order$pairwise$permuted, "PERMANOVA_observed_permuted_pvalue_order.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
write.table(per_order$statistic, "PERMANOVA_observed_permuted_pvalue_order.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
# Calculate the general p-value
res_adonis<-adonis(mat_pca~merged_df$Description, permutations=999, method="euclidean")
write.table(res_adonis$aov.tab, "PERMANOVA_adonis_res_order.txt", sep="\t", row.names = TRUE, col.names = NA)
pdf("permanova_homogeneity_order.pdf")
plot(mod_disper, xlab = paste( "PCoA1 explains ", round((mod_disper$eig[1]/sum(mod_disper$eig)*100), digits = 5), "% of sample's dispersion", sep=""), ylab = paste( "PCoA2 explains ", round((mod_disper$eig[2]/sum(mod_disper$eig)*100), digits = 5), "% of sample's dispersion", sep=""), main = "order dispersion plot")
dev.off()
#
#
#
####Family###
family_no_Unclassified <- subset(input,!(input$Level_5=="Unclassified"))
family_no_Unclassified_edt <-cbind (family_no_Unclassified$Level_5, family_no_Unclassified[,!lapply(family_no_Unclassified, class) == "factor"])
# One row per taxon...
combine__rows_family_m <- aggregate(family_no_Unclassified_edt[,-1], by = list(family_no_Unclassified_edt[,1]), FUN = sum)
# Prepare to and transform data again... We are using the same variable names (overwriting) except for the final transformed file
# Retrieve only the values
numbers_to_norm <- combine__rows_family_m[,!lapply(combine__rows_family_m, class) == "factor"]
# Transform each column to proportion.
prop_otus<-t(t(numbers_to_norm)/rowSums(t(numbers_to_norm)))
# Check if each column's sum is 1.
colSums(prop_otus)
# Convert prop_otus to data frame format.
df_prop_otus<-as.data.frame(prop_otus)
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
end_file_family <- cbind (combine__rows_family_m[,lapply(combine__rows_family_m, class) == "factor"], transformed_df)
# Save the transformed output file in case you need it later.
write.table (end_file, 'family_log_ratio_transformed_data.txt', sep = '\t', row.names = FALSE)
#
#
#
# Create a barplot with top five taxa based on the first group ...
combine__rows_family_m2 <- t(combine__rows_family_m)
colnames(combine__rows_family_m2)<- combine__rows_family_m2[1,]
combine__rows_family_m2 <- data.frame(combine__rows_family_m2)
combine__rows_family_m2 <- combine__rows_family_m2[-1,]
combine__rows_family_m2[, ncol(combine__rows_family_m2)+1] <- NA
for (i in levels(map_groups$Description)){
  data_i <- as.vector(subset(map_groups, Description == i)[,1])
  for (name_col in data_i) {
    combine__rows_family_m2[which(row.names(combine__rows_family_m2) == name_col), ncol(combine__rows_family_m2)] <- i
  }
}
colnames(combine__rows_family_m2)[ncol(combine__rows_family_m2)] <- c("group_ID")
for (i in 1:(ncol(combine__rows_family_m2)-1)){
  combine__rows_family_m2[,i] <- as.numeric(as.character(combine__rows_family_m2[,i]))
}
combine__rows_family_m2 <- aggregate(combine__rows_family_m2[,-ncol(combine__rows_family_m2)], by = list(combine__rows_family_m2[,ncol(combine__rows_family_m2)]), FUN = mean)
t_combine__rows_family_m2 <- t(combine__rows_family_m2)
colnames(t_combine__rows_family_m2) <- t_combine__rows_family_m2[1,]
t_combine__rows_family_m2 <- t_combine__rows_family_m2[-1,]
t_combine__rows_family_m2 <- data.frame(t_combine__rows_family_m2)
for (i in 1:ncol(t_combine__rows_family_m2)){
  t_combine__rows_family_m2[,i] <- as.numeric(as.character(t_combine__rows_family_m2[,i]))
}
prop_family_data <- t(t(t_combine__rows_family_m2)/rowSums(t(t_combine__rows_family_m2)))
colSums(prop_family_data)
prop_family_data <- prop_family_data[order(prop_family_data[,1], decreasing = TRUE),]
pdf('barplot_top5_family.pdf')
par(xpd = NA)
barplot(prop_family_data[1:5,],main = "Top 5 most representative Family",
        xlab = "Groups", cex.names = .3,  ylab = "Proportion", yaxt="n",
        col = color_data[1:5], ylim=c(0,1.2))
ticks<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
axis(2,at=ticks,labels=ticks)
legend("topright", row.names(prop_family_data[1:5,]), cex=0.5, fill = color_data[1:5])
dev.off()
#
#
#
# Normality tests and comparisons
# Transpose the end_file_taxa Level table and create the transposed data frame
transposed_taxa <- t(end_file_family)
colnames(transposed_taxa) <- transposed_taxa[1,]
transposed_taxa <- transposed_taxa[-1,]
transposed_taxa <- as.data.frame(transposed_taxa)
# Make sure that values in data frame are numeric
for (i in 1:ncol(transposed_taxa)){
  transposed_taxa[,i] <- as.numeric(as.character(transposed_taxa[,i]))
}
# Change samples per groups names
merged_df <- merge(map_groups, transposed_taxa, by.x='X.SampleID', by.y= 'row.names')
transposed_taxa <- merged_df[,c(4:ncol(merged_df))]
# Perform normality test and the adequate comparisons
family_taxa <- c()
family_t <- c()
family_fdr <- c()
family_taxa_wil <- c()
family_wil_w <- c()
family_wil_fdr <- c()
for (i in 2:ncol(transposed_taxa)){
  z<-colnames(transposed_taxa[i])
  x <- shapiro.test(transposed_taxa[,i])
  if (length(levels(transposed_taxa$Description)) <= 2){
    warning("You don't have more than two groups")} else if (x$p.value >= 0.05 && length(levels(transposed_taxa$Description)) >= 3){
      tuk_res <- TukeyHSD(aov(transposed_taxa[,i] ~ transposed_taxa$Description))
      write.table(z, file = "comparisons_family_tukey.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
      tuk_res$`transposed_taxa$Description`[,4] <- p.adjust(tuk_res$`transposed_taxa$Description`[,4],method = 'fdr')
      write.table(tuk_res$`transposed_taxa$Description`, file = "comparisons_family_tukey.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    } else if (x$p.value < 0.05 && length(levels(transposed_taxa$Description)) >= 3){
      krusk_res <- kwAllPairsDunnTest(x=transposed_taxa[,i], g=transposed_taxa$Description, p.adjust.method="fdr")
      write.table(z, file = "comparisons_family_kruskal.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
      write.table(krusk_res$method,file = "comparisons_family_kruskal.txt", append = TRUE, row.names = TRUE, col.names = TRUE, sep='\t')
      write.table(krusk_res$p.value,file = "comparisons_family_kruskal.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    }
  if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (x$p.value >= 0.05 && length(levels(transposed_taxa$Description)) == 2){
      ttest_res <- t.test(transposed_taxa[,i] ~ transposed_taxa$Description)
      family_taxa <- c(family_taxa, z)
      family_t <- c(family_t, ttest_res[1])
      family_fdr <- c(family_fdr, ttest_res[3])
      family_ttest_sum <- cbind(family_taxa, family_t, family_fdr)
      colnames(family_ttest_sum)[1:3] <- c('family', 't_value', 'FDR')
    } else if (x$p.value < 0.05 && length(levels(transposed_taxa$Description)) == 2){
      wilcox_res <-wilcox.test(transposed_taxa[,i] ~ transposed_taxa$Description)
      family_taxa_wil <- c(family_taxa_wil, z)
      family_wil_w <- c(family_wil_w, wilcox_res[1])
      family_wil_fdr <- c(family_wil_fdr, wilcox_res[3])
      family_wil_sum <- cbind(family_taxa_wil, family_wil_w, family_wil_fdr)
      colnames(family_wil_sum)[1:3] <- c('family', 'w_value', 'FDR')
    }
}
#
if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (exists('family_ttest_sum') == TRUE){
        family_ttest_sum[,3] <- p.adjust(family_ttest_sum[,3], method = 'fdr')
        write.table(family_ttest_sum, file = "comparisons_family_t_test.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
    }
if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (exists('family_wil_sum') == TRUE){
        family_wil_sum[,3] <- p.adjust(family_wil_sum[,3], method = 'fdr')
        write.table(family_wil_sum, file = "comparisons_family_wilcox.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
    }
#
#
# Create the boxplot for each taxa per comparison group
for_boxplot<-cbind (combine__rows_family_m[,lapply(combine__rows_family_m, class) == "factor"], df_prop_otus)
transposed_for_boxplot <- t(for_boxplot)
colnames(transposed_for_boxplot) <- transposed_for_boxplot[1,]
transposed_for_boxplot <- transposed_for_boxplot[-1,]
transposed_for_boxplot <- as.data.frame(transposed_for_boxplot)
transposed_for_boxplot <- merge (map_groups, transposed_for_boxplot,by.x = "X.SampleID", by.y="row.names")
transposed_for_boxplot <- transposed_for_boxplot[,c(4,5:ncol(transposed_for_boxplot))]
for (i in 2:ncol(transposed_for_boxplot)){
  name = paste("boxplot_family_",i-1,".pdf",sep="")
  pdf(name, width=10, height=4)
  boxplot(as.numeric(as.character(transposed_for_boxplot[,i])) ~ Description, data = transposed_for_boxplot, main = colnames(transposed_for_boxplot[i]), cex.axis=.5, ylab = "Proportion", xlab = "groups", col="gray")
  dev.off()
}
#
#
#
# Plot the heat map of taxa proportion per groups
for (i in 2:ncol(transposed_for_boxplot)){
  transposed_for_boxplot[,i] <- as.numeric(as.character(transposed_for_boxplot[,i]))
}
for_samp_prop <- aggregate(transposed_for_boxplot[,-1], by = list(transposed_for_boxplot[,1]), FUN = sum)
row.names(for_samp_prop) <- for_samp_prop[,1]
for_samp_prop <- for_samp_prop[,-1]
for_samp_prop <- for_samp_prop/rowSums(for_samp_prop)
mat_samp <- as.matrix(for_samp_prop)
pdf('plot_heat_map_proportion_family_groups.pdf')
heatmap.2(mat_samp,dendrogram = c("both"),key.title = "Proportion",key.xlab = "Proportion",cexRow=0.3,Colv="Rowv",
          cexCol=0.3,keysize = 1,col= mypalette, main = "Proportion 16S family between groups",notecol="black", 
          density.info="none", trace="none")
dev.off()
# log 10 proportion plot
pdf('plot_heat_map_proportion_family_groups_log_scale.pdf')
heatmap.2(log10(mat_samp +1),dendrogram = c("both"),key.title = "Proportion log scale +1",
          key.xlab = "Proportion log scale",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, 
          main = "Proportion log scale 16S family between groups",notecol="black", density.info="none", trace="none")
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
pdf('plot_heat_map_zscore_Family_groups.pdf')
heatmap.2(mat_samp_zscore,dendrogram = c("both"),key.title = "Z-score",
          key.xlab = "Z-score",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, 
          main = "Z-score Family between groups",notecol="black", density.info="none", trace="none")
dev.off()
#
# Sample's correlation per family (heat map correlation)
cor_family_between_samples <- cor(transformed_df)
pdf('plot_heat_map_correlation_between_samples_family.pdf')
heatmap.2(cor_family_between_samples,dendrogram = c("both"),key.title = "Correlation",key.xlab = "Correlation",
          cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, main = "Correlation 16S family between samples",
          notecol="black", density.info="none", trace="none",ColSideColors = colors_in, 
          RowSideColors = colors_in)
dev.off()
# Samples distances per family (pvclust of euclidean distances)
result_family_pvclust <- pvclust(end_file_family[,-1], method.dist = "euclidean", method.hclust="complete", nboot=1000)
pdf('rplot_pvclust_family.pdf')
plot(result_family_pvclust, main="family taxa plot", cex = .5, cex.pv=.3)
dev.off()
# family's association (Perform all comparisons family x family and save association tables' results)
# Rename columns in sense to make file understandable by lmGC function
colnames(transposed_taxa) <- sub("f__", "", colnames(transposed_taxa))
colnames(transposed_taxa) <- sub("\\[", "", colnames(transposed_taxa))
colnames(transposed_taxa) <- sub("\\]", "", colnames(transposed_taxa))
for (i in 2:ncol(transposed_taxa)){
  for (j in 2:ncol(transposed_taxa)){
    cname <- c(colnames(transposed_taxa[i]),colnames(transposed_taxa[j]))
    write.table(cname, file = "regression_family.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
    reg_value <- summary(lm(transposed_taxa[,i]~ transposed_taxa[,j] ,data=transposed_taxa))
    write.table(reg_value$coefficients,file = "regression_family.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    write.table(reg_value[9],file = "regression_family.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    write.table(reg_value[10],file = "regression_family.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
  }
}
#
#
#
# PCA plot and contribution analysis
# Create a data frame with the variavles (taxa names) as columns and the groups (samples) in lines
mat_pca<-merged_df[,c(5:ncol(merged_df))]
# PCA calculation
comp_Data<-prcomp(mat_pca)
# The PCA plot
pca_plot <- ggbiplot( comp_Data, groups = merged_df$Description,  obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = TRUE)
pca_plot <- pca_plot + scale_color_discrete(name = '')
pca_plot <- pca_plot + theme(legend.direction = 'horizontal', legend.position = 'top')
pdf('pca_plot_family.pdf')
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
taxa_contribution<-loadings %*% pca_comp_value # use the symbol "%" to separate the multiplication signal from the names of matices and vectors
# See each taxon contribution
taxa_contribution
# Check of the sum is 1 (or really close to it since the numbers approximation)
sum (taxa_contribution)
# Compare the value and the visual of the PCA plot. The largest contribution values should be represented as the longer arraows in the plot.
# Save the taxa contributions to a file
write.table(taxa_contribution, "taxa_contribution_family.txt",sep="\t", col.names = NA)
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
write.table(mult_tukey_res$group, "PERMANOVA_mult_tukey_res_family.txt", sep="\t", row.names = TRUE, col.names = NA)
# Now the permutation comparison to check groups homogeneity
per_family<-permutest(mod_disper,pairwise=TRUE, control=permControl(nperm=999))
per_family
write.table(per_family$pairwise$observed, "PERMANOVA_observed_permuted_pvalue_family.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
write.table(per_family$pairwise$permuted, "PERMANOVA_observed_permuted_pvalue_family.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
write.table(per_family$statistic, "PERMANOVA_observed_permuted_pvalue_family.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
# Calculate the general p-value
res_adonis<-adonis(mat_pca~merged_df$Description, permutations=999, method="euclidean")
write.table(res_adonis$aov.tab, "PERMANOVA_adonis_res_family.txt", sep="\t", row.names = TRUE, col.names = NA)
pdf("permanova_homogeneity_family.pdf")
plot(mod_disper, xlab = paste( "PCoA1 explains ", round((mod_disper$eig[1]/sum(mod_disper$eig)*100), digits = 5), "% of sample's dispersion", sep=""), ylab = paste( "PCoA2 explains ", round((mod_disper$eig[2]/sum(mod_disper$eig)*100), digits = 5), "% of sample's dispersion", sep=""), main = "family dispersion plot")
dev.off()
#
#
#
####Genus###
genus_no_Unclassified <- subset(input,!(input$Level_6=="Unclassified"))
genus_no_Unclassified_edt <-cbind (genus_no_Unclassified$Level_6, genus_no_Unclassified[,!lapply(genus_no_Unclassified, class) == "factor"])
# One row per taxon...
combine__rows_genus_m <- aggregate(genus_no_Unclassified_edt[,-1], by = list(genus_no_Unclassified_edt[,1]), FUN = sum)
# Prepare to and transform data again... We are using the same variable names (overwriting) except for the final transformed file
# Retrieve only the values
numbers_to_norm <- combine__rows_genus_m[,!lapply(combine__rows_genus_m, class) == "factor"]
# Transform each column to proportion.
prop_otus<-t(t(numbers_to_norm)/rowSums(t(numbers_to_norm)))
# Check if each column's sum is 1.
colSums(prop_otus)
# Convert prop_otus to data frame format.
df_prop_otus<-as.data.frame(prop_otus)
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
end_file_genus <- cbind (combine__rows_genus_m[,lapply(combine__rows_genus_m, class) == "factor"], transformed_df)
# Save the transformed output file in case you need it later.
write.table (end_file, 'genus_log_ratio_transformed_data.txt', sep = '\t', row.names = FALSE)
#
#
#
# Create a barplot with top five taxa based on the first group ...
combine__rows_genus_m2 <- t(combine__rows_genus_m)
colnames(combine__rows_genus_m2)<- combine__rows_genus_m2[1,]
combine__rows_genus_m2 <- data.frame(combine__rows_genus_m2)
combine__rows_genus_m2 <- combine__rows_genus_m2[-1,]
combine__rows_genus_m2[, ncol(combine__rows_genus_m2)+1] <- NA
for (i in levels(map_groups$Description)){
  data_i <- as.vector(subset(map_groups, Description == i)[,1])
  for (name_col in data_i) {
    combine__rows_genus_m2[which(row.names(combine__rows_genus_m2) == name_col), ncol(combine__rows_genus_m2)] <- i
  }
}
colnames(combine__rows_genus_m2)[ncol(combine__rows_genus_m2)] <- c("group_ID")
for (i in 1:(ncol(combine__rows_genus_m2)-1)){
  combine__rows_genus_m2[,i] <- as.numeric(as.character(combine__rows_genus_m2[,i]))
}
combine__rows_genus_m2 <- aggregate(combine__rows_genus_m2[,-ncol(combine__rows_genus_m2)], by = list(combine__rows_genus_m2[,ncol(combine__rows_genus_m2)]), FUN = mean)
t_combine__rows_genus_m2 <- t(combine__rows_genus_m2)
colnames(t_combine__rows_genus_m2) <- t_combine__rows_genus_m2[1,]
t_combine__rows_genus_m2 <- t_combine__rows_genus_m2[-1,]
t_combine__rows_genus_m2 <- data.frame(t_combine__rows_genus_m2)
for (i in 1:ncol(t_combine__rows_genus_m2)){
  t_combine__rows_genus_m2[,i] <- as.numeric(as.character(t_combine__rows_genus_m2[,i]))
}
prop_genus_data <- t(t(t_combine__rows_genus_m2)/rowSums(t(t_combine__rows_genus_m2)))
colSums(prop_genus_data)
prop_genus_data <- prop_genus_data[order(prop_genus_data[,1], decreasing = TRUE),]
pdf('barplot_top5_genus.pdf')
par(xpd = NA)
barplot(prop_genus_data[1:5,],main = "Top 5 most representative Genus",
        xlab = "Groups", cex.names = .3,  ylab = "Proportion", yaxt="n",
        col = color_data[1:5], ylim=c(0,1.2))
ticks<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
axis(2,at=ticks,labels=ticks)
legend("topright", row.names(prop_genus_data[1:5,]), cex=0.5, fill = color_data[1:5])
dev.off()
#
#
#
# Normality tests and comparisons
# Transpose the end_file_taxa Level table and create the transposed data frame
transposed_taxa <- t(end_file_genus)
colnames(transposed_taxa) <- transposed_taxa[1,]
transposed_taxa <- transposed_taxa[-1,]
transposed_taxa <- as.data.frame(transposed_taxa)
# Make sure that values in data frame are numeric
for (i in 1:ncol(transposed_taxa)){
  transposed_taxa[,i] <- as.numeric(as.character(transposed_taxa[,i]))
}
# Change samples per groups names
merged_df <- merge(map_groups, transposed_taxa, by.x='X.SampleID', by.y= 'row.names')
transposed_taxa <- merged_df[,c(4:ncol(merged_df))]
# Perform normality test and the adequate comparisons
genus_taxa <- c()
genus_t <- c()
genus_fdr <- c()
genus_taxa_wil <- c()
genus_wil_w <- c()
genus_wil_fdr <- c()
for (i in 2:ncol(transposed_taxa)){
  z<-colnames(transposed_taxa[i])
  x <- shapiro.test(transposed_taxa[,i])
  if (length(levels(transposed_taxa$Description)) <= 2){
    warning("You don't have more than two groups")} else if (x$p.value >= 0.05 && length(levels(transposed_taxa$Description)) >= 3){
      tuk_res <- TukeyHSD(aov(transposed_taxa[,i] ~ transposed_taxa$Description))
      write.table(z, file = "comparisons_genus_tukey.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
      tuk_res$`transposed_taxa$Description`[,4] <- p.adjust(tuk_res$`transposed_taxa$Description`[,4],method = 'fdr')
      write.table(tuk_res$`transposed_taxa$Description`, file = "comparisons_genus_tukey.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    } else if (x$p.value < 0.05 && length(levels(transposed_taxa$Description)) >= 3){
      krusk_res <- kwAllPairsDunnTest(x=transposed_taxa[,i], g=transposed_taxa$Description, p.adjust.method="fdr")
      write.table(z, file = "comparisons_genus_kruskal.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
      write.table(krusk_res$method,file = "comparisons_genus_kruskal.txt", append = TRUE, row.names = TRUE, col.names = TRUE, sep='\t')
      write.table(krusk_res$p.value,file = "comparisons_genus_kruskal.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    }
  if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (x$p.value >= 0.05 && length(levels(transposed_taxa$Description)) == 2){
      ttest_res <- t.test(transposed_taxa[,i] ~ transposed_taxa$Description)
      genus_taxa <- c(genus_taxa, z)
      genus_t <- c(genus_t, ttest_res[1])
      genus_fdr <- c(genus_fdr, ttest_res[3])
      genus_ttest_sum <- cbind(genus_taxa, genus_t, genus_fdr)
      colnames(genus_ttest_sum)[1:3] <- c('genus', 't_value', 'FDR')
    } else if (x$p.value < 0.05 && length(levels(transposed_taxa$Description)) == 2){
      wilcox_res <-wilcox.test(transposed_taxa[,i] ~ transposed_taxa$Description)
      genus_taxa_wil <- c(genus_taxa_wil, z)
      genus_wil_w <- c(genus_wil_w, wilcox_res[1])
      genus_wil_fdr <- c(genus_wil_fdr, wilcox_res[3])
      genus_wil_sum <- cbind(genus_taxa_wil, genus_wil_w, genus_wil_fdr)
      colnames(genus_wil_sum)[1:3] <- c('genus', 'w_value', 'FDR')
    }
}
#
if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (exists('genus_ttest_sum') == TRUE){
        genus_ttest_sum[,3] <- p.adjust(genus_ttest_sum[,3], method = 'fdr')
        write.table(genus_ttest_sum, file = "comparisons_genus_t_test.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
    }
if (length(levels(transposed_taxa$Description)) != 2){
    warning("You don't have two groups")} else if (exists('genus_wil_sum') == TRUE){
        genus_wil_sum[,3] <- p.adjust(genus_wil_sum[,3], method = 'fdr')
        write.table(genus_wil_sum, file = "comparisons_genus_wilcox.txt", append = FALSE, row.names = FALSE, col.names = TRUE, sep='\t')
    }
#
#
# Create the boxplot for each taxa per comparison group
for_boxplot<-cbind (combine__rows_genus_m[,lapply(combine__rows_genus_m, class) == "factor"], df_prop_otus)
transposed_for_boxplot <- t(for_boxplot)
colnames(transposed_for_boxplot) <- transposed_for_boxplot[1,]
transposed_for_boxplot <- transposed_for_boxplot[-1,]
transposed_for_boxplot <- as.data.frame(transposed_for_boxplot)
transposed_for_boxplot <- merge (map_groups, transposed_for_boxplot,by.x = "X.SampleID", by.y="row.names")
transposed_for_boxplot <- transposed_for_boxplot[,c(4,5:ncol(transposed_for_boxplot))]
for (i in 2:ncol(transposed_for_boxplot)){
  name = paste("boxplot_genus_",i-1,".pdf",sep="")
  pdf(name, width=10, height=4)
  boxplot(as.numeric(as.character(transposed_for_boxplot[,i])) ~ Description, data = transposed_for_boxplot, main = colnames(transposed_for_boxplot[i]), cex.axis=.5, ylab = "Proportion", xlab = "groups", col="gray")
  dev.off()
}
#
#
#
# Plot the heat map of taxa proportion per groups
for (i in 2:ncol(transposed_for_boxplot)){
  transposed_for_boxplot[,i] <- as.numeric(as.character(transposed_for_boxplot[,i]))
}
for_samp_prop <- aggregate(transposed_for_boxplot[,-1], by = list(transposed_for_boxplot[,1]), FUN = sum)
row.names(for_samp_prop) <- for_samp_prop[,1]
for_samp_prop <- for_samp_prop[,-1]
for_samp_prop <- for_samp_prop/rowSums(for_samp_prop)
mat_samp <- as.matrix(for_samp_prop)
pdf('plot_heat_map_proportion_genus_groups.pdf')
heatmap.2(mat_samp,dendrogram = c("both"),key.title = "Proportion",key.xlab = "Proportion",cexRow=0.3,Colv="Rowv",
          cexCol=0.3,keysize = 1,col= mypalette, main = "Proportion 16S genus between groups",notecol="black", 
          density.info="none", trace="none")
dev.off()
# log 10 proportion plot
pdf('plot_heat_map_proportion_genus_groups_log_scale.pdf')
heatmap.2(log10(mat_samp +1),dendrogram = c("both"),key.title = "Proportion log scale +1",
          key.xlab = "Proportion log scale",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, 
          main = "Proportion log scale 16S genus between groups",notecol="black", density.info="none", trace="none")
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
pdf('plot_heat_map_zscore_Genus_groups.pdf')
heatmap.2(mat_samp_zscore,dendrogram = c("both"),key.title = "Z-score",
          key.xlab = "Z-score",cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, 
          main = "Z-score Genus between groups",notecol="black", density.info="none", trace="none")
dev.off()
#
# Sample's correlation per genus (heat map correlation)
cor_genus_between_samples <- cor(transformed_df)
pdf('plot_heat_map_correlation_between_samples_genus.pdf')
heatmap.2(cor_genus_between_samples,dendrogram = c("both"),key.title = "Correlation",key.xlab = "Correlation",
          cexRow=0.3,Colv="Rowv",cexCol=0.3,keysize = 1,col= mypalette, main = "Correlation 16S genus between samples",
          notecol="black", density.info="none", trace="none", ColSideColors = colors_in, 
          RowSideColors = colors_in)
dev.off()
# Samples distances per genus (pvclust of euclidean distances)
result_genus_pvclust <- pvclust(end_file_genus[,-1], method.dist = "euclidean", method.hclust="complete", nboot=1000)
pdf('rplot_pvclust_genus.pdf')
plot(result_genus_pvclust, main="genus taxa plot", cex = .5, cex.pv=.3)
dev.off()
# genus's association (Perform all comparisons genus x genus and save association tables' results)
# Rename columns in sense to make file understandable by lmGC function
colnames(transposed_taxa) <- sub("g__", "", colnames(transposed_taxa))
colnames(transposed_taxa) <- sub("\\[", "", colnames(transposed_taxa))
colnames(transposed_taxa) <- sub("\\]", "", colnames(transposed_taxa))
for (i in 2:ncol(transposed_taxa)){
  for (j in 2:ncol(transposed_taxa)){
    cname <- c(colnames(transposed_taxa[i]),colnames(transposed_taxa[j]))
    write.table(cname, file = "regression_genus.txt", append = TRUE, row.names = FALSE, col.names = FALSE, sep='\t')
    reg_value <- summary(lm(transposed_taxa[,i]~ transposed_taxa[,j] ,data=transposed_taxa))
    write.table(reg_value$coefficients,file = "regression_genus.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    write.table(reg_value[9],file = "regression_genus.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
    write.table(reg_value[10],file = "regression_genus.txt", append = TRUE, row.names = TRUE, col.names = NA, sep='\t')
  }
}
#
#
#
# PCA plot and contribution analysis
# Create a data frame with the variavles (taxa names) as columns and the groups (samples) in lines
mat_pca<-merged_df[,c(5:ncol(merged_df))]
# PCA calculation
comp_Data<-prcomp(mat_pca)
# The PCA plot
pca_plot <- ggbiplot( comp_Data, groups = merged_df$Description,  obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = TRUE)
pca_plot <- pca_plot + scale_color_discrete(name = '')
pca_plot <- pca_plot + theme(legend.direction = 'horizontal', legend.position = 'top')
pdf('pca_plot_genus.pdf')
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
taxa_contribution<-loadings %*% pca_comp_value # use the symbol "%" to separate the multiplication signal from the names of matices and vectors
# See each taxon contribution
taxa_contribution
# Check of the sum is 1 (or really close to it since the numbers approximation)
sum (taxa_contribution)
# Compare the value and the visual of the PCA plot. The largest contribution values should be represented as the longer arraows in the plot.
# Save the taxa contributions to a file
write.table(taxa_contribution, "taxa_contribution_genus.txt",sep="\t", col.names = NA)
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
write.table(mult_tukey_res$group, "PERMANOVA_mult_tukey_res_genus.txt", sep="\t", row.names = TRUE, col.names = NA)
# Now the permutation comparison to check groups homogeneity
per_genus<-permutest(mod_disper,pairwise=TRUE, control=permControl(nperm=999))
per_genus
write.table(per_genus$pairwise$observed, "PERMANOVA_observed_permuted_pvalue_genus.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
write.table(per_genus$pairwise$permuted, "PERMANOVA_observed_permuted_pvalue_genus.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
write.table(per_genus$statistic, "PERMANOVA_observed_permuted_pvalue_genus.txt", append = TRUE, sep="\t", row.names = TRUE, col.names = NA)
# Calculate the general p-value
res_adonis<-adonis(mat_pca~merged_df$Description, permutations=999, method="euclidean")
write.table(res_adonis$aov.tab, "PERMANOVA_adonis_res_genus.txt", sep="\t", row.names = TRUE, col.names = NA)
pdf("permanova_homogeneity_genus.pdf")
plot(mod_disper, xlab = paste( "PCoA1 explains ", round((mod_disper$eig[1]/sum(mod_disper$eig)*100), digits = 5), "% of sample's dispersion", sep=""), ylab = paste( "PCoA2 explains ", round((mod_disper$eig[2]/sum(mod_disper$eig)*100), digits = 5), "% of sample's dispersion", sep=""), main = "genus dispersion plot")
dev.off()
save.image("metagenomic_results.RData")
#
end_time <- Sys.time()
end_time - start_time
