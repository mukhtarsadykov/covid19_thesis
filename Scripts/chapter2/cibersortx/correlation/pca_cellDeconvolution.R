### Put into order ICU and NON-ICU samples for cell deconvolution ###
# I used bRNAseq all 80 samples for this purpose
a1 <- read.csv("~/Desktop/rnaseq/last/raw/cts_total.csv", header = T, row.names = 1)

ensg <- sub("\\..*", "", rownames(a1))  # remove version number in case you have it
sym <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys=ensg,
                             column="SYMBOL", keytype="GENEID") # Unable to map 3533 of 60663 requested IDs.
gene <- data.frame(ENSGID=ensg, SYMBOL=sym, stringsAsFactors=F)
rownames(gene) <- rownames(a1)
a2<-data.frame(gene, a1)
a2<-a2[,-1]
a2<-a2[complete.cases(a2), ]
dim(a2) #57131    81
a2<-a2[!duplicated(a2$SYMBOL), ] #55418    33
dim(a2) #55418    81
rownames(a2) <- a2$SYMBOL

meta <- read.csv("~/Desktop/rnaseq/last/raw/meta_1_2_3_4_5.csv", header = T)

a3 <- write_delim(a2[, c("SYMBOL", meta$SampleID)], delim = "\t",
                 "~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/cts_ordered_80s.txt")
# Date: 2023-11-04 05:43:47
# Job type: Impute Cell Fractions
# Signature matrix file: CIBERSORTx_Job6_scRNA_count_cybersrotx_all_typeCells_1_inferred_phenoclasses.CIBERSORTx_Job6_scRNA_count_cybersrotx_all_typeCells_1_inferred_refsample.bm.K999.txt
# Mixture file: cts_ordered_80s.txt
# Batch correction: enabled
# Batch correction mode: B-mode
# Disable quantile normalization: true
# Run mode (relative or absolute): relative
# Permutations: 500

## Change SampleID to ID (ICU_01...)

a4 <- read.csv("~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/CIBERSORTx_Job18_Results.csv",
               header = T)
a5<-merge(a4, meta, by = "SampleID")
write.csv(a5, "~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/CIBERSORx_merged_meta.csv")

# PCA for cell proportions #

a6 <- read.csv("~/Desktop/phd/Thesis/ch2/cibersortx/results/CIBERSORTx_Job7_Results-2.csv", header = T)

a6  <- a6[, colSums(a6 != 0) > 0]


a6.pca <- prcomp(a6[,c(2:33)], 
                   center = TRUE, 
                   scale. = TRUE) 

# summary of the  
# prcomp object 
summary(a6.pca)
# loading library 
library(ggfortify) 
icu_nonICU.plot <- autoplot(a6.pca, 
                            data = a6, 
                            colour = 'ICU') 
icu_nonICU.plot
ggsave("~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/pca_icuVSnonicu.pdf", 
       icu_nonICU.plot, dpi = 300, width = 10, height = 10) 

# Load necessary libraries
library(ggfortify)
library(ggplot2)

# Generate the PCA plot with larger fonts and enhanced style
icu_nonICU.plot <- autoplot(a6.pca, 
                            data = a6, 
                            colour = 'ICU',
                            size = 4) +
  theme_minimal() +  # Cleaner background
  theme(
    text = element_text(size = 14),  # Increase font size
    axis.title = element_text(size = 16, face = "bold"),  # Axis title size and boldness
    axis.text = element_text(size = 14),  # Axis text size
    legend.title = element_text(size = 14),  # Legend title size
    legend.text = element_text(size = 12),  # Legend text size
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # Centered and bold title
  ) +
  labs(
    #title = "PCA Plot: ICU vs Non-ICU",  # Add a title
    #x = "Principal Component 1",         # Customize x-axis label
    #y = "Principal Component 2",         # Customize y-axis label
    colour = "Group"                     # Customize legend title
  ) +
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF"))  # Customize color palette

# Print the plot
print(icu_nonICU.plot)
ggsave("~/Desktop/phd/Thesis/ch2/cibersortx/pca_icuVSnonicu_11.pdf", 
       icu_nonICU.plot, dpi = 300, width = 5, height = 5) 


## PCA with 67 samples with new metadata file ##
a7 <- read.csv("~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/corrected_meta_cibersortx_67.csv", header = T)
a7.pca <- prcomp(a7[,c(4:41)], 
                 center = TRUE, 
                 scale. = TRUE) 
# summary of the  
# prcomp object 
summary(a7.pca)
icu_nonICU.plot <- autoplot(a7.pca, 
                            data = a7, 
                            colour = 'ICU') 
icu_nonICU.plot
ggsave("~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/pca_icuVSnonicu_67correctedMeta.pdf", 
       icu_nonICU.plot, dpi = 300, width = 10, height = 10) 



# library
library(ggplot2)
library(tidyr)
library(reshape)

a6 <- read.csv("~/Desktop/phd/Thesis/ch2/cibersortx/results/CIBERSORTx_Job7_Results-2.csv")


a6  <- a6[, colSums(a6 != 0) > 0]

data <- a6[ ,c(34,2:33)]
data <- data[order(data$ICU,decreasing=TRUE),]


# melt - to reshape and elongate the data frames in a user-defined manner
# cast is oposite to cast, it takes long-format data and casts it into wide-format data
dat.m = melt(data, id.var=c("ICU"))

# grouped boxplot
ggplot(dat.m, aes(x=ICU, y=value, fill=variable)) + 
  geom_boxplot()

# one box per variety
ggplot(dat.m, aes(x=ICU, y=value, fill=variable)) + 
  geom_boxplot() +
  facet_wrap(~variable, scale="free", ncol=5, strip.position = "left")

#save plot when you're happy with it
ggsave("~/Desktop/phd/Thesis/ch2/cibersortx/boxplots_11.pdf", width = 20,
       height = 30)

#########   CALCULATE STATS #########
#####################################

# upload the cell deconvoluted df
data <- a6[ ,c(34,2:33)]
data <- data[order(data$ICU,decreasing=TRUE),] # ICU group: Yes, No and the cell proportions

# Split the data into ICU and NON groups
icu_data <- subset(data, grepl("^ICU", ICU))
non_data <- subset(data, grepl("^non", ICU))

# Create an empty list to store results
results_list <- list()

# Perform Wilcoxon rank sum test for each variable and store results in the list
variables <- colnames(data[,2:33])

for (variable in variables) {
  test_result <- wilcox.test(icu_data[[variable]], non_data[[variable]])
  results_list[[variable]] <- c(variable, test_result$statistic, test_result$p.value)
}

# Convert the list to a data frame
results_df <- as.data.frame(do.call(rbind, results_list), stringsAsFactors = FALSE, row.names = 0)

# Set correct column names
colnames(results_df) <- c("Variable", "W_statistic", "p_value")

# Adjust p-values using Benjamini-Hochberg method
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")

# save the results data frame
write.csv(results_df, "~/Desktop/phd/Thesis/ch2/cibersortx/stats_11_wilcoxonText.csv")

ordered_stats <- results_df[order(results_df$p_value,decreasing=FALSE),]

# Filter dataframe to include only rows where the 'Name' column matches the name_vector
filtered_df <- dat.m %>%
  filter(variable %in% ordered_stats$Variable[1:10])

ggplot(filtered_df, aes(x=ICU, y=value, fill=variable)) + 
  geom_boxplot() +
  facet_wrap(~variable, scale="free", ncol=5, strip.position = "left")

filtered_df <- filtered_df %>%
  mutate(variable = factor(variable, levels = ordered_stats$Variable[1:10]))

# Create the plot with wrapped text in facet labels
library(ggplot2)
library(dplyr)
library(ggpubr)  # For adding statistical significance

# Assuming filtered_df is your dataframe
ggplot(filtered_df, aes(x = ICU, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +  # More compact boxplots and hiding outliers
  facet_wrap(~variable, scales = "free", ncol = 5, strip.position = "left") +
  theme_minimal(base_size = 14) +  # Clean theme with larger base font
  theme(
    strip.text = element_text(size = 8, face = "bold"),  # Larger facet labels with bold text
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1),  # Larger and rotated x-axis labels
    axis.text.y = element_text(size = 12),  # Larger y-axis labels
    axis.title.x = element_text(size = 16),  # Larger x-axis title
    axis.title.y = element_text(size = 16),  # Larger y-axis title
    legend.position = "none",  # Remove legend for compactness
    panel.spacing = unit(0.5, "lines"),  # Reduce spacing between panels
    strip.placement = "outside",  # Place facet labels outside the plot
    plot.margin = unit(c(0.1, 1, 1, 1), "cm"),  # Add margins around the plot
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)  # Add black border/frame around plots
  ) +
  labs(x = NULL, y = NULL) +
  scale_fill_brewer(palette = "Set3") +  # Use a nicer color palette
  stat_compare_means(
    aes(group = ICU),  # Grouping by ICU/non-ICU
    method = "wilcox.test",  # Choose appropriate test, e.g., t.test or wilcox.test
    label = "p.signif",  # Add stars for significance
    bracket.size = 10,
    hide.ns = FALSE  # Hide 'not significant' markers
  )

ggsave("~/Desktop/phd/Thesis/ch2/cibersortx/boxplots_11_withStats.pdf", width = 10,
       height = 10)


# CORRELATION ANALYSIS between gene expression and cell type proportions
# Inflammatory MФ; Interferone stimulating ciliated and secretory cells
# Let's assume your gene expression data is in a dataframe called gene_expression_df
# with columns "SampleID" and individual gene names (e.g., "Gene1", "Gene2", ...).
# And your Cibersortx results are in a dataframe called cibersortx_df with columns "SampleID"
# and cell type proportions for different cell types (e.g., "T_cell", "B_cell", ...).
# Make sure that the "SampleID" columns match in both dataframes, so you can merge them
# later based on this column.


library(tidyverse)
library(ggplot2)
library(ggfortify)
library(corrplot)
library(pheatmap)


########## Stats on cell proportions in ICU vs NON ########## 
#############################################################

## calculate the stats for cell proportions between ICU vs NON
# upload the cell deconvoluted df
data = read.csv("~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/corrected_meta_cibersortx_80.csv", header = T)
data <- data[ ,c(1:39)] # ICU group: Yes, No and the cell proportions
# Split the data into ICU and NON groups
icu_data <- subset(data, grepl("^ICU", SampleID))
non_data <- subset(data, grepl("^NON", SampleID))
# Create an empty list to store results
results_list <- list()
# Perform Wilcoxon rank sum test for each variable and store results in the list
variables <- colnames(data[,2:39])
for (variable in variables) {
  test_result <- wilcox.test(icu_data[[variable]], non_data[[variable]])
  results_list[[variable]] <- c(variable, test_result$statistic, test_result$p.value)
}
# Convert the list to a data frame
results_df <- as.data.frame(do.call(rbind, results_list), stringsAsFactors = FALSE, row.names = 0)
# Set correct column names
colnames(results_df) <- c("Variable", "W_statistic", "p_value")
# Adjust p-values using Benjamini-Hochberg method
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")
# Print the results data frame
print(results_df)
# save the results data frame
write.csv(results_df, "~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/DMRs_stats_80s_wilcoxonText.csv")
write.csv(correlation_results, "~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/DMRs_correlation80.csv")


##########################################################
##########################################################

save_pheatmap_pdf <- function(x, filename, width=9, height=11) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

### Batch normalized the data for 80 rnaseq samples
cts  <- read.csv("~/Desktop/rnaseq/last/raw/cts_total.csv", row.names = 1)
rownames_cts <- row.names(cts)
study <- read.csv("~/Desktop/rnaseq/last/raw/meta_1_2_3_4_5.csv", row.names=1) # dim is 80 x 19
col_by_icu <- row.names(study)
cts <- cts[col_by_icu]

## Remove low expression genes for batch analysis ##
## Do DEGs ##
raw <- edgeR::DGEList(counts=cts, samples=study)
sum(keep_raw <- edgeR::filterByExpr(raw, group=raw$samples$ICU, min.count = 1))  ## 60663 -> 29790
toDoBN <- raw[keep_raw, , keep.lib.sizes=TRUE]
# count matrix for BN
toDoBN_cts <- toDoBN$counts

# Check for bias in sequencing depth by KO types.
p<-ggplot(toDoBN$samples, aes(x=ICU, y=lib.size/1000000)) + 
  geom_violin(trim=FALSE)+geom_boxplot(width=0.1) + theme_minimal()+ ggtitle("Sequencing library size (millions)")

p <- p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)

df<-data.frame(colnames(toDoBN$counts), toDoBN$samples$ID, toDoBN$samples$lib.size)
# df[order(df$toDoBN.samples.lib.size, decreasing = F),][1:4,]
# ICU_62, ICU_13, ICU_58, NON_15  - < 10M
# S5054_13a, S5052_9, S5044_5a, S5049_144
write.csv(df, "~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/correlation/library_sizes_IDs.csv")
drops <- c("S5054_13a", "S5052_9", "S5044_5a", "S5049_144") # library sizes lower than 10M 
cts <- cts[, !(names(cts) %in% drops)]
study <- study[!(rownames(study) %in% drops),]
col_by_icu <- row.names(study)
cts <- cts[col_by_icu]

raw <- edgeR::DGEList(counts=cts, samples=study)
sum(keep_raw <- edgeR::filterByExpr(raw, group=raw$samples$ICU, min.count = 1))  ## 60663 -> 30106
toDoBN <- raw[keep_raw, , keep.lib.sizes=FALSE]
# count matrix for BN
toDoBN_cts <- toDoBN$counts

## Batch normalize using COMBAT-Seq ##
conditions = study$ICU
library_methods = study$Batch
replicates = study$Replicates
sample_names = names(cts)[1:length(names(cts))]

groups = sapply(as.character(conditions), switch, "No" = 1, "Yes" = 2)
batches = sapply(as.character(library_methods), switch, "Batch1" = 1, "Batch2" = 2,
                 "Batch3" = 3, "Batch4" = 4, "Batch5" =5, "Batch6"=6, USE.NAMES = F)
corrected_data = sva::ComBat_seq(counts = toDoBN_cts, batch = library_methods, group = groups)
write.csv(corrected_data, "~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/correlation/cts_batch_normalized_76.csv")

## add genes to the DGE list ##
ensg <- sub("\\..*", "", rownames(corrected_data))  # remove version number in case you have it
sym <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys=ensg,
                             column="SYMBOL", keytype="GENEID") # Unable to map 3533 of 60663 requested IDs.
gene <- data.frame(ENSGID=ensg, SYMBOL=sym, stringsAsFactors=F)
rownames(gene) <- rownames(corrected_data)
bn_dge <- edgeR::DGEList(counts=corrected_data, samples=study, genes = gene)
bn_dge <- edgeR::calcNormFactors(bn_dge, method="TMM") ## to estimate the effective library size

## DGEs ##
mm <- model.matrix( ~ 0 + ICU, data=bn_dge$samples)
colnames(mm) <- gsub("ICU", "", colnames(mm))
vm <- limma::voomWithQualityWeights(bn_dge, design=mm, plot=T)
cm <- limma::makeContrasts(
  ICUvsNON = Yes - No,
  levels= mm )
fit <- limma::lmFit(vm[ , rownames(mm) ], mm) # lmFit fits a linear model using weighted least squares for each gene
fit <- limma::contrasts.fit(fit, cm) # Estimate contrast for each gene
fit <- limma::eBayes(fit, trend=TRUE) 
summary(limma::decideTests(fit)) # Default p valut is 0.05.
tt <- topTable(fit, coef="ICUvsNON", adjust.method="fdr", n=Inf) %>% 
  rownames_to_column("ID") %>% 
  arrange(P.Value) %>% 
  mutate(sig0  = sign(logFC)*( abs(logFC) > log2(1.5) & adj.P.Val < 0.05 ),
         sig1  = sign(logFC)*( abs(logFC) > log2(2) & adj.P.Val < 0.05 ),
         sig2 = sign(logFC)*( abs(logFC) > log2(4) & adj.P.Val < 0.01 )) %>% 
  dplyr::select(ID, ENSGID, SYMBOL, 
                ICUvsNON_LFC=logFC, ICUvsNON_P=P.Value, ICUvsNON_FDR=adj.P.Val,
                ICUvsNON_sig0=sig0, ICUvsNON_sig1=sig1,ICUvsNON_sig2=sig2)




genesToSelectCorr <- dplyr::filter(tt, ICUvsNON_sig1 == 1 | ICUvsNON_sig1 == -1)

toDoBN_cts <- toDoBN$counts[genesToSelectCorr$ID,]

####
#gene_expression_df <- read.csv("~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/correlation/cts_batch_normalized_80.csv", header = T)
#gene_expression_df <- read.csv("~/Desktop/rnaseq/last/raw/80s_new_metadata/batch_normalized_67.csv", header = T)
gene_expression_df <- read.csv("~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/correlation/degs_67_bn.csv", header = T)
colnames(gene_expression_df)[1] <- "Geneid"
gene_expression_df <- gene_expression_df[,1:68]
rownames(gene_expression_df) <- gene_expression_df$Geneid
# Transpose data with setting First Column as Heading 
# setNames: making a dataframe with transposed gene_exp_df
gene_expression_df_t <- setNames(data.frame(t(gene_expression_df[ , - 1])), gene_expression_df[ , 1])

# Meta to rename my data file which is next
meta <- read.csv("~/Desktop/rnaseq/last/raw/80s_new_metadata/meta_80s.csv")
keep_samples <- colnames(gene_expression_df)[2:length(gene_expression_df)]
rownames(meta) <- meta$SampleID
meta <- meta[meta$SampleID %in% keep_samples, ]

col_by_icu <- row.names(meta)
gene_expression_df <- gene_expression_df[,-1]
gene_expression_df <- gene_expression_df[col_by_icu]
colnames(gene_expression_df) <- meta$AltName
gene_expression_df$Geneid <- rownames(gene_expression_df)
gene_expression_df <- gene_expression_df[,c(ncol(gene_expression_df),1:(ncol(gene_expression_df)-1))]
gene_expression_df_t <- setNames(data.frame(t(gene_expression_df[ , - 1])), gene_expression_df[ , 1])



# upload the cell deconvoluted df
data = read.csv("~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/corrected_meta_cibersortx_80.csv", header = T)
data <- data[ ,c(1:39)] # ICU group: Yes, No and the cell proportions

# Keep only significant cell types
keep <- c("SampleID","ITGAX.high.Macrophages","AZGP1.SCGB3A1.LTF.high.Goblet.Cells","VEGFA.high.Squamous.Cells",
          "SCGB1A1.high.Goblet.Cells","Interferon.Responsive.Cytotoxic.CD8.T.Cells","Dendritic.Cells",
          "Interferon.Responsive.Ciliated.Cells","B.Cells","Interferon.Responsive.Macrophages",
          "Inflammatory.Macrophages","Early.Response.FOXJ1.high.Ciliated.Cells","Ionocytes")
data <- data[, colnames(data) %in% keep]
keep_samples <- colnames(gene_expression_df)[2:length(gene_expression_df)]
data <- data[data$SampleID %in% keep_samples, ]

# Perform correlation analysis
correlation_results <- cor(gene_expression_df_t, data[,2:length(data)], method = "spearman", use = "complete.obs")
write.csv(correlation_results, "~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/DEGs/correlation_table_67.csv")

# calculate the p-values for each correlation
p_values <- matrix(NA, nrow = ncol(gene_expression_df_t), ncol = ncol(data)-1)

data <- data[,2:length(data)]

# Loop through the columns of the correlation matrix and perform correlation tests
for (i in 1:ncol(data)) {
  for (j in 1:ncol(gene_expression_df_t)) {
    correlation_test <- cor.test(gene_expression_df_t[, j], data[, i], method = "spearman", use = "complete.obs", exact = FALSE)
    p_values[j, i] <- correlation_test$p.value
  }
}

colnames(p_values) <- colnames(data)
rownames(p_values) <- colnames(gene_expression_df_t)

write.csv(p_values, "~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/DEGs/DEGs_stats_80s_wilcoxonText_p_values.csv")

df <- na.omit(correlation_results)
heatmap1<-pheatmap(df, show_rownames = F, cutree_rows = 2)

save_pheatmap_pdf(heatmap1, "~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/correlation/heatmap_bn_all_genes_corrected_sig_cellTypes_sigDEGs_67.pdf")

# Plots for correlation coefficients in histogram
correlation_results <- as.data.frame(correlation_results)
for (column in names(correlation_results)) {
  # Create a ggplot object for each column
  p <- ggplot(correlation_results, aes(x = !!sym(column))) +
    geom_density(fill = "lightblue", color = "black", alpha = 0.7, aes(y = ..scaled..)) +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "red") +
    geom_rug() +
    labs(title = paste("Density Plot of", column),
         x = "Correlation Coefficient",
         y = "Density") +
    theme_minimal()
  
  # Save each plot to a file (e.g., PNG)
  ggsave(paste0(column, "_density_plot_RNAseq.png"), plot = p, height = 5, width = 7)
}


#2 groups
three_clusters <- sort(cutree(heatmap1$tree_row, k=3))
print(table(three_clusters))
#three_clusters
#1     2     3 
#5641 11988  2996 
three_clusters <- t(data.frame(as.list(three_clusters)))
ID <- row.names(three_clusters)
three_clusters <- data.frame(ID, three_clusters)

write.csv(three_clusters, '~/Desktop/dnam/rnaseq/cybersortX/correlation_cybersortX/three_clusters_heatmap.csv')

#####################################
## Correlation DMRs and Cell Types ##
#####################################

dmr_betaValues_df <- read.delim("~/Desktop/dss/DMRichR/correlation_cell_type_DMRs/DMR_individual_smoothed_methylation_80s.txt", header = T)
dmr_betaValues_df <- dmr_betaValues_df[,c(1,18:97)]
# Transpose data with setting First Column as Heading 
# setNames: making a dataframe with transposed gene_exp_df
dmr_betaValues_df_t <- setNames(data.frame(t(dmr_betaValues_df[ , - 1])), dmr_betaValues_df[ , 1])

# upload the cell deconvoluted df
data = read.csv("~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/corrected_meta_cibersortx_80.csv", header = T)
data <- data[ ,c(1:39)] # ICU group: Yes, No and the cell proportions

# Keep only significant cell types
keep <- c("SampleID","ITGAX.high.Macrophages","AZGP1.SCGB3A1.LTF.high.Goblet.Cells","VEGFA.high.Squamous.Cells",
          "SCGB1A1.high.Goblet.Cells","Interferon.Responsive.Cytotoxic.CD8.T.Cells","Dendritic.Cells",
          "Interferon.Responsive.Ciliated.Cells","B.Cells","Interferon.Responsive.Macrophages",
          "Inflammatory.Macrophages","Early.Response.FOXJ1.high.Ciliated.Cells","Ionocytes")
data <- data[, colnames(data) %in% keep]

# Perform correlation analysis
correlation_results <- cor(dmr_betaValues_df_t, data[,2:length(data)], method = "spearman", use = "complete.obs")

df <- na.omit(correlation_results)

heatmap1<-pheatmap(df, show_rownames = F, cutree_rows = 2)

save_pheatmap_pdf(heatmap1, "~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/correlation/heatmap_bn_all_genes_corrected_sig_cellTypes_sigDMRs.pdf")

# Plot the correlation density plot
# Plots for correlation coefficients in histogram
correlation_results <- as.data.frame(correlation_results)
for (column in names(correlation_results)) {
  # Create a ggplot object for each column
  p <- ggplot(correlation_results, aes(x = !!sym(column))) +
    geom_density(fill = "lightblue", color = "black", alpha = 0.7, aes(y = ..scaled..)) +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "red") +
    geom_rug() +
    labs(title = paste("Density Plot of", column),
         x = "Correlation Coefficient",
         y = "Density") +
    theme_minimal()
  
  # Save each plot to a file (e.g., PNG)
  ggsave(paste0(column, "_density_plot.png"), plot = p, height = 5, width = 7)
}


# calculate the p-values for each correlation
p_values <- matrix(NA, nrow = ncol(dmr_betaValues_df_t), ncol = ncol(data)-1)

data <- data[,2:length(data)]

# Loop through the columns of the correlation matrix and perform correlation tests
for (i in 1:ncol(data)) {
  for (j in 1:ncol(dmr_betaValues_df_t)) {
    correlation_test <- cor.test(dmr_betaValues_df_t[, j], data[, i], method = "spearman", use = "complete.obs", exact = FALSE)
    p_values[j, i] <- correlation_test$p.value
  }
}

colnames(p_values) <- colnames(data)
rownames(p_values) <- colnames(dmr_betaValues_df_t)

write.csv(p_values, "~/Desktop/dnam/rnaseq/cybersortX/all_bRNA_samples_80/stats_80s_wilcoxonText_p_values.csv")



