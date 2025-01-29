rm(list=ls())
pacman::p_load(limma, Glimma, edgeR, AnnotationDbi,org.Hs.eg.db, EnsDb.Hsapiens.v86, tidyverse, ggplot2, ... = gridExtra, ggrepel, reshape2, EnhancedVolcano, GGally, sva, pheatmap, clusterProfiler, viridis, devtools, GeneSetCluster, readxl, gplots, ggpubr, wesanderson, factoextra)

save_pheatmap_pdf <- function(x, filename, width=9, height=11) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# upload the data
cts  <- read.csv("~/Desktop/phd/Thesis/ch2/cts_total.csv", row.names = 1)
study <- read.csv("~/Desktop/phd/Thesis/ch2/meta_11.csv", row.names=1) # dim is 80 x 19
study <- study[1:11,]
study$AltName <- c("ICU_01", "ICU_02", "ICU_03", "ICU_04", "NON_01", "NON_02",
                   "NON_03", "NON_04", "NON_05", "NON_06", "NON_07")
cts <- cts[,row.names(study)]
names(cts) <- study$AltName
row.names(study) <- study$AltName
stopifnot( identical( colnames(cts), rownames(study) ) )

## Remove low expression genes for batch analysis ##
## Do DEGs ##
ensg <- sub("\\..*", "", rownames(raw_cts))  # remove version number in case you have it
sym <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys=ensg,
                             column="SYMBOL", keytype="GENEID") # Unable to map 1023 of 27147 requested IDs.
gene <- data.frame(ENSGID=ensg, SYMBOL=sym, stringsAsFactors=F)
rownames(gene) <- rownames(raw_cts)

raw <- edgeR::DGEList(counts=raw_cts, samples=study, genes=gene)
sum(keep_raw <- edgeR::filterByExpr(raw, group=raw$samples$ICU, min.count = 1))  ## 60663 -> 30565
raw <- raw[keep_raw, , keep.lib.sizes=TRUE]

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
gene_expression_df <- read.csv("~/Desktop/phd/Thesis/ch2/batch_normalized_11.csv", header = T)
colnames(gene_expression_df)[1] <- "Geneid"
gene_expression_df <- gene_expression_df[,1:12]
names(gene_expression_df) <- c("Geneid", "ICU_01", "ICU_02", "ICU_03", "ICU_04", "NON_01", "NON_02", "NON_03", "NON_04", "NON_05", "NON_06", "NON_07")
rownames(gene_expression_df) <- gene_expression_df$Geneid
# Transpose data with setting First Column as Heading 
# setNames: making a dataframe with transposed gene_exp_df
gene_expression_df_t <- setNames(data.frame(t(gene_expression_df[ , - 1])), gene_expression_df[ , 1])


# upload the cell deconvoluted df
a6 <- read.csv("~/Desktop/phd/Thesis/ch2/cibersortx/results/CIBERSORTx_Job7_Results-2.csv", header = T)

data = read.csv("~/Desktop/phd/Thesis/ch2/cibersortx/results/CIBERSORTx_Job7_Results-2.csv", header = T)
data <- data[1:11, ]
data  <- data[, colSums(data != 0) > 0]
data <- data[ ,c(1:33)] # ICU group: Yes, No and the cell proportions
data$Mixture <- c("ICU_01", "ICU_02", "ICU_03", "ICU_04", "NON_01", "NON_02", "NON_03", "NON_04", "NON_05", "NON_06", "NON_07")
names(data)[1] <- "SampleID"
# Keep only significant cell types
keep <- c("SampleID","ITGAX.high.Macrophages","AZGP1.SCGB3A1.LTF.high.Goblet.Cells","VEGFA.high.Squamous.Cells",
          "Interferon.Responsive.Ciliated.Cells", "Basal.Cells", "Mast.Cells",
          "Inflammatory.Macrophages", "Interferon.Responsive.Secretory.Cells", "HOPX.high.Squamous.Cells", 
          "Interferon.Responsive.Cytotoxic.CD8.T.Cells")

data <- data[, colnames(data) %in% keep]
keep_samples <- colnames(gene_expression_df)[2:length(gene_expression_df)]
data <- data[data$SampleID %in% keep_samples, ]
row.names(data) <- data$SampleID
data <- data[,2:length(data)]
data <- data.frame(lapply(data, function(x) as.numeric(as.character(x))))



# Perform correlation analysis
correlation_results <- cor(gene_expression_df_t, data, method = "spearman", use = "complete.obs")
write.csv(correlation_results, "~/Desktop/phd/Thesis/ch2/cibersortx/correlation/correlation_all_genes/correlation_table_11.csv")

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

write.csv(p_values, "~/Desktop/phd/Thesis/ch2/cibersortx/correlation/correlation_all_genes/Allgenes_stats_11s_wilcoxonText_p_values.csv")

df <- na.omit(correlation_results)
heatmap1<-pheatmap(df, show_rownames = F, cutree_rows = 2)

save_pheatmap_pdf(heatmap1, "~/Desktop/phd/Thesis/ch2/cibersortx/correlation/heatmap_all_genes_corrected_sig_cellTypes_11.pdf")

# Plots for correlation coefficients in histogram
correlation_results <- as.data.frame(correlation_results)
for (column in names(correlation_results)) {
  # Create a ggplot object for each column
  p <- ggplot(correlation_results, aes(x = !!sym(column))) +
    geom_density(fill = "lightblue", color = "black", alpha = 0.7, aes(y = ..scaled..), linewidth = 2) +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "red", linewidth = 2) +
    geom_rug() +
    theme_minimal() +
    labs(
      #x = "Correlation Coefficient",
      y = "Density") +
    theme(
      axis.text.x = element_text(size = 20),  # Increase axis text size
      axis.text.y = element_text(size = 20),
      axis.title.y = element_text(size = 22),                                      # Increase axis title size
      axis.title.x = element_text(size = 22),
      axis.line = element_line(size = 1.2),                                        # Increase axis line thickness
      axis.ticks = element_line(size = 1.2),                                       # Thicken axis ticks
      #legend.position = "right",
      legend.text = element_text(size = 20),                                       # Increase legend text size
      legend.title = element_text(size = 20)                                       # Increase legend title size
    )
  
  # Save each plot to a file (e.g., PNG)
  ggsave(paste0(column, "_density_plot_RNAseq_allGenes.png"), plot = p, height = 5, width = 7)
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



