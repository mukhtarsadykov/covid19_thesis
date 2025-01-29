### ONLY for DEGs ###

####
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

gene_expression_df <- read.csv("~/Desktop/phd/Thesis/ch2/cibersortx/correlation/correlation_DEGs/counts_degs.csv", header = T)
colnames(gene_expression_df)[1] <- "Geneid"
gene_expression_df <- gene_expression_df[,1:12]
names(gene_expression_df) <- c("Geneid", "ICU_01", "ICU_02", "ICU_03", "ICU_04", "NON_01", "NON_02", "NON_03", "NON_04", "NON_05", "NON_06", "NON_07")
rownames(gene_expression_df) <- gene_expression_df$Geneid
# Transpose data with setting First Column as Heading 
# setNames: making a dataframe with transposed gene_exp_df
gene_expression_df_t <- setNames(data.frame(t(gene_expression_df[ , - 1])), gene_expression_df[ , 1])

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
write.csv(correlation_results, "~/Desktop/phd/Thesis/ch2/cibersortx/correlation/correlation_DEGs/correlation_table_DEGs_11.csv")

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

# Assuming your data frame is named pval_df
# Check the data types
sapply(p_values, class)
# Filter genes where any p-value is less than 0.05
filtered_pval_df <- p_values[apply(p_values, 1, function(x) any(x < 0.05)), ]
# Filter genes where all p-values are less than 0.05
filtered_pval_df <- p_values[apply(p_values, 1, function(x) all(x < 0.05)), ]


write.csv(p_values, "~/Desktop/phd/Thesis/ch2/cibersortx/correlation/correlation_DEGs/onlyDEGs_stats_11s_wilcoxonText_p_values.csv")

### NON-CORRELATION GENES in all cell types ###
# Filter genes where all absolute correlations are <= 0.3
filtered_genes <- correlation_results[apply(correlation_results, 1, function(row) all(abs(row) <= 0.5)), ]

# Define the minimum number of cell types with p-value < 0.05
N <- 7

# Filter genes where p-value is less than 0.05 in at least N cell types
filtered_genes <- correlation_results[apply(correlation_results, 1, function(row) sum(abs(row) <= 0.35) >= 6), ]
write.csv(filtered_genes, "~/Desktop/phd/Thesis/ch2/cibersortx/correlation/correlation_DEGs/non_correlation_genes_0.35.csv")

non_corr <- rownames(filtered_genes)

tt_sig <- read.csv("~/Desktop/phd/Thesis/ch2/topTable_sig1_11.csv", row.names = 1)

i <- which(tt_sig$ID %in% non_corr)
non_corr_genes <- tt_sig[i,]
write.csv(non_corr_genes, "~/Desktop/phd/Thesis/ch2/cibersortx/correlation/correlation_DEGs/non_correlation_genes_tt_table.csv")

df <- na.omit(correlation_results)
heatmap1<-pheatmap(df, show_rownames = F, cutree_rows = 2)

save_pheatmap_pdf <- function(x, filename, width=9, height=11) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(heatmap1, "~/Desktop/phd/Thesis/ch2/cibersortx/correlation/heatmap_DEG_genes_corrected_sig_cellTypes_sigDEGs_11.pdf")

# Plots for correlation coefficients in histogram
correlation_results <- as.data.frame(correlation_results)
for (column in names(correlation_results)) {
  # Create a ggplot object for each column
  p <- ggplot(correlation_results, aes(x = !!sym("Mast.Cells"))) +
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
  ggsave(paste0(column, "_density_plot_RNAseq_DEGs.png"), plot = p, height = 5, width = 7)
}

#2 groups
two_clusters <- sort(cutree(heatmap1$tree_row, k=2))
print(table(two_clusters))
#three_clusters
#1     2  
#2434 1686   
two_clusters <- t(data.frame(as.list(two_clusters)))
ID <- row.names(two_clusters)
two_clusters <- data.frame(ID, two_clusters)

write.csv(two_clusters, '~/Desktop/phd/Thesis/ch2/cibersortx/correlation/correlation_DEGs/two_clusters_heatmap_11_degs.csv')

library(VennDiagram)

# Assuming two_clusters is a data frame and has a column 'two_clusters'
degs <- read.csv("~/Desktop/phd/Thesis/ch2/two_clusters_heatmap11.csv")
rows_with_1_corr <- two_clusters[two_clusters$two_clusters == 1, ] # from correlation analysis
rows_with_1_degs <- degs[degs$two_clusters==1, ] # from DEGs

# Define the sets
correlation <- rows_with_1_corr$ID
DGEs <- rows_with_1_degs$ID

# Get the overlapping genes (intersection)
overlap_genes <- intersect(correlation, DGEs)

# Get the unique genes in correlation that are not in DGEs
unique_to_correlation <- setdiff(correlation, DGEs)

# Get the unique genes in DGEs that are not in correlation
unique_to_DEGs <- setdiff(DGEs, correlation)

# Generate a Venn diagram
venn.plot <- venn.diagram(
  x = list(Correlation = correlation, Original = DGEs), # Define set names
  category.names = c("Cell proportions", "DEGs"), # Set labels
  filename = NULL, # Don't save it automatically, keep it in R
  output = TRUE,   # Output to R's plot window
  resolution = 300, # High-resolution image setting
  imagetype = "png", # Image type (optional)
  fill = c("blue", "green"), # Colors for the circles
  alpha = 0.5, # Transparency of circles
  cex = 2, # Size of the text
  cat.cex = 2, # Size of the category names
  cat.pos = c(-20, 20) # Adjusts label positions
)

# Plot it
grid.draw(venn.plot)

library(clusterProfiler)

cluster1_larger <- dplyr::filter(two_clusters, two_clusters=="1")
cluster1_genes <- cluster1_larger$ID

cluster1 <- enrichGO(
  keyType = "ENSEMBL",
  gene = unique_to_correlation,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

cluster2_smaller <- dplyr::filter(two_clusters, two_clusters=="2")
cluster2_genes <- cluster2_smaller$ID

cluster2 <- enrichGO(
  keyType = "ENSEMBL",
  gene = unique_to_DEGs,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

write.csv(cluster1@result, "~/Desktop/phd/Thesis/ch2/clusterProfiler_cluster1.csv")
write.csv(cluster2@result, "~/Desktop/phd/Thesis/ch2/clusterProfiler_cluster2.csv")

p1 <- dotplot(cluster1, showCategory = 20 ) + viridis::scale_fill_viridis()
p2 <- dotplot(cluster2, showCategory = 20 ) + viridis::scale_fill_viridis()

pdf(file="~/Desktop/phd/Thesis/ch2/clusterProfiler_heatmapClusters_selected_11.pdf", height=10, width=13)
cowplot::plot_grid(p1, p2, ncol=2) 
dev.off()



