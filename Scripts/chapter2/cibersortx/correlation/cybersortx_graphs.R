# library
library(ggplot2)

data = read.csv("~/Desktop/dnam/rnaseq/cybersortX/tmpp/toPlot copy.csv", header = T)

# melt - to reshape and elongate the data frames in a user-defined manner
# cast is oposite to cast, it takes long-format data and casts it into wide-format data
dat.m = melt(data, id.var=c("Group"))

# grouped boxplot
ggplot(dat.m, aes(x=Group, y=value, fill=variable)) + 
  geom_boxplot()

# one box per variety
ggplot(dat.m, aes(x=Group, y=value, fill=variable)) + 
  geom_boxplot() +
  facet_wrap(~variable, scale="free", strip.position = "left")

#save plot when you're happy with it
ggsave("~/Desktop/dnam/rnaseq/cybersortX/tmpp/boxplots.pdf", width = 10,
       height = 20)


# Excluded
#ICU_22
#ICU_10
#ICU_38
#ICU_23
#ICU_19
#ICU_31

# calculate stats

df <- data.frame(Cell= character(0), p_value = numeric(0), ci_l = numeric(0), ci_u = numeric(0))

for (i in 2:31){
  a <- t.test(data[1:14,i], data[c(17:30),i], paired = T)
  output <- c(colnames(data)[i], a[["p.value"]], a[["conf.int"]][1], a[["conf.int"]][2])
  df <- rbind(df, output)
}

write.csv(df, "~/Desktop/dnam/rnaseq/cybersortX/tmpp/stats_30s.csv")

t.test(data[1:14,2], data[c(17:30),2], paired = T) #developing_ciliated_cells: p-value = 0.007657
t.test(data[1:14,3], data[c(17:30),3], paired = T) #SPRR2D high Squamous Cells: p-value = 0.1293
t.test(data[1:14,4], data[c(17:30),4], paired = T) #KRT24 KRT13 high Secretory Cells: p-value = 0.05069
t.test(data[1:14,5], data[c(17:30),5], paired = T) #ITGAX high Macrophages: p-value = 1.992e-07
t.test(data[1:14,6], data[c(17:30),6], paired = T) #B Cells: p-value = 0.009354
t.test(data[1:14,7], data[c(17:30),7], paired = T) #ITGAX_macrophages_cells: p-value = 0.000464
t.test(data[1:14,8], data[c(17:30),8], paired = T) #MSR1C_macrophages_cells: p-value = 0.2018
t.test(data[1:14,9], data[c(17:30),9], paired = T) #dendritic_cells: p-value = 0.05443

t.test(data[1:14,13], data[c(15,18:25,28:30,32),13], paired = T) #IFN_responsive_secretory_cells: p-value = 0.01767

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
library(corrplot)

save_pheatmap_pdf <- function(x, filename, width=9, height=11) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

gene_expression_df <- read.csv("~/Desktop/dnam/rnaseq/cybersortX/correlation_cybersortX/cts_batch_normalized.csv", header = T)
colnames(gene_expression_df)[1] <- "Group"

# Transpose data with setting First Column as Heading 
gene_expression_df_t <- setNames(data.frame(t(gene_expression_df[ , - 1])), gene_expression_df[ , 1])

# Perform correlation analysis
correlation_results <- cor(gene_expression_df_t, data, method = "pearson", use = "complete.obs")
df <- na.omit(correlation_results[2:20626,])
heatmap1<-pheatmap(df, show_rownames = F, cutree_rows = 3)

save_pheatmap_pdf(heatmap1, "~/Desktop/dnam/rnaseq/cybersortX/correlation_cybersortX/heatmap_bn_3clusters.pdf")

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





