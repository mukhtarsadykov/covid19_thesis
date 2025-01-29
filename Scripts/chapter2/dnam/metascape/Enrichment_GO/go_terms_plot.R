# Load required library
library(ggplot2)

a1 <- read.csv("~/Desktop/phd/Thesis/ch2/dnam/metascape/Enrichment_GO/GO_membership_parent1.csv")
# Create the data
# Load required library
library(ggplot2)

# Create the data
data <- data.frame(
  FDR_p_value = c(5.83036E-08, 7.6351E-07, 2.26526E-06, 2.54683E-06, 1.21677E-05, 
                  1.34306E-05, 2.59646E-05, 0.000114992, 0.000309698, 0.00050597, 
                  0.000632684, 0.000698949, 0.001627747, 0.003349074, 0.00673618),
  GO_List = c("cell morphogenesis",
              "regulation of plasma membrane bounded cell projection organization",
              "actin cytoskeleton organization",
              "import into cell",
              "regulation of GTPase activity",
              "enzyme-linked receptor protein signaling pathway",
              "negative regulation of cell population proliferation",
              "regulation of epithelial cell migration",
              "regulation of leukocyte differentiation",
              "protein phosphorylation",
              "neuromuscular process",
              "positive regulation of cell projection organization",
              "viral genome replication",
              "neuron projection extension",
              "regulation of pattern recognition receptor signaling pathway"),
  genes = c(28, 23, 25, 25, 22, 21, 21, 20, 22, 25, 19, 14, 19, 23, 19)
)

p1 <- ggplot(data, aes(x = genes, y = reorder(GO_List, genes), fill = FDR_p_value)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "darkblue", high = "lightblue") +
  labs(fill = "FDR\np-value") +  # Custom legend title
  xlab("Number of Genes") +
  ylab("GO Term") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),       # X-axis label size
    axis.title.y = element_text(size = 16),       # Y-axis label size
    axis.text.x = element_text(size = 14),        # X-axis text size
    axis.text.y = element_text(size = 14),        # Y-axis text size
    legend.title = element_text(size = 14),       # Legend title size
    legend.text = element_text(size = 12)         # Legend text size
  )

# Display the plot
pdf(file="~/Desktop/phd/Thesis/ch2/dnam/metascape/Enrichment_GO/GO_terms_numberGenes.pdf", height=8, width=12)
print(p1)
dev.off()
