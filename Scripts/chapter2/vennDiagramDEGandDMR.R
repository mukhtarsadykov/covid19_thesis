library(VennDiagram)

# Assuming two_clusters is a data frame and has a column 'two_clusters'
degs <- read.csv("~/Desktop/phd/Thesis/ch2/dnam/metascape/promoter_exon_intron_5utr/hypo_hyper_up_down.csv")

hypo <- degs$Hypomethylated[1:12] 
hyper <- degs$Hypermethylated[1:160]
up <- degs$Upregulated
down <- degs$Downregulated[1:1686]

# Get the overlapping genes (intersection)
#overlap_genes <- intersect(hypo, hyper, up, down)

# Get the unique genes in correlation that are not in DGEs
#unique_to_correlation <- setdiff(correlation, DGEs)

# Get the unique genes in DGEs that are not in correlation
#unique_to_DEGs <- setdiff(DGEs, correlation)

# Generate a Venn diagram
# Generate a 4-set Venn diagram
venn.plot <- venn.diagram(
  x = list(
    Hypomethylated = hypo,
    Hypermethylated = hyper,
    Upregulated = up,
    Downregulated = down
  ),
  category.names = c("Hypo", "Hyper", "Up", "Down"), # Set labels
  filename = NULL,   # Don't save automatically, display in R
  output = TRUE,     # Plot the Venn diagram in the R plot window
  resolution = 300,  # High-resolution
  imagetype = "png", # Optional: can also be "jpeg", "tiff"
  fill = c("red", "blue", "green", "yellow"), # Colors for each circle
  alpha = 0.5,       # Transparency of the circles
  cex = 1.5,         # Font size for numbers
  cat.cex = 1.5,     # Font size for labels
  cat.pos = c(-20, 20, 120, -120), # Adjusts label positions
  cat.dist = c(0.05, 0.05, 0.05, 0.05), # Distance of labels from the circles
  lwd = 2,           # Line width of the circles
  margin = 0.1       # Adjusts the margin around the diagram
)

# Plot it
grid.draw(venn.plot)

# Intersection of hypo and hyper
intersect_hypo_hyper <- intersect(hypo, hyper)

# Intersection of hypo and up
intersect_hypo_up <- intersect(hypo, up)

# Intersection of hypo and down
intersect_hypo_down <- intersect(hypo, down)

# Similarly for other pairs
intersect_hyper_up <- intersect(hyper, up)
intersect_hyper_down <- intersect(hyper, down)
intersect_up_down <- intersect(up, down)

venn.plot <- venn.diagram(
  x = list(
    Hypomethylated = hypo,
    Upregulated = up
  ),
  category.names = c("Hypo", "Up"), # Set labels
  filename = NULL,   # Don't save automatically, display in R
  output = TRUE,     # Plot the Venn diagram in the R plot window
  resolution = 300,  # High-resolution
  imagetype = "png", # Optional: can also be "jpeg", "tiff"
 # fill = c("red", "blue", "green", "yellow"), # Colors for each circle
  alpha = 0.5,       # Transparency of the circles
  cex = 1.5,         # Font size for numbers
  cat.cex = 1.5,     # Font size for labels
  #cat.pos = c(-20, 20, 120, -120), # Adjusts label positions
  #cat.dist = c(0.05, 0.05, 0.05, 0.05), # Distance of labels from the circles
  lwd = 2,           # Line width of the circles
  margin = 0.1       # Adjusts the margin around the diagram
)
# Plot it
grid.draw(venn.plot)
