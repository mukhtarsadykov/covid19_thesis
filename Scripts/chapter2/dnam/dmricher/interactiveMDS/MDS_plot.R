# Load necessary libraries
library(ggplot2)
library(ggrepel)

# Create the data frame with your provided data
df <- data.frame(
  dim1 = c(-0.2462, -0.2806, -0.2318, -0.1771, 0.2666, 0.0451, 0.1371, 0.12, 0.09587, -0.009653, 0.2807),
  dim2 = c(-0.2277, 0.2171, 0.1357, 0.02414, 0.08596, -0.05162, 0.08547, 0.00598, 0.03366, -0.3134, 0.004655),
  dim3 = c(-0.132, -0.1909, 0.3878, -0.06848, 0.003603, 0.001156, -0.06461, -0.01907, -0.04868, 0.1043, 0.02681),
  dim4 = c(0.2858, -0.175, 0.07799, -0.03482, 0.1244, -0.05819, 0.04795, -0.01423, 0.004965, -0.241, -0.01788),
  dim5 = c(0.01539, -0.02849, 0.009801, -0.05197, -0.2001, 0.3213, -0.01984, -0.009724, -0.002064, -0.1427, 0.1084),
  dim6 = c(-0.06034, -0.1615, -0.001811, 0.2088, -0.1978, -0.09017, 0.1403, 0.09928, 0.1433, -0.0372, -0.04301),
  dim7 = c(0.06985, 0.07437, 0.02639, -0.09164, -0.1503, -0.201, -0.0122, -0.003119, 0.01053, -0.001089, 0.2782),
  dim8 = c(0.02248, 0.04545, 0.02596, -0.2391, -0.06757, -0.00442, 0.1443, 0.09059, 0.07881, 0.03771, -0.1342),
  label = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
  Name = c("ICU_01",	"ICU_02",	"ICU_03",	"ICU_04",	"NON_01",	"NON_02",	"NON_03",	"NON_04",	"NON_05",	"NON_06",	"NON_07"),
  Group = c("ICU_ABS", "ICU_ABS", "ICU_ABS", "ICU_ABS", "NONICU", "NONICU", "NONICU", "NONICU", "NONICU", "NONICU", "NONICU"),
  Age = c(62, 47, 63, 45, 72, 62, 45, 44, 65, 50, 70),
  Sex = c("Male", "Female", "Female", "Male", "Male", "Female", "Female", "Male", "Male", "Female", "Male"),
  Batch = c("Batch4", "Batch6", "Batch7", "Batch11", "Batch1", "Batch7", "Batch3", "Batch3", "Batch3", "Batch4", "Batch4")
)

# Create 'ICU' variable based on 'Group'
df$ICU <- ifelse(df$Group == "ICU_ABS", "Yes", "No")

# Eigenvalues (variance explained)
eigenvalues <- c(0.23, 0.13, 0.12, 0.11, 0.10, 0.09, 0.09, 0.07)
total_variance <- sum(eigenvalues)
var_exp <- eigenvalues / total_variance

# Create the MDS plot with specified format and style
p1 <- ggplot(df, aes(x = dim1, y = dim2, label = Name, col = ICU)) +
  geom_point(size = 7, alpha = 0.7) +  # Increase point size and add transparency
  geom_text_repel(
    box.padding = 0.6,             # Increase padding around the text box
    point.padding = 0.5,           # Increase padding around the data point
    min.segment.length = 0.3,      # Adjust the length of the segment lines
    size = 6                       # Increase text size
  ) +
  scale_color_manual(values = c("Yes" = "#ECA39A", "No" = "#74CFD4")) +
  labs(
    x = paste0("Component 1 (var. explained ", round(var_exp[1] * 100, 2), "%)"),
    y = paste0("Component 2 (var. explained ", round(var_exp[2] * 100, 2), "%)")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20),    # Increase axis text size
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 22),   # Increase axis title size
    axis.title.y = element_text(size = 22),
    axis.line = element_line(size = 1.2),     # Increase axis line thickness
    axis.ticks = element_line(size = 1.2),    # Thicken axis ticks
    legend.position = "right",
    legend.text = element_text(size = 20),    # Increase legend text size
    legend.title = element_text(size = 20)    # Increase legend title size
  )

# Display the plot
pdf(file="dmricher/interactiveMDS/MDS_11_newNames.pdf", height=8, width=12)
print(p1)
dev.off()
