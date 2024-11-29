library(ggplot2)

a1 <- readxl::read_excel("~/Desktop/phd/Thesis/trinucleotide_mutations.xlsx")

pdf("~/Desktop/phd/Thesis/trinucleotide_mutations.pdf", width = 20, height = 10)
p <- ggplot(a1, aes(x=codons, y=APOBEC_signature, fill=mutation)) + 
  geom_bar(stat="identity", position=position_dodge())
p + scale_fill_brewer(palette="Paired") + 
  theme_minimal() + 
  labs(x = "codons", y = "prop. of mutations", fill = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # Rotate x-axis labels and increase text size
        axis.text.y = element_text(size = 16),  # Increase y-axis text size
        axis.title.x = element_text(size = 16),  # Increase x-axis title text size
        axis.title.y = element_text(size = 16),  # Increase y-axis title text size
        axis.ticks.x = element_line(),
        legend.position = "top",  # Move legend to the top
        legend.text = element_text(size = 16),  # Increase legend text size
        legend.title = element_text(size = 16)) + 
  guides(fill = guide_legend(nrow = 2)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1))
dev.off()

pdf("~/Desktop/phd/Thesis/trinucleotide_mutations_separate.pdf", width = 15, height = 10)
p <- ggplot(a1, aes(x=codons, y=APOBEC_signature, fill=mutation)) + 
  geom_bar(stat="identity", position=position_dodge())

p + scale_fill_brewer(palette="Paired") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
  facet_wrap(~mutation)
dev.off()


library(dplyr)

pdf("~/Desktop/phd/Thesis/trinucleotide_mutations_separate_nozero.pdf", width = 15, height = 10)
p <- ggplot(a1 %>% filter(APOBEC_signature != 0), aes(x=codons, y=APOBEC_signature, fill=mutation)) + 
  geom_bar(stat="identity", position=position_dodge())

p + scale_fill_brewer(palette="Paired") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
  facet_wrap(~mutation)
dev.off()
theme_minimal() +  # Clean, minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # Rotate x-axis labels and increase text size
    axis.text.y = element_text(size = 16),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis title text size
    axis.title.y = element_text(size = 16),  # Increase y-axis title text size
    axis.ticks.x = element_line(),
    legend.position = "top",  # Move legend to the top
    legend.text = element_text(size = 16),  # Increase legend text size
    legend.title = element_text(size = 16)) # Increase legend title size

# Load necessary libraries
library(ggplot2)
library(ggalluvial)
library(tidyverse)
library(tidyr)  # For the gather function

data <- read.csv("Desktop/phd/Thesis/analysis/Updated_APOBEC_Signature_and_Others_Data.csv")
data$Date <- as.Date(data$Date, format = "%m/%d/%y")  # Adjust date format to match your data
data_long <- gather(data, key = "signature", value = "Exposure", -Date)

# Reorder factor levels for the signature column to control the stacking order
data_long$signature <- factor(data_long$signature, levels = c("others_signatures", "APOBECs_signature"))

# Plot the stacked area chart
pdf("~/Desktop/phd/Thesis/apobecs_signature_temporal.pdf", width = 20, height = 10)
ggplot(data_long, aes(x = Date, y = Exposure, fill = signature)) +
  geom_area(position = "fill") +  # Stacked and normalized to 100%
  labs(x = "", y = "", fill = "") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +  # Show y-axis as percentage
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months", expand = c(0, 0)) +  # Show only month and year on x-axis
  theme_minimal() +  # Clean, minimal theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),  # Rotate x-axis labels and increase text size
    axis.text.y = element_text(size = 16),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis title text size
    axis.title.y = element_text(size = 16),  # Increase y-axis title text size
    axis.ticks.x = element_line(),
    legend.position = "top",  # Move legend to the top
    legend.text = element_text(size = 16),  # Increase legend text size
    legend.title = element_text(size = 16)  # Increase legend title size
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) # Center legend title
dev.off()



a1 <- readxl::read_excel("~/Desktop/phd/Thesis/trinucleotide_mutations.xlsx")
a1 <- a1[,1:4]

##Go through each row and determine if a value is zero
numeric_cols <- sapply(a1, is.numeric) # Identify numeric columns
tmp1 <- a1[!apply(a1[, numeric_cols], 1, function(row) all(row == 0)), ]
##Subset as usual
names(tmp1)[2] <- "mutation"
names(tmp1)[3] <- "prop_of_mutations"
names(tmp1)[4] <- "source"

pdf("~/Desktop/phd/Thesis/apobecs_signature_alluvial_full.pdf", width = 10, height = 16)
ggplot(tmp1,
       aes(y = prop_of_mutations, axis1 = codons, axis2 = mutation)) +
  geom_alluvium(aes(fill = source), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Trinucleotide", "Mutation"), expand = c(.2, .2)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),       # Hide y-axis text
    axis.title.x = element_text(size = 16),
    axis.title.y = element_blank(),       # Hide y-axis title
    axis.ticks.y = element_blank(),       # Hide y-axis ticks
    axis.ticks.x = element_line(),
    legend.position = "top",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )
dev.off()

# frequency table

a2 <- read.csv("~/Desktop/phd/Thesis/Trinucleotide_Frequencies_in_SARS-CoV-2_Genome.csv")
# lock in factor level order
a2$codons <- factor(a2$codons, levels = a2$codons)

# Change barplot line colors by groups

pdf("~/Desktop/phd/Thesis/trinucleotide_frequencies_with_sources.pdf", width = 20, height = 10)
ggplot(a2, aes(x=codons, y=count, fill=source)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#C0E0D4", "#F7CBB8", "#CBD2E6", "#D9D9D9"))+
  theme_minimal()+
  labs(x = "trinucleotides", y = "counts", fill = "sources: ") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 2.2, vjust = 0.5, size = 16),
    axis.text.y = element_text(size = 16),       # Hide y-axis text
    axis.title.y = element_text(size = 16),       # Hide y-axis title
    axis.ticks.y = element_line(),       # Hide y-axis ticks
    #axis.ticks.x = element_line(),
    axis.title.x = element_text(size = 16),
    legend.position = "top",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )
dev.off()


