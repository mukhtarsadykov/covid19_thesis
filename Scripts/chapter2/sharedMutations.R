# Load necessary libraries
library(ggplot2)

# Load your data
data <- read.csv("~/Desktop/m1.csv")

# Convert true_match and Batch_dna to factors
data$true_match <- as.factor(data$true_match)
data$Batch_dna <- as.factor(data$Batch_dna)

n_shapes <- length(unique(data$Batch_dna))
shapes <- c(0:13)

# Perform a logistic regression to check the relationship between the predictors and true_match
model <- glm(true_match ~ shared_mutations + Batch_dna, data = data, family = binomial)

# Output the summary of the logistic regression
summary(model)

# Plot with custom shapes
ggplot(data, aes(x = shared_mutations, y = perc_match, color = true_match, shape = Batch_dna)) +
  geom_point(size=3) +  # Adjust size as needed
  scale_shape_manual(values = shapes[1:n_shapes]) +  # Assign unique shapes
  labs(title = "Shared Mutations vs Percentage Match with True Match and Batch DNA",
       x = "Shared Mutations", y = "Percentage Match") +
  theme_minimal()

ggplot(data, aes(x = shared_mutations, y = perc_match, color = true_match)) +
  geom_point(size=0.5) + 
  facet_wrap(~ Batch_dna) +  # Create a small plot for each Batch_dna
  labs(title = "Shared Mutations vs Percentage Match by Batch DNA",
       x = "Shared Mutations", y = "Percentage Match") +
  theme_minimal()
