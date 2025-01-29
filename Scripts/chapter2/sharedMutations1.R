# Load necessary libraries
library(ggplot2)

meta_dna <- read.csv("~/Desktop/phd/meta_67s.csv")
# Load your data
data <- read.csv("~/Desktop/phd/rna_dna_mismatch/mix.csv")

m1 <- merge(data, meta_dna, by = "dnam")

write.csv(m1, "~/Desktop/m1.csv")

# Convert true_match to a factor
data$true_match <- as.factor(data$true_match)

# Perform a logistic regression to check the relationship between the predictors and true_match
model <- glm(true_match ~ shared_mutations + perc_match, data = data, family = binomial)

# Output the summary of the logistic regression
summary(model)

# Optionally, you can plot the relationship between true_match and shared_mutations/perc_match
ggplot(data, aes(x = shared_mutations, y = perc_match, color = true_match)) +
  geom_point() +
  labs(title = "Shared Mutations vs Percentage Match with True Match", x = "Shared Mutations", y = "Percentage Match")
