library(ggplot2)
library(tidyr)
library(dplyr)

# Read the CSV
data <- read.csv("Desktop/phd/Thesis/analysis/Updated_APOBEC_Signature_and_Others_Data.csv")

# Convert Date column to Date class
data$Date <- as.Date(data$Date, format = "%m/%d/%y")  # Adjust date format to match your data

# Reshape data for ggplot
data_long <- gather(data, key = "signature", value = "Exposure", -Date)

# Plot the stacked area chart
ggplot(data_long, aes(x = Date, y = Exposure, fill = signature)) +
  geom_area(position = "fill") +  # Stacked and normalized to 100%
  labs(x = "Date", y = "Exposure", fill = "Signature") +
  scale_y_continuous(labels = scales::percent_format()) +  # Show y-axis as percentage
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +  # Show only month and year on x-axis
  theme_minimal() +  # Clean, minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

# Load necessary libraries
library(ggplot2)
library(tidyr)  # For the gather function

# Read the CSV
data <- read.csv("Desktop/phd/Thesis/analysis/Updated_APOBEC_Signature_and_Others_Data.csv")

# Convert Date column to Date class
data$Date <- as.Date(data$Date, format = "%m/%d/%y")  # Adjust date format to match your data

# Reshape data for ggplot
data_long <- gather(data, key = "signature", value = "Exposure", -Date)

# Reorder factor levels for the signature column to control the stacking order
data_long$signature <- factor(data_long$signature, levels = c("others_signatures", "APOBECs_signature"))

# Plot the stacked area chart
ggplot(data_long, aes(x = Date, y = Exposure, fill = signature)) +
  geom_area(position = "fill") +  # Stacked and normalized to 100%
  labs(x = "Date", y = "", fill = "Signature") +
  scale_y_continuous(labels = scales::percent_format()) +  # Show y-axis as percentage
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +  # Show only month and year on x-axis
  theme_minimal() +  # Clean, minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
