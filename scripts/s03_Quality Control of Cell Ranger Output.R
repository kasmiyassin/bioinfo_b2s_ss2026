# 3. Quality Control of Cell Ranger Output ----
## Metrics evaluation ----


# Add your custom directory to the list of places R looks for gems
.libPaths(c("/courses/software/R_libs", .libPaths()))

# Now you can load them normally

# Load libraries
library(tidyverse)
library(reshape2)
library(plyr)

# Names of samples (same name as folders stored in data)
samples <- c("ctrl", "stim")

# Loop over each sample and read the metrics summary in
metrics <- list()
for (sample in samples) {
  path_csv <- paste0("data/", sample, "_metrics_summary.csv")
  df <- read.csv(path_csv)
  rownames(df) <- sample
  df$sample <- sample
  metrics[[sample]] <- df
}
# Concatenate each sample metrics together
metrics <- ldply(metrics, rbind)

# Remove periods and percentages to make the values numeric
metrics <- metrics %>%
  column_to_rownames(".id") %>%
  mutate(across(everything(), parse_number(str_replace(., ",", "")))) %>%
  mutate(across(everything(), parse_number(str_replace(., "%", ""))))

#  columns_in_metrics_data
colnames(metrics)


# Mapping percentage of reads to different regions of the genome.
# Columns of interest
cols <- c("Reads.Mapped.Confidently.to.Intergenic.Regions",
          "Reads.Mapped.Confidently.to.Intronic.Regions",
          "Reads.Mapped.Confidently.to.Exonic.Regions",
          "sample")

# Data wrangling to sculpt dataframe in a ggplot friendly manner
df <- metrics %>%
  select(all_of(cols)) %>%
  melt(id.vars="sample") %>%
  mutate(variable = str_replace_all(variable, "Reads.Mapped.Confidently.to.", "")) %>%
  mutate(variable = str_replace_all(variable, ".Regions", "")) %>%
  mutate(value = str_replace_all(value, "%", "")) %>% # Remove % from value
  mutate(value = as.numeric(value)) # Make percentage numeric

# ggplot code to make a barplot
df %>% ggplot() +
  geom_bar(
    aes(x = sample, y = value, fill = variable),
    position = "stack",
    stat = "identity") +
  coord_flip() +
  labs(
    x = "Sample",
    y = "Percentage of Reads",
    title = "Percent of Reads Mapped to Each Region",
    fill = "Region")


