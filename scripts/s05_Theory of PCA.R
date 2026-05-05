# Theory of PCA ####


# Add your custom directory to the list of places R looks for gems
.libPaths(c("/courses/software/R_libs", .libPaths()))

# Now you can load them normally

library(tidyverse)
## Creating an example data set

# Create a vector for Cell IDs
cells <- c("Cell_1", "Cell_2", "Cell_3", "Cell_4")
# Create a vector to hold expression values for Gene A across all of the cells
Gene_A <- c(0, 12, 65, 23)
# Create a vector to hold expression values for Gene B across all of the cells
Gene_B <- c(4, 30, 57, 18)

# Create a tibble to hold the cell names and expression values
expression_tibble <- tibble(cells, Gene_A, Gene_B)
 
 # View the expression tibble
 expression_tibble
 
 
 # Let’s get an idea of what our data looks like by plotting it:
 
 # Create a plot to view the raw expression data
 ggplot(expression_tibble, aes(x = Gene_A, y = Gene_B, label = cells)) +
   geom_point(color = "cornflowerblue") +
   geom_text(hjust = 0, vjust = -1) +
   xlim(0, 80) +
   ylim(0, 80) +
   theme_bw() +
   xlab("Gene A") +
   ylab("Gene B") +
   ggtitle("Example Expression Values from Four Cells") +
   theme(plot.title = element_text(hjust = 0.5))
 
 ## Re-centering the dataset
 
 # Determine the center of the data by:
 # Finding the average expression of gene A
 Gene_A_mean <- mean(expression_tibble$Gene_A)
 # Finding the average expression of gene B
   Gene_B_mean <- mean(expression_tibble$Gene_B)
   
# Create a vector to hold the center of the data
center_of_data <- c(Gene_A_mean, Gene_B_mean)
# Assign names to the components of the vector
names(center_of_data) <- c("Gene_A", "Gene_B")
# Print out the center_of_data vector
  center_of_data

  
  # View where the center of the data is located
ggplot(expression_tibble, aes(x = Gene_A, y = Gene_B, label = cells)) +
  geom_point(color = "cornflowerblue") +
  annotate("point", x = Gene_A_mean, y = Gene_B_mean, color = "red", size = 3) +
  geom_text(hjust = 0, vjust = -1) +
  annotate("text", x = Gene_A_mean, y = Gene_B_mean, color = "red", label="New Center", hjust = 0, vjust = -1) +
  xlim(0, 80) +
  ylim(0, 80) +
  theme_bw() +
  xlab("Gene A") +
  ylab("Gene B") +
  ggtitle("Example Expression Values from Four Cells") +
  theme(plot.title = element_text(hjust = 0.5))

## 2. Translate the points from their raw expression coordinates to coordinates centered around the center of the data
# Shift the data points so that they data is centered on the origin
recentered_expression_tibble <- expression_tibble %>%
  mutate(
    Gene_A = Gene_A - Gene_A_mean,
    Gene_B = Gene_B - Gene_B_mean
    )
  
# Print out the re-centered data
recentered_expression_tibble
  
  
# Plot the raw data after it has been shifted to have the center of the data align with the origin
ggplot(recentered_expression_tibble, aes(x = Gene_A, y = Gene_B, label = cells)) +
  geom_point(color = "cornflowerblue") +
  annotate("point", x = 0, y = 0, color = "red", size = 3) +
  geom_text(hjust = 0, vjust = -1) +
  annotate("text", x = 0, y = 0, color = "red", label="New Center", hjust = 0, vjust = -1) +
  xlim(-50, 50) +
  ylim(-50, 50) +
  theme_bw() +
  xlab("Gene A") +
  ylab("Gene B") +
  ggtitle("Example Re-centered Expression Values from Four Cells") +
  theme(plot.title = element_text(hjust = 0.5))

