# Install required packages
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("SNPRelate", "vcfR"))

library(SNPRelate)


# Close all open GDS files to reset the handles
snpgdsClose(genofile) 

# If that variable is lost, use this nuclear option:
showfile.gds(closeall=TRUE)

# Now try opening it again
genofile <- snpgdsOpen("~/vcf_treamt/multisample_300.gds")

setwd("~/vcf_treamt")
# Define paths (referencing your session02vc directory)
vcf_fn <- "/courses/master_b2s/raw_data/s02_vc/treatment/treatment_02.vcf"
gds_fn <- "~/vcf_treamt/multisample_300.gds"


# Corrected Conversion Command
snpgdsVCF2GDS(vcf_fn, gds_fn, method="biallelic.only")
# Open the newly created GDS file
genofile <- snpgdsOpen(gds_fn)

# Convert VCF to GDS (Genomic Data Structure)
#snpgdsVCF2GDS(vcf_fn, gds_fn, method="ALL")
#genofile <- snpgdsOpen(gds_fn)

# Calculate IBS Matrix
ibs <- snpgdsIBS(genofile, num.thread=2)

# Visualize similarity as a Heatmap
image(ibs$ibs, main="Haplotype Similarity Heatmap")

# Perform PCA
pca <- snpgdsPCA(genofile, num.thread=2)

# Create a data frame for plotting
pc_df <- data.frame(sample.id = pca$sample.id,
                    EV1 = pca$eigenvect[,1], # Principal Component 1
                    EV2 = pca$eigenvect[,2]) # Principal Component 2

# Plot the 300 samples
plot(pc_df$EV1, pc_df$EV2, xlab="PC1", ylab="PC2", main="Haplotype Clustering PCA")


# Cluster the samples based on the IBS distance
set.seed(123)
ibs_dist <- snpgdsHCluster(ibs)



# **********************************
# Use the standard R cutree function on the 'hclust' element
# Replace 3 with 10 if you want your 10 original groups
groups <- cutree(ibs_dist$hclust, k = 8)

# Create a data frame to see which Sample ID belongs to which Cluster
results <- data.frame(
  sample.id = ibs_dist$sample.id,
  cluster = groups
)

# View the first few rows
head(results)
# **********************************************++



# Merge the cluster groups into your PCA data frame
pc_df$cluster <- as.factor(results$cluster)

# Plot with colors based on the clusters
# 'pch = 19' makes solid circles, 'col' assigns colors by group
plot(pc_df$EV1, pc_df$EV2, 
     col = pc_df$cluster, 
     pch = 19, 
     xlab = "PC1", ylab = "PC2", 
     main = "PCA Colored by Haplotype Clusters")

# Add a legend so you know which color is which group
legend("topright", 
       legend = levels(pc_df$cluster), 
       col = 1:nlevels(pc_df$cluster), 
       pch = 19, 
       title = "Cluster ID")




# ************************
# Install and load dendextend if you haven't already
# .libPaths(c("/courses/software/R/lib/4.4", .libPaths()))
# install.packages("dendextend")
library(dendextend)

# Convert your hclust object to a dendrogram
dend <- as.dendrogram(ibs_dist$hclust)

# Color the branches based on your 8 groups (k = 8)
dend <- color_branches(dend, k = 8)

# Plot the colored tree
# 'leaflab = "none"' hides sample IDs if the plot is too crowded
plot(dend, main = "Haplotype Clustering Dendrogram", leaflab = "none")






# Check if directory exists, if not, create it
if(!dir.exists("~/vcf_treamt")) dir.create("~/vcf_treamt")
setwd("~/vcf_treamt")








# *************************************
# In R: Export Sample IDs for each group
for(i in 1:8) {
  write.table(results$sample.id[results$cluster == i], 
              file = paste0("group_", i, "_samples.txt"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}


# Define your 8 groups (1=Healthy, 2-8=Effects)
# Ensure 'results' is matched with your GDS sample order
phenotype <- results$cluster 

# Open your GDS file (from previous steps)
genofile <- snpgdsOpen("~/vcf_treamt/multisample_300.gds")


# Example: Finding variants key to Group 8 (Heart Problems)



# 1. Compare Group 8 (Heart) vs Group 1 (Healthy)
case_ids <- results$sample.id[results$cluster == 8]
control_ids <- results$sample.id[results$cluster == 1]

# 2. Calculate Allele Frequencies using the correct function: snpgdsSNPRateFreq
# This returns a data frame with a column called 'AlleleFreq'
af_case_data <- snpgdsSNPRateFreq(genofile, sample.id = case_ids)
af_control_data <- snpgdsSNPRateFreq(genofile, sample.id = control_ids)

# 3. Extract the frequencies
af_case <- af_case_data$AlleleFreq
af_control <- af_control_data$AlleleFreq

# 4. Calculate the Difference (Delta AF)
# A high difference indicates a mutation potentially linked to the group's effect
delta_af <- abs(af_case - af_control)

# 5. Create the Key Mutation Table
snps_info <- snpgdsSNPList(genofile)
association_results <- data.frame(
  snp.id = snps_info$snp.id,
  chr = snps_info$chromosome,
  pos = snps_info$position,
  af_healthy = af_control,
  af_heart = af_case,
  diff = delta_af
)

# 6. View the Top 10 "Key" Mutations for Heart Problems
key_mutations_heart <- association_results[order(-association_results$diff), ]
head(key_mutations_heart, 10)



# ******************************
# 1. Calculate Allele Frequency for all 8 groups

# Close all open handles first to avoid errors
showfile.gds(closeall=TRUE)
genofile <- snpgdsOpen("~/vcf_treamt/multisample_300.gds")

# Calculate AF for all 8 groups
af_results <- list()
for(i in 1:8) {
  sample_ids <- results$sample.id[results$cluster == i]
  af_results[[paste0("G", i)]] <- snpgdsSNPRateFreq(genofile, sample.id = sample_ids)$AlleleFreq
}

# Create master matrix
af_matrix <- as.data.frame(do.call(cbind, af_results))
snps <- snpgdsSNPList(genofile)
af_matrix <- cbind(snps[,c("chromosome", "position")], af_matrix)


# Efficacy markers: High in Group 8, Low in Group 1
efficacy_markers <- af_matrix[af_matrix$G8 > 0.40 & af_matrix$G1 < 0.05, ]
write.csv(efficacy_markers, "Efficacy_Markers_8vs1.csv")


# Example: Find markers for Group 7 (Bleeding)
# Logic: High in G7, but Low in G1 (Healthy) AND Low in G8 (Success)
toxicity_bleeding <- af_matrix[af_matrix$G7 > 0.30 & 
                                 af_matrix$G1 < 0.05 & 
                                 af_matrix$G8 < 0.05, ]

# Generalizing for all "Negative" groups (2-7)
# This identifies mutations distinguishing negative groups from the 'Control+Success' baseline
af_matrix$baseline_max <- apply(af_matrix[, c("G1", "G8")], 1, max)

# Extract variants unique to each side-effect group
for(i in 2:7) {
  group_name <- paste0("G", i)
  unique_side_effect <- af_matrix[af_matrix[[group_name]] > 0.30 & af_matrix$baseline_max < 0.05, ]
  write.csv(unique_side_effect, paste0("Unique_Markers_Group_", i, ".csv"))
}

##  **************
# 1. Define the Good Baseline (Average of G1 and G8)
af_matrix$good_baseline <- (af_matrix$G1 + af_matrix$G8) / 2

# 2. Identify Distinguishing Markers for Negative Groups
# We look for a "Difference" > 0.40 compared to the Good Baseline
# Example: Group 7 (Bleeding)
key_bleeding <- af_matrix[abs(af_matrix$G7 - af_matrix$good_baseline) > 0.40, ]

# 3. Analyze all Negative Groups (G2-G7) simultaneously
# This finds mutations where a group differs from the healthy/success standard
for(i in 2:7) {
  group_col <- paste0("G", i)
  
  # Find variants where the specific group deviates significantly from G1/G8
  distinguishing_vars <- af_matrix[abs(af_matrix[[group_col]] - af_matrix$good_baseline) > 0.50, ]
  
  # Export the results
  write.csv(distinguishing_vars, paste0("Distinguishing_Markers_", group_col, ".csv"), row.names=FALSE)
}