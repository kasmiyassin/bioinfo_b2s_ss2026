library(SNPRelate)
library(dendextend) # For colored trees

setwd("~/vcf_treamt/treatment01/")
# --- STEP A: Data Loading & Conversion ---
vcf_gz <- "/courses/master_b2s/raw_data/s02_vc/treatment/simulated.vcf.gz"
gds_fn <- "~/vcf_treamt/treatment01/multisample_12G.gds"

# Close any previous open files to prevent errors
showfile.gds(closeall=TRUE)

# snpgdsVCF2GDS handles .vcf.gz files directly
snpgdsVCF2GDS(vcf_gz, gds_fn, method="biallelic.only")
genofile <- snpgdsOpen(gds_fn)

# --- STEP B: Haplotype Similarity & Clustering ---
ibs <- snpgdsIBS(genofile, num.thread=2)
# Update IBS calculation
# 
ibs <- snpgdsIBS(genofile, num.thread = 2, autosome.only = FALSE)
ibs_dist <- snpgdsHCluster(ibs)

# PCA
pca <- snpgdsPCA(genofile, num.thread=2)
pca <- snpgdsPCA(genofile, num.thread = 2, autosome.only = FALSE)

# *******************
# method Elbow
# determine number of groups
# 1. Use the PCA eigenvectors for clustering (faster than the whole VCF)
# We use the first 10 PCs as they contain most of the genetic variance
matrix_for_clustering <- pca$eigenvect[, 1:10]

# 2. Calculate WSS for k = 1 to 15
wss <- sapply(1:40, function(k){
  kmeans(matrix_for_clustering, centers = k, nstart = 20)$tot.withinss
})

# 3. Plot the Elbow
plot(1:40, wss, type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares",
     main="Elbow Method for Optimal K")

# method 2 The Silhouette
library(cluster)

# Calculate average silhouette width for k = 2 to 12
avg_sil <- sapply(2:40, function(k){
  km <- kmeans(matrix_for_clustering, centers = k, nstart = 20)
  ss <- silhouette(km$cluster, dist(matrix_for_clustering))
  mean(ss[, 3])
})

# Plot the results
plot(2:40, avg_sil, type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Average Silhouette Width",
     main="Silhouette Method for Optimal K")
# *******************************

# Cut into 8 groups based on your clinical observations
groups <- cutree(ibs_dist$hclust, k = 10)
results <- data.frame(sample.id = ibs_dist$sample.id, cluster = groups)

# --- STEP C: PCA Visualization ---

pc_df <- data.frame(sample.id = pca$sample.id, EV1 = pca$eigenvect[,1], 
                    EV2 = pca$eigenvect[,2], cluster = as.factor(groups))

plot(pc_df$EV1, pc_df$EV2, col = pc_df$cluster, pch = 19, 
     main = "WGS Haplotype Clusters (300 Samples)", xlab = "PC1", ylab = "PC2")
legend("topright", legend = levels(pc_df$cluster), col = 1:11, pch = 19, title = "Group")



# --- STEP D: Allele Frequency Matrix ---
af_list <- list()
for(i in 1:11) {
  ids <- results$sample.id[results$cluster == i]
  af_list[[paste0("G", i)]] <- snpgdsSNPRateFreq(genofile, sample.id = ids)$AlleleFreq
}

af_matrix <- as.data.frame(do.call(cbind, af_list))
snps <- snpgdsSNPList(genofile)
af_matrix <- cbind(snps[,c("chromosome", "position")], af_matrix)

# --- STEP E: Distinguishing Negative Effects (G2-G7) from (G1 & G8) ---
# Create a baseline of "Good Outcomes" (Healthy + Successful)
af_matrix$good_outcome <- (af_matrix$G1 + af_matrix$G8) / 2

# Example: Finding markers unique to Group 5 (Diarrhea)
# We look for variants where G5 is DIFFERENT from G1 and G8
key_diarrhea <- af_matrix[abs(af_matrix$G5 - af_matrix$good_outcome) > 0.50, ]

# Export for clinical report
write.csv(key_diarrhea, "Key_Markers_Diarrhea_G5.csv", row.names=FALSE)
