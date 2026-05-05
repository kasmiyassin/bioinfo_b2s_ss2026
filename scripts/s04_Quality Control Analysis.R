# 4. Quality Control Analysis ####


# Add your custom directory to the list of places R looks for gems
.libPaths(c("/courses/software/R_libs", .libPaths()))

# Now you can load them normally

# Load merged seurat object
library(Seurat)
library(tidyverse)
library(R)
# merged_seurat <- readRDS("data/merged_seurat.RDS")

## Generating quality metrics ####
# Explore merged metadata
View(merged_seurat@meta.data)


## Novelty score ####

# create_novelty_score
# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)


## Mitochondrial Ratio ####
# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, 
                                                pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100


## Additional metadata columns ####

# Create metadata dataframe
metadata <- merged_seurat@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Create sample column
metadata$sample <- metadata$orig.ident


# Rename columns
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="data/merged_filtered_seurat_1.RData")


## Cell counts ####
# Barplot of the number of cells per sample.
# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


## UMI counts (transcripts) per cell ####

# Density plot showing the distribution of UMI counts per cell for each sample
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)


## Genes detected per cell ####

# Density plot showing the distribution of genes detected per cell for each sample
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)


## Complexity ####

# Density plot showing the overall complexity of gene expression per cell for each sample

# Visualize the overall complexity of the gene expression by visualizing 
# the genes detected per UMI (novelty score)
metadata %>%
  	ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8)


## Mitochondrial counts ratio ####

# Density plot showing the distribution of mitochondrial gene expression detected per cell for each sample.
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

## Reads per cell ####
### Joint filtering effects ####

# Scatterplot contrasting nUMI and nGenes while coloring each cell by mitochondrial ratio.
# Visualize the correlation between genes detected and number of UMIs and 
# determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

## Filtering ####
### Cell-level filtering ####

# Filter out low quality cells using selected 
# thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset = (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))
filtered_seurat


### Gene-level filtering ####

# extract_non_zero_counts
# Extract counts
counts <- GetAssayData(object = filtered_seurat, layer = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# retain_expressed_genes
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# reassign_filtered counts
# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, 
                                      meta.data = filtered_seurat@meta.data)
filtered_seurat



### Re-assess QC metrics ####

# exercise_extract_metadata
# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data
metadata_clean


### Saving filtered cells ####

# label: save_filtered_seurat_object
# Create .RData object to load at any time
save(filtered_seurat, file="data/seurat_filtered_1.RData")



