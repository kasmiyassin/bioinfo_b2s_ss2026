# Clustering Analysis ####

## Clustering cells based on top PCs (metagenes) ####
### Set up ####

# Add your custom directory to the list of places R looks for gems
.libPaths(c("/courses/software/R_libs", .libPaths()))

# Now you can load them normally


# Single-cell RNA-seq - clustering

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)


## Identify significant PCs ####
# Heatmap of the first 9 principal components, showing expression across 500 cells for the top genes (negative and positive) for the component.
# Explore heatmap of PCs
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)

# try also 1000 and 1500
# not significant change

# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# Plot showing the standard deviation represented with each component to determine a cutoff value for the number of PCs.
# Plot the elbow plot
ElbowPlot(object = seurat_integrated, 
          ndims = 40)



## Cluster the cells ####
### Find neighbors ####

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

### Find clusters ####

# Determine the clusters for various resolutions                              
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.5, 0.6, 0.8, 1.0, 1.4))

## Visualize clusters of cells ####

# Explore resolutions
seurat_integrated@meta.data %>% 
  View()


# Explore resolutions
seurat_integrated@meta.data %>% 
  head(n = 5) %>% 
  remove_rownames() %>%
  relocate(cells) %>%
  DT::datatable() %>% 
  DT::formatStyle("cells", 
                  "white-space" = "nowrap")


# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"


# Calculate UMAP code
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)


# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

 
# if not identifical 
load(bzfile("data/seurat_integrated.RData.bz2"))

Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)



