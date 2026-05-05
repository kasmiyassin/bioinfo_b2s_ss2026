# Clustering Quality Control ####
## Importing data to Rstudio ####

# Add your custom directory to the list of places R looks for gems
.libPaths(c("/courses/software/R_libs", .libPaths()))

# Now you can load them normally
# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(cowplot)

load(bzfile("data/seurat_integrated.RData.bz2"))


## Exploration of quality control metrics #####
### Segregation of clusters by sample ####

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample)

# Barplot of number of cells per cluster by sample
ggplot(n_cells, aes(x=ident, y=n, fill=sample)) +
  geom_bar(position=position_dodge(), stat="identity") +
  theme_classic() +
  geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))



# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()

# Barplot of proportion of cells in each cluster by sample
ggplot(seurat_integrated@meta.data) +
  geom_bar(aes(x=integrated_snn_res.0.8, fill=sample), 
           position=position_fill())  +
  theme_classic()


### Segregation of clusters by cell cycle phase ####
# Explore whether clusters segregate by cell cycle phase
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()


### Segregation of clusters by various sources of uninteresting variation #####

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)


# Boxplot of nGene per cluster
ggplot(seurat_integrated@meta.data) +
  geom_boxplot(aes(x=integrated_snn_res.0.8, y=nGene, 
                   fill=integrated_snn_res.0.8)) +
  theme_classic() +
  NoLegend()

### Exploration of the PCs driving the different clusters ####

# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:18),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)
head(pc_data)


# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:18),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)

# ********************************************************
# Extract the UMAP coordinates for the first 10 cells
seurat_integrated@reductions$umap@cell.embeddings[1:10, 1:2]

# Defining the information in the Seurat object of interest
columns <- c(paste0("PC_", 1:18),
             "ident",
             "UMAP_1", "UMAP_2")
columns

# Defining the information in the seurat object of interest
# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:18),
             "ident",
             "umap_1", "umap_2")  # <-- THIS IS THE FIX (LOWERCASE)

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)

# Look at the data to confirm it worked
head(pc_data)
# *******************************************************

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  dplyr::summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(umap_1, umap_2))+
#         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)


# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)

## Exploring known cell type markers ####

# UMAP plot, representing each cell as a colored point corresponding with the cluster identified at resolution 0.8.
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE) + NoLegend()

# SCTransform dimensions**
dim(seurat_integrated[["RNA"]])
dim(seurat_integrated[["integrated"]])



# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
seurat_integrated




# **CD14+ monocyte markers**
  
# label: fig-CD4_FeaturePlot
# fig-cap: UMAP `FeaturePlot()` of top CD14+ monocyte markers.
# fig.width: 8
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD14", "LYZ"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


# **FCGR3A+ monocyte markers**

# label: fig-FCGR3A_monocyte_plot
# fig-cap: UMAP `FeaturePlot()` of top FCGR3A+ monocyte markers.
# fig.width: 8
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("FCGR3A", "MS4A7"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)




# **Macrophages**


# label: fig-Macrophages_plot
# fig-cap: UMAP `FeaturePlot()` of top FCGR3A+ Macrophages markers.
# fig.width: 8
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("MARCO", "ITGAM", "ADGRE1"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)



# label: fig-Macrophages_plot_vln
# fig-cap: Violin plot of top FCGR3A+ Macrophages markers.
# fig.height: 6
VlnPlot(seurat_integrated,
        c("MARCO", "ITGAM", "ADGRE1"),
        ncol = 1)

# **Conventional dendritic cell markers**

# label: fig-conventional_dendritic cell_plot
# fig-cap: UMAP `FeaturePlot()` of top Conventional dendritic cell markers.
# fig-width: 8
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("FCER1A", "CST3"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# **Plasmacytoid dendritic cell markers**
  
# label: fig-plasmacytoid_dendritic_cell_plot
# fig-cap: UMAP `FeaturePlot()` of top Plasmacytoid dendritic cell markers.
# fig.width: 8
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)



# List of known celltype markers
markers <- list()
markers[["CD14+ monocytes"]] <- c("CD14", "LYZ")
markers[["FCGR3A+ monocyte"]] <- c("FCGR3A", "MS4A7")
markers[["Macrophages"]] <- c("MARCO", "ITGAM", "ADGRE1")
markers[["Conventional dendritic"]] <- c("FCER1A", "CST3")
markers[["Plasmacytoid dendritic"]] <- c("IL3RA", "GZMB", "SERPINF1", "ITM2C")

# Create dotplot based on RNA expression
DotPlot(seurat_integrated, markers, assay="RNA")


# B cell
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD79A", "MS4A1"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# T cells
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD3D"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


# CD4+ T cells
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD3D", "IL7R", "CCR7"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# CD8+ T cells
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD3D", "CD8A"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# NK cells

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("GNLY", "NKG7"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


# Megakaryocytes

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("PPBP"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


# Erythrocytes
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("HBB", "HBA2"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


