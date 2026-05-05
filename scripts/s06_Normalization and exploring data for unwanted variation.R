# Normalization and exploring data for unwanted variation ####
## Methods for scRNA-seq normalization ####
## Explore sources of unwanted variation ####


# Add your custom directory to the list of places R looks for gems
.libPaths(c("/courses/software/R_libs", .libPaths()))

# Now you can load them normally

# Load libraries
library(Seurat)
library(tidyverse)
library(knitr)
#load("data/seurat_filtered.RData")


### Set-up ####
# Single-cell RNA-seq - normalization
# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)


# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)

### Evaluating effects of cell cycle  ####

# Load cell cycle markers
load("data/cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells
View(seurat_phase@meta.data)    


# View cell cycle scores and phases assigned to cells
seurat_phase@meta.data %>% 
  head(n = 5) %>% 
  remove_rownames() %>%
  relocate(cells) %>%
  DT::datatable() %>% 
  DT::formatStyle("cells", 
                  "white-space" = "nowrap")



### Using PCA to evaluate the effects of cell cycle ####


# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Identify the 15 most highly variable genes
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:15]

# Plot the average expression and variance of these genes
# With labels to indicate which genes are in the top 15
p <- VariableFeaturePlot(seurat_phase)
LabelPoints(plot = p, points = top_genes, repel = TRUE)


# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase")

DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

# to do regression https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette.html
# marrow <- ScaleData(seurat_phase, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat_phase))
# Now, a PCA on the variable genes no longer returns components associated with cell cycle
#  marrow <- RunPCA(seurat_phase, features = VariableFeatures(seurat_phase), nfeatures.print = 10)



## Normalization and regressing out sources of unwanted variation using SCTransform ####
### Regressing out covariates ####
### Iterating over samples in a dataset ####

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")
split_seurat


# increase_memory
options(future.globals.maxSize = 4000 * 1024^2)


#  label: SCTransform_all_samples
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], 
                                   vars.to.regress = c("mitoRatio"),
                                   vst.flavor = "v2")
}


# Check which assays are stored in objects
split_seurat$ctrl@assays


# Save the split seurat object
saveRDS(split_seurat, "data/split_seurat_1.rds")





