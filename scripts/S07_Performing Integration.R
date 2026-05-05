# Performing Integration  ####

# Add your custom directory to the list of places R looks for gems
.libPaths(c("/courses/software/R_libs", .libPaths()))

# Now you can load them normally

## Running CCA ####
# Load the split seurat object into the environment
# split_seurat <- readRDS("data/split_seurat.rds")

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

# Find best buddies - can take a while to run
# note that the progress bar in your console will stay at 0%, but know that it is actually running.
# around 40 min
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

# 
saveRDS(integ_anchors, "data/integ_anchors_1.rds")
# integ_anchors <- readRDS("data/integ_anchors.rds")

# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
#
saveRDS(seurat_integrated, "data/seurat_integrated_1.rds")

# Rejoin the layers in the RNA assay that we split earlier
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])



### UMAP visualization ####
# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
PCAPlot(seurat_integrated,
        group.by = "sample",
        split.by = "sample")

# UMAP of dataset after integration with CCA.
# Set seed
set.seed(123456)

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP                             
DimPlot(seurat_integrated,
        group.by = "sample") 

# Plot UMAP split by sample
DimPlot(seurat_integrated,
        group.by = "sample",
        split.by = "sample")  

# Save integrated Seurat object
saveRDS(seurat_integrated, "data/integrated_seurat_2.rds")


# *****************************************************************
## unintergaTED data
# 1. Merge the split objects back into one
# This combines them but does NOT align them (no integration)
seurat_unintegrated <- merge(split_seurat[[1]], 
                             y = split_seurat[2:length(split_seurat)], 
                             add.cell.ids = names(split_seurat))

# 2. Set the variable features to the ones you picked for integration
VariableFeatures(seurat_unintegrated) <- integ_features

# 3. Run PCA and UMAP on the unintegrated data
seurat_unintegrated <- RunPCA(seurat_unintegrated, verbose = FALSE)
seurat_unintegrated <- RunUMAP(seurat_unintegrated, dims = 1:40)

# 4. Plot the results
DimPlot(seurat_unintegrated, group.by = "sample") + 
  ggtitle("Unintegrated Data (Batch Effect Visible)")


# ********************************************************

