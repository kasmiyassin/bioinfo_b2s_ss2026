# 7. Unbiased Marker Identification ####

# Add your custom directory to the list of places R looks for gems
.libPaths(c("/courses/software/R_libs", .libPaths()))

# Load necessary libraries (Make sure DESeq2 is installed for the pseudobulk step)
# BiocManager::install("DESeq2", lib="/courses/software/R_libs")
library(Seurat)
library(tidyverse)
library(cowplot)
library(DESeq2) 

BiocManager::install(version = "3.23",dependencies = TRUE,lib="/courses/software/R_libs")
BiocManager::install("DESeq2",dependencies = TRUE, lib="/courses/software/R_libs")
BiocManager::install("EnhancedVolcano",dependencies = TRUE, lib="/courses/software/R_libs")
BiocManager::install("SingleCellExperiment",dependencies = TRUE, lib="/courses/software/R_libs")
BiocManager::install("miloR",dependencies = TRUE, lib="/courses/software/R_libs")
BiocManager::install("clusterProfiler",dependencies = TRUE, lib="/courses/software/R_libs")
BiocManager::install("org.Mm.eg.db",dependencies = TRUE, lib="/courses/software/R_libs")
BiocManager::install("sccomp",dependencies = TRUE, lib="/courses/software/R_libs")
BiocManager::install("speckle",dependencies = TRUE, lib="/courses/software/R_libs")


install.packages("Seurat",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("tidyverse",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("pheatmap",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("RColorBrewer",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("cowplot",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("dplyr",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("ggalluvial",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("msigdbr",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("ggvenn",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("rlang", lib="/courses/software/R_libs", type="source",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages(c("rlang","cli","vctrs","pillar",dependencies = TRUE, lib="/courses/software/R_libs"),
                 lib="/courses/software/R_libs",
                 type="source")
# ... [Your previous code ends with FeaturePlot for Macrophages] ...

# 7. Unbiased Marker Identification ####

# Add your custom directory to the list of places R looks for gems
.libPaths(c("/courses/software/R_libs", .libPaths()))

# Load necessary libraries (Make sure DESeq2 is installed for the pseudobulk step)
# BiocManager::install("DESeq2", lib="/courses/software/R_libs")
library(Seurat, lib="/courses/software/R_libs")
library(tidyverse)
library(cowplot)
library(DESeq2) 

# Ensure we are using the RNA assay for marker identification (not integrated)
DefaultAssay(seurat_integrated) <- "RNA"

# IMPORTANT FOR SEURAT V5: Ensure layers are joined before running FindMarkers
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])

# Ensure your active identity is set to the clustering resolution you chose
# (Let's assume you settled on resolution 0.8 from your previous steps)
Idents(seurat_integrated) <- "integrated_snn_res.0.8"

## Find markers for every cluster compared to all remaining cells
# This uses the 'presto' package you installed for incredible speed
all_markers <- FindAllMarkers(object = seurat_integrated, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

# Extract top 5 markers per cluster to inspect
top5_markers <- all_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

print(top5_markers)

# Visualize top markers for a specific cluster (e.g., Cluster 0)
VlnPlot(seurat_integrated, features = top5_markers$gene[1:3], pt.size = 0)
FeaturePlot(seurat_integrated, features = top5_markers$gene[1:3])


## Find Conserved Markers across conditions (ctrl vs stim)
# This finds markers that define a cell type regardless of treatment
# Example for cluster 0:
cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE)
head(cluster0_conserved_markers)


# 8. Cell Type Annotation ####

# Once you've analyzed the markers (using literature, databases like CellMarker, or your top markers),
# you will rename your numbered clusters to actual cell types.

# NOTE: You will need to change these names based on YOUR actual marker results!
# This is a template based on standard PBMC naming.
seurat_annotated <- RenameIdents(object = seurat_integrated, 
                                 "0" = "CD14+ Monocytes",
                                 "1" = "CD4 T cells",
                                 "2" = "CD8 T cells",
                                 "3" = "B cells",
                                 "4" = "FCGR3A+ Monocytes",
                                 "5" = "NK cells",
                                 "6" = "Dendritic cells",
                                 "7" = "Megakaryocytes") # Add all your clusters

# Save the cluster names as a new metadata column
seurat_annotated$cell_type <- Idents(seurat_annotated)

# Plot UMAP with new cell type annotations
DimPlot(seurat_annotated, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 4,
        repel = TRUE) + 
  ggtitle("Annotated Cell Types")

# Save your annotated object!
saveRDS(seurat_annotated, "data/seurat_annotated.rds")


# 9. Pseudobulk Preparation (for Differential Expression) ####
# Following the HBC Pseudobulk tutorial logic

# We want to compare ctrl vs stim WITHIN a specific cell type. 
# First, we aggregate the counts to the sample level for each cell type.

# Create pseudobulk counts using Seurat's AggregateExpression
# This sums the counts of all cells of the same cell_type within the same sample
pseudo_data <- AggregateExpression(seurat_annotated, 
                                   assays = "RNA", 
                                   return.seurat = FALSE, 
                                   group.by = c("cell_type", "sample"))

# Extract the RNA count matrix
pb_counts <- pseudo_data$RNA

# View the dimensions and first few rows (columns are now celltype_sample)
dim(pb_counts)
head(pb_counts[, 1:4])

# 10. Differential Expression using DESeq2 ####

# Let's say we want to perform DE analysis on "CD14+ Monocytes" 
# to see what genes change between 'ctrl' and 'stim'

target_cell_type <- "CD14+ Monocytes"

# 1. Subset the pseudobulk count matrix to only include our target cell type
# We use regex to find columns that start with our target cell type
cols_to_keep <- grep(paste0("^", target_cell_type, "_"), colnames(pb_counts))
counts_subset <- pb_counts[, cols_to_keep]

# 2. Create sample metadata for DESeq2
# Extract the sample names from the column names
col_names <- colnames(counts_subset)
# Assuming format is "CellType_Sample", e.g., "CD14+ Monocytes_ctrl"
samples <- str_replace(col_names, paste0(target_cell_type, "_"), "")

metadata_pb <- data.frame(
  row.names = col_names,
  sample_id = samples,
  condition = factor(samples, levels = c("ctrl", "stim")) # Make sure 'ctrl' is first (reference)
)

# 3. Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                              colData = metadata_pb,
                              design = ~ condition)

# 4. Filter out genes with very low counts (optional but recommended)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 5. Run DESeq2
dds <- DESeq(dds)

# 6. Extract results
# Contrast: Condition stim vs ctrl
res <- results(dds, contrast = c("condition", "stim", "ctrl"))

# Convert to dataframe and sort by adjusted p-value
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "gene") %>%
  arrange(padj)

# View top differentially expressed genes
head(res_df)

# Save the results
write.csv(res_df, paste0("data/DE_results_", target_cell_type, "_stim_vs_ctrl.csv"), row.names = FALSE)

# 7. Volcano Plot of DE results (using EnhancedVolcano which you installed earlier)
library(EnhancedVolcano)

EnhancedVolcano(res_df,
                lab = res_df$gene,
                x = 'log2FoldChange',
                y = 'padj',
                title = paste0('Stim vs Ctrl in ', target_cell_type),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 6.0)