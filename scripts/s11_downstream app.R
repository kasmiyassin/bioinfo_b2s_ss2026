# Marker identification ####

# Add your custom directory to the list of places R looks for gems
.libPaths(c("/courses/software/R_libs", .libPaths()))

# Now you can load them normally
# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(knitr)
library(DT)






# Add celltype annotation as a column in meta.data 
seurat_subset_labeled$celltype <- Idents(seurat_subset_labeled)

# Compute number of cells per celltype
n_cells <- FetchData(seurat_subset_labeled, 
                     vars = c("celltype", "sample")) %>%
  dplyr::count(celltype, sample)

# Barplot of number of cells per celltype by sample
ggplot(n_cells, aes(x=celltype, y=n, fill=sample)) +
  geom_bar(position=position_dodge(), stat="identity") +
  theme_classic() +
  geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))


# Subset seurat object to just B cells
seurat_b_cells <- subset(seurat_subset_labeled, subset = (celltype == "B cells"))

# Run a wilcox test to compare ctrl vs stim
Idents(seurat_b_cells) <- "sample"
b_markers <- FindMarkers(seurat_b_cells,
                         ident.1 = "ctrl",
                         ident.2 = "stim",
                         grouping.var = "sample",
                         only.pos = FALSE,
                         logfc.threshold = 0.25)


library(EnhancedVolcano)
EnhancedVolcano(b_markers,
                row.names(b_markers),
                x="avg_log2FC",
                y="p_val_adj",
                title="B Cells",
                subtitle="Stim vs. Ctrl"
)



# ***************************************

# Subset seurat object to just B cells
seurat_m_cells <- subset(seurat_subset_labeled, subset = (celltype == "CD14+ monocytes"))

# Run a wilcox test to compare ctrl vs stim
Idents(seurat_m_cells) <- "sample"
m_markers <- FindMarkers(seurat_m_cells,
                         ident.1 = "ctrl",
                         ident.2 = "stim",
                         grouping.var = "sample",
                         only.pos = FALSE,
                         logfc.threshold = 0.25)


library(EnhancedVolcano)
EnhancedVolcano(m_markers,
                row.names(m_markers),
                x="avg_log2FC",
                y="p_val_adj",
                title="Monocyzes Cells",
                subtitle="Stim vs. Ctrl"
)

