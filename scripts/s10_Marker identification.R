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

load(bzfile("data/seurat_integrated.RData.bz2"))



## Identification of conserved markers in all conditions ####

DefaultAssay(seurat_integrated) <- "RNA"
Idents(seurat_integrated) <- "integrated_snn_res.0.8"



# label: FindConservedMarkers_cluster_0
# eval: false
cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)

# Inspect the output of FindConservedMarkers
View(cluster0_conserved_markers)



# label: tbl-_FindConservedMarkers_cluster_0
# tbl-cap: "Output from FindConservedMarkers"
cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)

cluster0_conserved_markers %>% 
  head(n = 10) %>% 
  kable()


### Adding Gene Annotations ####

# label: read_in_annotations
annotations <- read.csv("data/annotation.csv")

# Combine markers with gene descriptions 
cluster0_ann_markers <- cluster0_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

head(cluster0_ann_markers)



# Combine markers with gene descriptions 
cluster0_ann_markers <- cluster0_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

cluster0_ann_markers %>%
  head() %>%
  knitr::kable()


### Running on multiple samples ####

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}


# Iterate function across desired clusters
conserved_markers <- map_dfr(c(4,0,6,2), get_conserved)
head(conserved_markers)


# Iterate function across desired clusters
conserved_markers <- map_dfr(c(4,0,6,2), get_conserved)
conserved_markers %>%
  head() %>%
  DT::datatable() %>% 
  DT::formatStyle("description", 
                  "white-space" = "nowrap")



# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (ctrl_avg_log2FC + stim_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)


# Visualize top 10 markers per cluster
View(top10)


# Visualize top 10 markers per cluster
top10 %>%
  head() %>%
  DT::datatable() %>% 
  DT::formatStyle("description", 
                  "white-space" = "nowrap")


### Visualizing marker genes ####


# Plot interesting marker gene expression for cluster 4
FeaturePlot(object = seurat_integrated, 
            features = c("HSPH1", "HSPE1", "DNAJB1"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)


# Vln plot - cluster 4
VlnPlot(object = seurat_integrated, 
        features = c("HSPH1", "HSPE1", "DNAJB1"))


## Identifying gene markers for each cluster ####

# Determine differentiating markers for CD4+ T cell
cd4_tcells <- FindMarkers(seurat_integrated,
                          ident.1 = 2,
                          ident.2 = c(0,4,6))                  

# Add gene symbols to the DE table
cd4_tcells <- cd4_tcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# Reorder columns and sort by padj      
cd4_tcells <- cd4_tcells[, c(1, 3:5,2,6:7)]

cd4_tcells <- cd4_tcells %>%
  dplyr::arrange(p_val_adj) 


# View data
View(cd4_tcells)


# "Inspecting FindMarkers output"
cd4_tcells %>% 
  head() %>%
  DT::datatable() %>% 
  DT::formatStyle("description", 
                  "white-space" = "nowrap")


# Rename all identities
seurat_integrated <- RenameIdents(object = seurat_integrated, 
                                  "0" = "Naive or memory CD4+ T cells",
                                  "1" = "CD14+ monocytes",
                                  "2" = "Activated T cells",
                                  "3" = "CD14+ monocytes",
                                  "4" = "Stressed cells / Unknown",
                                  "5" = "CD8+ T cells",
                                  "6" = "Naive or memory CD4+ T cells",
                                  "7" = "B cells",
                                  "8" = "NK cells",
                                  "9" = "CD8+ T cells",
                                  "10" = "FCGR3A+ monocytes",
                                  "11" = "B cells",
                                  "12" = "NK cells",
                                  "13" = "B cells",
                                  "14" = "Conventional dendritic cells",
                                  "15" = "Megakaryocytes",
                                  "16" = "Plasmacytoid dendritic cells")


# Plot the UMAP
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)


# Remove the stressed or dying cells
seurat_subset_labeled <- subset(seurat_integrated,
                                idents = "Stressed cells / Unknown", invert = TRUE)

# Re-visualize the clusters
DimPlot(object = seurat_subset_labeled, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)

