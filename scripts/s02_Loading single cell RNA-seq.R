# 2. Loading single cell RNA-seq data into Seurat ----
## Importing data to Rstudio ----

# Add your custom directory to the list of places R looks for gems
.libPaths(c("/courses/software/R/lib/4.4/", .libPaths()))
.libPaths()

# Now you can load them normally

# mkdir -p ~/session03scRNA
# Note: showWarnings = FALSE prevents an error if the folder already exists
dir.create("~/session03scRNA_test", showWarnings = FALSE)

# cd ~/session03scRNA
setwd("~/session03scRNA")

# cp -r /courses/master_b2s/raw_data/s03_scrna/single_cell_rnaseq/ ./*

source_path <- "/courses/master_b2s/raw_data/s03_scrna/"
files_to_copy <- list.files(source_path, full.names = TRUE)

file.copy(from = files_to_copy, to = ".", recursive = TRUE)

# tree
list.files(recursive = TRUE)

# Or for a pretty tree-like view:
# install.packages("fs")
library(fs)

fs::dir_tree()


# zcat data/ctrl_raw_feature_bc_matrix/barcodes.tsv.gz | head
# Use gzfile to handle the compression and readLines to limit the output
con <- gzfile("data/ctrl_raw_feature_bc_matrix/barcodes.tsv.gz")
readLines(con, n = 100)
close(con)



## Loading libraries ----

# Add your custom directory to the list of places R looks for gems
.libPaths(c("/courses/software/R_libs", .libPaths()))

# Now you can load them normally

# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

## Loading single-cell RNA-seq count data ----

### Single Sample or Group ----

# How to read in 10X data for a single sample (output is a sparse matrix)
ctrl_counts <- Read10X(data.dir = "data/ctrl_raw_feature_bc_matrix")

# Turn count matrix into a Seurat object (output is a Seurat object)
ctrl <- CreateSeuratObject(counts = ctrl_counts,
                           min.features = 100)
ctrl


# Explore the metadata
head(ctrl@meta.data)


# *****************************************************
### Reading in multiple samples with a `for loop` ----
# *****************************************************

# Step 1: Specify inputs

sample_names <- c("ctrl", "stim")

# Empty list to populate seurat object for each sample
list_seurat <- list()

for (sample in sample_names) {
  # Path to data directory
  data_dir <- paste0("data/", sample, "_raw_feature_bc_matrix")
  
  # Create a Seurat object for each sample
  # Step 2: Read in data for the input
  
  seurat_data <- Read10X(data.dir = data_dir)
  
  # Step 3: Create Seurat object from the 10X count data
  
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.features = 100,
                                   project = sample)
  
  # Save seurat object to list
  # Step 4: Assign Seurat object to a new variable based on sample
  
  list_seurat[[sample]] <- seurat_obj
}

# list_seurat_objects
list_seurat


# Create a merged Seurat object
merged_seurat <- merge(x = list_seurat[["ctrl"]], 
                       y = list_seurat[["stim"]], 
                       add.cell.id = c("ctrl", "stim"))

merged_seurat


wir <- print("We have problem >> 2 layers present: counts.ctrl, counts.stim!!!! Stop stop")

JoinLayers(wir)

# Concatenate the count matrices of both samples together
merged_seurat <- JoinLayers(merged_seurat)
merged_seurat

# Check that the merged object has the appropriate sample-specific prefixes
# first lines
head(merged_seurat@meta.data)

# Check that last lines: the merged object has the appropriate sample-specific prefixes
tail(merged_seurat@meta.data)


# Last several rows of a the `merged_seurat@meta.data`
merged_seurat@meta.data %>% tail() %>% R::kable()


