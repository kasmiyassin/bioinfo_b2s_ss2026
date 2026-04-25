## Why scRNA
The scRNA method can be used to:

* Explore which cell types are present in a tissue
* Identify unknown/rare cell types or states
* Elucidate the changes in gene expression during differentiation processes or across time or states
* Identify genes that are differentially expressed in particular cell types between conditions (e.g. treatment or disease)
* Explore changes in expression among a cell type while incorporating spatial, regulatory, and/or protein information


## 1. install R 

* Setup R environment for single-cell analysis.
* Establish good data and metadata management.
* Demonstrate how to import data and set up project for upcoming quality control analysis.


```bash
# update indices
sudo apt update -qq

sudo apt update
sudo apt install build-essential -y

# install two helper packages we need
sudo apt install --no-install-recommends software-properties-common dirmngr
# add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc 
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the repo from CRAN -- lsb_release adjusts to 'noble' or 'jammy' or ... as needed
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
# install R itself
sudo apt install --no-install-recommends r-base

# change lib to course
echo 'export R_LIBS_SITE="/courses/software/R_libs"' | sudo tee /etc/profile.d/course_r_libs.sh

```

## 2. install Rstudio

```bash
sudo apt-get install gdebi-core

wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2026.01.2-418-amd64.deb

sudo gdebi rstudio-server-2026.01.2-418-amd64.deb
```
## 3. install packages



install new package do

```r
# install one by one not all together
#(1) Install the 10 packages listed below from CRAN using the install.packages() function.
install.packages("ggplot2",dependencies = TRUE, lib="/courses/software/R_libs")

install.packages("tidyverse",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("Matrix",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("RCurl",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("scales",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("cowplot",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("BiocManager",dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("Seurat",dependencies = TRUE, lib="/courses/software/R_libs")


BiocManager::install("multtest", lib = "/courses/software/R_libs")
install.packages("metap", dependencies = TRUE, lib = "/courses/software/R_libs")
install.packages("reshape2", dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("plyr", dependencies = TRUE, lib="/courses/software/R_libs")
install.packages("devtools", dependencies = TRUE, lib="/courses/software/R_libs")

# (2) Install the 4 packages listed below from Bioconductor using the the BiocManager::install() function.

library(BiocManager)
BiocManager::install("AnnotationHub",dependencies = TRUE, lib="/courses/software/R_libs")
BiocManager::install("ensembldb",dependencies = TRUE, lib="/courses/software/R_libs")
BiocManager::install("multtest",dependencies = TRUE, lib="/courses/software/R_libs")
BiocManager::install("glmGamPoi",dependencies = TRUE, lib="/courses/software/R_libs")

# (3) Install Presto from GitHub using the devtools::install_github() function:
library(devtools)
devtools::install_github("immunogenomics/presto")

# (4) Finally, please check that all the packages were installed successfully by loading them one at a time using the library() function.

library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(cowplot)
library(BiocManager)
library(metap)
library(reshape2)
library(plyr)
library(devtools)
library(AnnotationHub)
library(ensembldb)
library(multtest)
library(glmGamPoi)
library(presto)

```





---
---

# Loading single cell RNA-seq data into Seurat

## Learning Objectives:

- Setup R environment for single-cell analysis.
- Establish good data and metadata management.
- Demonstrate how to import data and set up project for upcoming quality control analysis.




## Quality Control of Cell Ranger Output

* Describe how cellranger is run and what the ouputs are
* Review the cellranger generated QC report (web summary HTML)
* Create plots with cellranger metrics

### Metrics evaluation

```R
## Cellranger metrics evaluation

# Load libraries
library(tidyverse)
library(reshape2)
library(plyr)

# Names of samples (same name as folders stored in data)
samples <- c("ctrl", "stim")

# Loop over each sample and read the metrics summary in
metrics <- list()
for (sample in samples) {
    path_csv <- paste0("data/", sample, "_metrics_summary.csv")
    df <- read.csv(path_csv)
    rownames(df) <- sample
    df$sample <- sample
    metrics[[sample]] <- df
}
# Concatenate each sample metrics together
metrics <- ldply(metrics, rbind)

# Remove periods and percentages to make the values numeric
metrics <- metrics %>%
  column_to_rownames(".id") %>%
  mutate(across(everything(), parse_number(str_replace(., ",", "")))) %>%
  mutate(across(everything(), parse_number(str_replace(., "%", ""))))

colnames(metrics)


```

With all of this information available as a dataframe, we can use ggplot to visualize these values. As an example of how this information can be used, we can display what percentage of reads map to the various parts of the genome (Intergentic, Intronic, and Exonic).



```R

# Columns of interest
cols <- c("Reads.Mapped.Confidently.to.Intergenic.Regions",
          "Reads.Mapped.Confidently.to.Intronic.Regions",
          "Reads.Mapped.Confidently.to.Exonic.Regions",
          "sample")

# Data wrangling to sculpt dataframe in a ggplot friendly manner
df <- metrics %>%
    select(all_of(cols)) %>%
    melt(id.vars="sample") %>%
    mutate(variable = str_replace_all(variable, "Reads.Mapped.Confidently.to.", "")) %>%
    mutate(variable = str_replace_all(variable, ".Regions", "")) %>%
    mutate(value = str_replace_all(value, "%", "")) %>% # Remove % from value
    mutate(value = as.numeric(value)) # Make percentage numeric

# ggplot code to make a barplot
df %>% ggplot() +
    geom_bar(
        aes(x = sample, y = value, fill = variable),
        position = "stack",
        stat = "identity") +
    coord_flip() +
    labs(
        x = "Sample",
        y = "Percentage of Reads",
        title = "Percent of Reads Mapped to Each Region",
        fill = "Region")


```


The most important files that are generated during this cellranger run are the two matrix folders, which contain the count matrices from the experiment:

* raw_feature_bc_matrix
* filtered_feature_bc_matrix

In the next lesson, we will use raw_feature_bc_matrix to load the counts into Seurat. You can similarly do the same with filtered_feature_bc_matrix, the difference being that the filtered matrix has removed cells that cellranger determined as low quality using a variety of different tools. We chose to start with the raw counts matrix in this lesson so that you can better see what metrics are used to determine which cells are considered high quality.

## Theory of PCA

* Derive the covariance matrix used for Principal Components Analysis
* Explain the roles of eigenvectors and eigenvalues within a Principal Components Analysis
* Compare and contrast our Principal Components Analysis to the output of prcomp()

### Setting up to calculate Eigenvalues and Eigenvectors

```R
library(tidyverse)

## Creating an example data set

# Create a vector for Cell IDs
cells <- c("Cell_1", "Cell_2", "Cell_3", "Cell_4")
# Create a vector to hold expression values for Gene A across all of the cells
Gene_A <- c(0, 12, 65, 23)
# Create a vector to hold expression values for Gene B across all of the cells
Gene_B <- c(4, 30, 57, 18)

# Create a tibble to hold the cell names and expression values
expression_tibble <- tibble(cells, Gene_A, Gene_B)

# View the expression tibble
expression_tibble


# Let’s get an idea of what our data looks like by plotting it:

# Create a plot to view the raw expression data
ggplot(expression_tibble, aes(x = Gene_A, y = Gene_B, label = cells)) +
  geom_point(color = "cornflowerblue") +
  geom_text(hjust = 0, vjust = -1) +
  xlim(0, 80) +
  ylim(0, 80) +
  theme_bw() +
  xlab("Gene A") +
  ylab("Gene B") +
  ggtitle("Example Expression Values from Four Cells") +
  theme(plot.title = element_text(hjust = 0.5))

## Re-centering the dataset

# Determine the center of the data by:
# Finding the average expression of gene A
Gene_A_mean <- mean(expression_tibble$Gene_A)
# Finding the average expression of gene B
Gene_B_mean <- mean(expression_tibble$Gene_B)

# Create a vector to hold the center of the data
center_of_data <- c(Gene_A_mean, Gene_B_mean)
# Assign names to the components of the vector
names(center_of_data) <- c("Gene_A", "Gene_B")
# Print out the center_of_data vector
center_of_data


# View where the center of the data is located
ggplot(expression_tibble, aes(x = Gene_A, y = Gene_B, label = cells)) +
  geom_point(color = "cornflowerblue") +
    annotate("point", x = Gene_A_mean, y = Gene_B_mean, color = "red", size = 3) +
  geom_text(hjust = 0, vjust = -1) +
  annotate("text", x = Gene_A_mean, y = Gene_B_mean, color = "red", label="New Center", hjust = 0, vjust = -1) +
  xlim(0, 80) +
  ylim(0, 80) +
  theme_bw() +
  xlab("Gene A") +
  ylab("Gene B") +
  ggtitle("Example Expression Values from Four Cells") +
  theme(plot.title = element_text(hjust = 0.5))


## 2. Translate the points from their raw expression coordinates to coordinates centered around the center of the data

# Shift the data points so that they data is centered on the origin
recentered_expression_tibble <- expression_tibble %>%
  mutate(
    Gene_A = Gene_A - Gene_A_mean,
    Gene_B = Gene_B - Gene_B_mean
  )

# Print out the re-centered data
recentered_expression_tibble


# Plot the raw data after it has been shifted to have the center of the data align with the origin
ggplot(recentered_expression_tibble, aes(x = Gene_A, y = Gene_B, label = cells)) +
  geom_point(color = "cornflowerblue") +
  annotate("point", x = 0, y = 0, color = "red", size = 3) +
  geom_text(hjust = 0, vjust = -1) +
  annotate("text", x = 0, y = 0, color = "red", label="New Center", hjust = 0, vjust = -1) +
  xlim(-50, 50) +
  ylim(-50, 50) +
  theme_bw() +
  xlab("Gene A") +
  ylab("Gene B") +
  ggtitle("Example Re-centered Expression Values from Four Cells") +
  theme(plot.title = element_text(hjust = 0.5))


```



## Loading single cell RNA-seq data into Seurat


I start here: 

1. move the directory of work

```R 
# cd ~/session03scRNA
setwd("~/session03scRNA")

```
2. call R or rstudio in server

```bash 
R 
```

3. in R call the necessary libraries

```R
# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
```

4. create the qc file test in rstudio

5. check file in bash by writing

```bash 
cd ~/session03scRNA

# copy data
cp /courses/master_b2s/raw_data/s03_scrna/single_cell_rnaseq/ ./*
tree

# read file
zcat data/ctrl_raw_feature_bc_matrix/barcodes.tsv.gz  | head

```


to do it in R or Rstudio

```R
# mkdir ~/session03scRNA
# Note: showWarnings = FALSE prevents an error if the folder already exists
dir.create("~/session03scRNA", showWarnings = FALSE)

# cd ~/session03scRNA
setwd("~/session03scRNA")

# cp -r /courses/master_b2s/raw_data/s03_scrna/single_cell_rnaseq/ ./*

source_path <- "/courses/master_b2s/raw_data/s03_scrna/single_cell_rnaseq/"
files_to_copy <- list.files(source_path, full.names = TRUE)

file.copy(from = files_to_copy, to = ".", recursive = TRUE)

# tree
list.files(recursive = TRUE)

# Or for a pretty tree-like view:
# install.packages("fs")
# library(fs)
fs::dir_tree()

# zcat data/ctrl_raw_feature_bc_matrix/barcodes.tsv.gz | head

# Use gzfile to handle the compression and readLines to limit the output
con <- gzfile("data/ctrl_raw_feature_bc_matrix/barcodes.tsv.gz")
readLines(con, n = 10)
close(con)
```





## Reading in a single sample

After processing 10X data using its proprietary software Cell Ranger, you will have an outs directory (always). Within this directory you will find a number of different files including the files listed below:



```R
# How to read in 10X data for a single sample (output is a sparse matrix)
ctrl_counts <- Read10X(data.dir = "data/ctrl_raw_feature_bc_matrix")

# Turn count matrix into a Seurat object (output is a Seurat object)
ctrl <- CreateSeuratObject(counts = ctrl_counts,
                           min.features = 100)
ctrl

# Explore the metadata
# Seurat automatically creates some metadata for each of the cells when you use the Read10X() function to read in data. This information is stored in the meta.data slot within the Seurat object.


head(ctrl@meta.data)

```

## Reading in multiple samples with a for loop

Today we will use it to iterate over the two sample folders and execute two commands for each sample as we did above for a single sample -

Read in the count data (Read10X()) and
Create the Seurat objects from the read in data (CreateSeuratObject())
Go ahead and copy and paste the code below into your script and then run it.



```R

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


# Now that we have created both of these objects, let’s take a quick look at the list we just created. We should see that there are two Seurat objects in our list that correspond to each sample.

list_seurat


# Next, we need to merge these objects together into a single Seurat object. This will make it easier to run the QC steps for both sample groups together and enable us to easily compare the data quality for all the samples.

# Create a merged Seurat object
merged_seurat <- merge(x = list_seurat[["ctrl"]], 
                       y = list_seurat[["stim"]], 
                       add.cell.id = c("ctrl", "stim"))

merged_seurat

#This indicates that the matrices for each sample are being stored separately. Instead, we want to concatenate these together so that in the future when we run normalization, it is run on the entire dataset collectively. So to create one counts matrix, that is not sample/batch specific, we run the JoinLayers() function.

# Concatenate the count matrices of both samples together
merged_seurat <- JoinLayers(merged_seurat)
merged_seurat


# Check that the merged object has the appropriate sample-specific prefixes
# First several rows of a the merged_seurat@meta.data
head(merged_seurat@meta.data)

# Check that the merged object has the appropriate sample-specific prefixes
# Last several rows of a the merged_seurat@meta.data
tail(merged_seurat@meta.data)


```

## Quality Control Analysis

* Construct quality control metrics and visually evaluate the quality of the data
* Apply appropriate filters to remove low quality cells

### Generating quality metrics

```R 
# Explore merged metadata
View(merged_seurat@meta.data)

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Additional metadata columns
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

## Saving the updated metadata to our Seurat object

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata
                           
# Create .RData object to load at any time
save(merged_seurat, file="data/merged_filtered_seurat.RData")


## Assessing the quality metrics

# Visualize the number of cell counts per sample
metadata %>% 
    ggplot(aes(x=sample, fill=sample)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells")

## UMI counts (transcripts) per cell
# Visualize the number UMIs/transcripts per cell
metadata %>% 
    ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)

## Genes detected per cell
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
    ggplot(aes(color=sample, x=nGene, fill= sample)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 300)

## Complexity
# Visualize the overall complexity of the gene expression by visualizing 
# the genes detected per UMI (novelty score)
metadata %>%
    ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)


## Mitochondrial counts ratio
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
    ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.2)


## Joint filtering effects

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



## Filtering
### Cell-level filtering

# Filter out low quality cells using selected 
# thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                         subset = (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (mitoRatio < 0.20))
filtered_seurat



## Gene-level filtering
# Extract counts
counts <- GetAssayData(object = filtered_seurat, layer = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0


# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]


# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, 
                                      meta.data = filtered_seurat@meta.data)
filtered_seurat


## Re-assess QC metrics
# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data


## Saving filtered cells
# Create .RData object to load at any time
save(filtered_seurat, file="data/seurat_filtered.RData")

```


## Normalization and exploring data for unwanted variation

* Discuss why normalizing counts is necessary for accurate comparison between cells
* Describe different normalization approaches
* Evaluate the effects from any unwanted sources of variation and correct for them


### Explore sources of unwanted variation

```{r}
#| label: load_libraries
#| echo: false
# Load libraries
library(Seurat)
library(tidyverse)
library(knitr)
load("data/seurat_filtered.RData")
```

Let's start by creating a new script for the normalization and integration steps. Create a new script (File -> New File -> R script), and save it as `SCT_integration_analysis.R`.

For the remainder of the workflow we will be mainly using functions available in the Seurat package. Therefore, we need to load the Seurat library in addition to the tidyverse library and a few others listed below. 

```{r}
#| label: set-up_script
# Single-cell RNA-seq - normalization

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
```

Before we make any comparisons across cells, we will **apply a simple normalization.** This is solely for the purpose of exploring the sources of variation in our data.

The input for this analysis is a `seurat` object. We will use the one that we created in the QC lesson called `filtered_seurat`.

```{r}
#| label: normalize_data
# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)
```

Next, we take this normalized data and check to see if data correction methods are necessary. 


### Evaluating effects of cell cycle 

To assign each cell a score based on its expression of G2/M and S phase markers, we can use the Seuart function `CellCycleScoring()`. This function calculates cell cycle phase scores based on canonical markers that required as input.

We have provided a list of human cell cycle markers for you in the `data` folder as an Rdata file called `cycle.rda`. However, if you are not working with human data we have [additional materials](Aside_cell_cycle_scoring.qmd) detailing how to acquire cell cycle markers for other organisms of interest.

```{r}
#| label: cell_cycle
# Load cell cycle markers
load("data/cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)
```

Let's inspect our updated metadata:

```{r}
#| label: inspect_cell_cycle_scores
#| eval: false
# View cell cycle scores and phases assigned to cells
View(seurat_phase@meta.data)     
```

```{r}
#| label: tbl-cell_cycle_scores
#| tbl-cap: "Inspecting that the cell cycle scores and phases have been appropriately added to our metadata."
#| echo: false
# View cell cycle scores and phases assigned to cells
seurat_phase@meta.data %>% 
  head(n = 5) %>% 
  remove_rownames() %>%
  relocate(cells) %>%
  DT::datatable() %>% 
  DT::formatStyle("cells", 
                  "white-space" = "nowrap")
```

After scoring the cells for cell cycle, we would like to **determine whether cell cycle is a major source of variation in our dataset using PCA**. 

### PCA

Principal Component Analysis (PCA) is a technique used to emphasize variation as well as similarity, and to bring out strong patterns in a dataset; it is one of the methods used for *"dimensionality reduction"*. 

### Using PCA to evaluate the effects of cell cycle

To perform PCA, we need to **first choose the most variable features, then scale the data**. Since highly expressed genes exhibit the highest amount of variation and we don't want our 'highly variable genes' only to reflect high expression, we need to scale the data to scale variation with expression level. The Seurat `ScaleData()` function will scale the data by:

- Adjusting the expression of each gene to give a mean expression across cells to be 0
- Scaling expression of each gene to give a variance across cells to be 1

```R
#| label: find_variable_genes
# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
```


For the `selection.method` and `nfeatures` arguments the values specified are the default settings. Therefore, you do not necessarily need to include these in your code. We have included it here for transparency and inform you what you are using.
:::

Highly variable gene selection is extremely important since many downstream steps are computed only on these genes. Seurat allows us to access the ranked highly variable genes with the `VariableFeatures()` function. We can additionally visualize the dispersion of all genes using Seurat's `VariableFeaturePlot()`, which shows a gene's average expression across all cells on the x-axis and variance on the y-axis. Ideally we want to use genes that have high variance since this can indicate a change in expression depending on populations of cells. Adding labels using the `LabelPoints()` helps us understand which genes will be driving shape of our data.

```{r}
#| label: fig-top_15_variable_genes
#| fig-cap: Average expression and variance for each gene, labelling the top most variable genes.
# Identify the 15 most highly variable genes
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:15]

# Plot the average expression and variance of these genes
# With labels to indicate which genes are in the top 15
p <- VariableFeaturePlot(seurat_phase)
LabelPoints(plot = p, points = top_genes, repel = TRUE)
```

Now, we can perform the PCA analysis and plot the first two principal components against each other. We also split the figure by cell cycle phase, to evaluate similarities and/or differences. **We do not see large differences due to cell cycle phase. Based on this plot, we would not regress out the variation due to cell cycle.**

```{r}
#| label: fig-run_pca
#| fig-cap: Cells plotted on PC1 vs PC2 and colored by cell cycle phase.
# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
```

## Normalization and regressing out sources of unwanted variation using SCTransform
## Iterating over samples in a dataset

Since we have two samples in our dataset (from two conditions), we want to keep them as separate objects and transform them as that is what is required for integration. We will first split the cells in `seurat_phase` object into "Control" and "Stimulated":

```{r}
#| label: split_seurat_object
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")
split_seurat
```


Now we will **use a 'for loop'** to run the `SCTransform()` on each sample, and regress out mitochondrial expression by specifying in the `vars.to.regress` argument of the `SCTransform()` function.

Before we run this `for loop`, we know that the output can generate large R objects/variables in terms of memory. If we have a large dataset, then we might need to **adjust the limit for allowable object sizes within R** (*Default is 500 x 1024 ^2 = 500 Mb*) using the following code:

```{r}
#| label: increase_memory
options(future.globals.maxSize = 4000 * 1024^2)
```

Now, we run the following loop to **perform the SCTransform on all samples**. This may take some time (~10 minutes):

```{r}
#| label: SCTransform_all_samples
for (i in 1:length(split_seurat)) {
    split_seurat[[i]] <- SCTransform(split_seurat[[i]], 
                                     vars.to.regress = c("mitoRatio"),
                                     vst.flavor = "v2")
    }
```

Please note that in the for loop above, we specify that `vst.flavor = "v2"` to use the updated version of SCT. "v2" was introduced in early 2022, and is now commonly used. This update improves:

- Speed and memory consumption
- The stability of parameter estimates
- Variable feature identification in subsequent steps

For more information, please see the [Seurat vignette's section on SCTransform, v2 regularization](https://satijalab.org/seurat/articles/sctransform_v2_vignette.html). 

Note, the last line of output specifies **"Set default assay to SCT"**. This specifies that moving forward we would like to use the data after SCT was implemented. We can view the different assays that we have stored in our seurat object.

```{r}
#| label: check_assays
# Check which assays are stored in objects
split_seurat$ctrl@assays
```

Now we can see that in addition to the raw RNA counts, we now have a SCT component in our `assays` slot. The most variable features will be the only genes stored inside the SCT assay. As we move through the scRNA-seq analysis, we will choose the most appropriate assay to use for the different steps in the analysis. 



***

::: callout-tip
# [**Exercise 2**](07_SCT_normalization-Answer_key.qmd#exercise-2)

1. Are the same assays available for the "stim" samples within the `split_seurat` object? What is the code you used to check that?

2. Any observations for the genes or features listed under *"First 10 features:"* and the *"Top 10 variable features:"* for "ctrl" versus "stim"?
:::

***


### Save the object

Before finishing up, let's save this object to the `data/` folder. It can take a while to get back to this stage especially when working with large datasets, it is best practice to save the object as an easily loadable file locally.

```{r}
#| label: save_seurat_object
# Save the split seurat object
saveRDS(split_seurat, "data/split_seurat.rds")
```

::: callout-note
To load the `.rds` file back into your environment you would use the following code:
```{r}
#| label: load_seurat_object_example
#| eval: false
## DO NOT RUN
# Load the split seurat object into the environment
split_seurat <- readRDS("data/split_seurat.rds")
```
:::



---
---



```{r}
#| label: load_libraries
#| echo: false
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
load("data/seurat_filtered.RData")
seurat_phase <- NormalizeData(filtered_seurat)
load("data/cycle.rda")
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
seurat_phase <- ScaleData(seurat_phase)
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:15]
p <- VariableFeaturePlot(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)
```

# Exercise 1

Mitochondrial expression is another factor which can greatly influence clustering. Oftentimes, it is useful to regress out variation due to mitochondrial expression. However, if the differences in mitochondrial gene expression represent a biological phenomenon that may help to distinguish cell clusters, then we advise not regressing this out. In this exercise, we can perform a quick check similar to looking at cell cycle and decide whether or not we want to regress it out.

First, turn the mitochondrial ratio variable into a new categorical variable based on quartiles (using the code below):

```{r}
#| label: mito_ratio_wrangling
# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                   breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                   labels=c("Low","Medium","Medium high", "High"))
```

1. Next, plot the PCA similar to how we did with cell cycle regression. *Hint: use the new `mitoFr`variable to split cells and color them accordingly.*

```{r}
#| label: plot_pca
# Plot the PCA colored by mitoFr
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")
```

2. Evaluate the PCA plot generated above
	
	* Determine whether or not you observe an effect.
	
	Yes, there is an effect. 
	
	* Describe what you see. 
	
	Based on this plot, we can see that there is a different pattern of scatter for the plot containing cells with "High" mitochondrial expression. We observe that the lobe of cells on the left-hand side of the plot is where most of the cells with high mitochondrial expression are. For all other levels of mitochondrial expression we see a more even distribution of cells across the PCA plot.
	
	* Would you regress out mitochondrial fraction as a source of unwanted variation?
	
	Since we see this clear difference, we will regress out the 'mitoRatio' when we identify the most variant genes.

# Exercise 2

```{r}
#| label: split_object
split_seurat <- SplitObject(seurat_phase, split.by = "sample")
options(future.globals.maxSize = 4000 * 1024^2)
for (i in 1:length(split_seurat)) {
    split_seurat[[i]] <- SCTransform(split_seurat[[i]], 
                                     vars.to.regress = c("mitoRatio"),
                                     vst.flavor = "v2")
    }
```

1. Are the same assays available for the "stim" samples within the `split_seurat` object? What is the code you used to check that?

Yes they are available. The code use is:

```{r}
#| label: check_assays
split_seurat$stim@assays
```

2. Any observations for the genes or features listed under "First 10 features:" and the "Top 10 variable features:" for "ctrl" versus "stim"?

For the first 10 features, it appears that the same genes are present in both "ctrl" and "stim"

For the top 10 variable features, these are different in the the 2 conditions with some overlap between them.

---
---

# Integration

This part introduces participants to the concepts involved in integrating single-cell RNA-seq datasets using canonical correlation analysis (CCA) within the Seurat framework. Participants will learn when integration is appropriate, how it aligns shared cell types across conditions and why evaluating data before integration is essential for accurate downstream analyses.


### Learning Objectives

* Describe the theory of integration with CCA

## To integrate or not to integrate?

Generally, we always look at our clustering **without integration** before deciding whether we need to perform any alignment. It can be helpful to first run through clustering with samples from different sample classes together to see whether there are condition-specific clusters for cell types present in both conditions. Oftentimes, when clustering cells from multiple conditions there are condition-specific clusters and integration can help ensure the same cell types cluster together.

 **Do not just always perform integration because you think there might be differences - explore the data.** If we had performed the normalization on both conditions together in a Seurat object and visualized the similarity between cells, we would have seen in our dataset there is condition-specific clustering.
 
We will discuss the UMAP algorithm in more detail in a future section, but for now we can see that once we calculate our UMAP coordinates, there is a clear split based upon `sample`.

```{r}
#| label: fig-UMAP_phases_plot
#| fig-cap: Sample-specific clustering in UMAP space, before integration.
# Run UMAP
seurat_phase <- RunUMAP(seurat_phase,
                        dims = 1:40,
                        reduction = "pca")
# Plot UMAP
DimPlot(seurat_phase,
        group.by = "sample")   
```


Condition-specific clustering of the cells indicates that we need to integrate the cells across conditions to ensure that cells of the same cell type cluster together. 
_**If cells cluster by sample, condition, batch, dataset, modality, performing integration can help align cells across the groups to greatly improve the clustering and the downstream analyses**._

**Why is it important that cells of the same cell type cluster together?** 

We want to identify  _**cell types which are present in all samples/conditions/modalities**_ within our dataset, and therefore would like to observe a representation of cells from both samples/conditions/modalities in every cluster. This will enable more interpretable results downstream (i.e. DE analysis, ligand-receptor analysis, differential abundance analysis...).

So then how would you determine if integration is necessary? In this dataset, we know that the gene `LYZ` is a marker for monocytes. Even if the monocyte populations differ slightly between experimental conditions, we still expect monoctyes from all batches to be biologically similar. Because of this, these cells should be considered comparable and should occupy the same region in the UMAP embedding.

```{r}
#| label: fig-lyz_featureplot
#| fig-cap: Expression of gene LYZ on the unintegrated UMAP to show when integration is necessary.
FeaturePlot(seurat_phase,
            features = "LYZ")
```

However we can clearly see that there is a split in the grouping of monocytes that is driven by batch, as seen in the previous UMAP plot.

In this lesson, we will cover the integration of our samples across conditions, which is adapted from the [Seurat Guided Integration Tutorial](https://satijalab.org/seurat/articles/integration_introduction.html).
