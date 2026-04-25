
# Single Cell RNA seq session

## Why scRNA
The scRNA method can be used to:

* Explore which cell types are present in a tissue
* Identify unknown/rare cell types or states
* Elucidate the changes in gene expression during differentiation processes or across time or states
* Identify genes that are differentially expressed in particular cell types between conditions (e.g. treatment or disease)
* Explore changes in expression among a cell type while incorporating spatial, regulatory, and/or protein information


## Preparation of server or computer

#### 1. install R 

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

### 2. install Rstudio

```bash
sudo apt-get install gdebi-core

wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2026.01.2-418-amd64.deb

sudo gdebi rstudio-server-2026.01.2-418-amd64.deb
```
### 3. install packages



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

<br>
<br>

# Loading single cell RNA-seq data into Seurat

## Learning Objectives:

- Setup R environment for single-cell analysis.
- Establish good data and metadata management.
- Demonstrate how to import data and set up project for upcoming quality control analysis.

## Importing data to Rstudio

### Bash working - Dont Do
```bash 
# 1. move the directory of work
cd ~/session03scRNA

# 2. call R or rstudio in server
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

### Rstudio working space - To Do
But we will do all wrok in R or Rstudio

```R
# mkdir -p ~/session03scRNA
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
library(fs)

fs::dir_tree()
```
we should see something like 
```R
#| label: data_structure
#| eval: false
single_cell_rnaseq/
├── data
├── figures
├── results
├── scripts
└── single_cell_rnaseq.Rproj
```

```R
# zcat data/ctrl_raw_feature_bc_matrix/barcodes.tsv.gz | head
# Use gzfile to handle the compression and readLines to limit the output
con <- gzfile("data/ctrl_raw_feature_bc_matrix/barcodes.tsv.gz")
readLines(con, n = 10)
close(con)

```


### Loading libraries 

Now, we can load the necessary libraries:

```R
#| label: load_libraries
#| warning: false
# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
```

### Loading single-cell RNA-seq count data 

Regardless of the technology or pipeline used to process your raw single-cell RNA-seq sequence data, the output with quantified expression will generally be the same. That is, for each individual sample you will have the following **three files**:

1. A file with the **cell IDs**, representing all cells quantified
2. A file with the **gene IDs**, representing all genes quantified
3. A **matrix of counts** per gene for every cell

We can explore these files by clicking the `data/ctrl_raw_feature_bc_matrix` folder:

1. `barcodes.tsv` 
2. `features.tsv`
3. `matrix.mtx`


#### Reading in a single sample

>> Dont run me

After processing 10X data using its proprietary software Cell Ranger, you will have an `outs` directory (always). Within this directory you will find a number of different files.

 
 If we had a single sample, we could generate the count matrix and then subsequently create a Seurat object:

The Seurat object is a custom list-like object that has well-defined spaces to store specific information/data. You can find more information about the slots in the Seurat object [at this link](https://github.com/satijalab/seurat/wiki/Seurat).

```R
#| label: read_in_single_dataset
# How to read in 10X data for a single sample (output is a sparse matrix)

ctrl_counts <- Read10X(data.dir = "data/ctrl_raw_feature_bc_matrix")

# Turn count matrix into a Seurat object (output is a Seurat object)
ctrl <- CreateSeuratObject(counts = ctrl_counts,
                           min.features = 100)
ctrl
```

**Seurat automatically creates some metadata** for each of the cells when you use the `Read10X()` function to read in data. This information is stored in the `meta.data` slot within the Seurat object. 

```R
#| label: head_ctrl_dataset
#| eval: false
# Explore the metadata
head(ctrl@meta.data)
```

```R
#| label: tbl-ctrl_metadata
#| tbl-cap: First several rows of a the `ctr@meta.data`
#| echo: false
ctrl@meta.data %>% head() %>% R::kable()
```

What do the columns of metadata mean?

- `orig.ident`: this often contains the sample identity if known, but will default to "SeuratProject"
- `nCount_RNA`: number of UMIs per cell
- `nFeature_RNA`: number of genes detected per cell

#### Reading in multiple samples with a `for loop`

In practice, you will likely have several samples that you will need to read in data for, and that can get tedious and error-prone if you do it one at a time. So, to make the data import into R more efficient we can use a `for` loop, which will iterate over a series of commands for each of the inputs given and create seurat objects for each of our samples. 


Today we will use it to **iterate over the two sample folders** and execute two commands for each sample as we did above for a single sample - 

1. Read in the count data (`Read10X()`) and
2. Create the Seurat objects from the read in data (`CreateSeuratObject()`)

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

```

Now that we have created both of these objects, let's take a quick look at the list we just created. We should see that there are two Seurat objects in our list that correspond to each sample.

```R
#| label: list_seurat_objects
list_seurat
```

Next, we need to merge these objects together into a single Seurat object. This will make it easier to run the QC steps for both sample groups together and enable us to easily compare the data quality for all the samples.  

We can use the `merge()` function from the Seurat package to do this. Here, we also specify `add.cell.id` because the same cell IDs can be used for different samples so we add a **sample-specific prefix** to the cell IDs to ensure that they are unique.


```R
#| label: merge_seurat_objects
# Create a merged Seurat object
merged_seurat <- merge(x = list_seurat[["ctrl"]], 
                       y = list_seurat[["stim"]], 
                       add.cell.id = c("ctrl", "stim"))

merged_seurat
```

However you may notice that when we look at our seurat object that we have 2 `layers` (count matrices) when it says:

`2 layers present: counts.ctrl, counts.stim`

This indicates that the matrices for each sample are being stored separately. Instead, we want to concatenate these together so that in the future when we run normalization, it is run on the entire dataset collectively. So to create one `counts` matrix, that is not sample/batch specific, we run the `JoinLayers()` function.

```R
#| label: JoinLayers
# Concatenate the count matrices of both samples together
merged_seurat <- JoinLayers(merged_seurat)
merged_seurat
```

Now that we’ve merged and joined the dataset as we wanted, we can verify that it looks correct. In the merged object’s metadata, we should see the rowname prefixes added by `add.cell.id` and the `orig.ident` values set during `CreateSeuratObject` with `project = sample`:

First several rows of a the merged_seurat@meta.data

```R
#| label: inspect_seurat_object_head
#| eval: false
# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
```

```R
#| label: tbl-merged_seurat_metadata_head
#| tbl-cap: First several rows of a the `merged_seurat@meta.data`
#| echo: false
merged_seurat@meta.data %>% head() %>% R::kable()
```

Last several rows of a the merged_seurat@meta.data

```R
#| label: inspect_seurat_object_tail
#| eval: false
# Check that the merged object has the appropriate sample-specific prefixes
tail(merged_seurat@meta.data)
```

```R
#| label: tbl-merged_seurat_metadata_tail
#| tbl-cap: Last several rows of a the `merged_seurat@meta.data`
#| echo: false
merged_seurat@meta.data %>% tail() %>% R::kable()
```


---
<br>
<br>

# Quality Control of Cell Ranger Output


## Learning Objectives:

- Describe how cellranger is run and what the ouputs are
- Review the cellranger generated QC report (web summary HTML)
- Create plots with cellranger metrics

> YK >>> show **cellranger**


## Metrics evaluation

Many of the core pieces of information from the web summary are stored in the `metrics_summary.csv`. As this is a csv file, we can read it into R and generate plots to include in reports on the general quality of the samples.

We have included these csv files for the control and stimulated samples in your **`data` directory**:

- Control sample `metrics.csv` file (`ctrl_metrics_summary.csv`)
- Stimulated sample `metrics.csv` file (`stim_metrics_summary.csv`)

First, to read the files in:

```R
#| label: load_metrics_data
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
```

The information available in this file include: 

```R
#| label: columns_in_metrics_data
colnames(metrics)
```

With all of this information available as a dataframe, we can use ggplot to visualize these values. As an example of how this information can be used, we can display what percentage of reads map to the various parts of the genome (Intergentic, Intronic, and Exonic).



```R
#| label: fig-mapping_reads_to_regions_plot
#| fig-cap: Mapping percentage of reads to different regions of the genome.
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


## Matrix folders

The most important files that are generated during this cellranger run are the two matrix folders, which contain the count matrices from the experiment:

- `raw_feature_bc_matrix`
- `filtered_feature_bc_matrix`

In the **next step, we will use `raw_feature_bc_matrix`** to load the counts into Seurat. You can similarly do the same with `filtered_feature_bc_matrix`, the difference being that the filtered matrix has removed cells that cellranger determined as low quality using a variety of different tools. We chose to start with the raw counts matrix in this lesson so that you can better see what metrics are used to determine which cells are considered high quality.

---

<br>
<br>


# Quality Control Analysis

## Learning Objectives:

- Construct quality control metrics and visually evaluate the quality of the data
- Apply appropriate filters to remove low quality cells


Let us load the data if not already there
```R
#| label: load_libraries
#| echo: false
# Load merged seurat object
library(Seurat)
library(tidyverse)
library(R)
merged_seurat <- readRDS("data/merged_seurat.RDS")
```


## Generating quality metrics

When data is loaded into Seurat and the initial object is created, there is some basic metadata asssembled for each of the cells in the count matrix. To take a close look at this metadata, let's view the data frame stored in the `meta.data` slot of our `merged_seurat` object:

```R
#| label: inspect_metadata
#| eval: false
# Explore merged metadata
View(merged_seurat@meta.data)
```

```R
#| label: tbl-table_metadata
#| tbl-cap: "Inspecting our metadata table to get a sense for the type of information stored within it."
#| echo: false
merged_seurat@meta.data %>% 
  head(n = 10) %>% 
  kable()
```

There are three columns of information:

- `orig.ident`: this column will contain the sample identity if known. It will default to the value we provided for the `project` argument when loading in the data
- `nCount_RNA`: this column represents  the number of UMIs per cell
- `nFeature_RNA`: this column represents the number of genes detected per cell

In order to create the appropriate plots for the quality control analysis, we need to calculate some additional metrics. These include:

- **Number of genes detected per UMI:** this metric with give us an idea of the complexity of our dataset (more genes detected per UMI, more complex our data)
- **Mitochondrial ratio:** this metric will give us a percentage of cell reads originating from the mitochondrial genes

## Novelty score

This value is quite easy to calculate, as we take the log10 of the number of genes detected per cell and the log10 of the number of UMIs per cell, then divide the log10 number of genes by the log10 number of UMIs. The novelty score and how it relates to complexity of the RNA species, is described in more detail later in this lesson.

```R
#| label: create_novelty_score
# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
```


## Mitochondrial Ratio

Seurat has a convenient function that allows us to calculate the **proportion of transcripts mapping to mitochondrial genes**. The `PercentageFeatureSet()` function takes in a `pattern` argument and searches through all gene identifiers in the dataset for that pattern. Since we are looking for mitochondrial genes, we are searching any gene identifiers that begin with the pattern "MT-". For each cell, the function takes the sum of counts across all genes (features) belonging to the "Mt-" set, and then divides by the count sum for all genes (features). This value is multiplied by 100 to obtain a percentage value. 

For our analysis, rather than using a percentage value we would prefer to work with the ratio value. As such, we will reverse that last step performed by the function by taking the output value and dividing by 100.


```R
#| label: compute_mito_ratio
# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, 
                                                pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
```

> **MT genes for organisms other than human**
> 
> The pattern provided ("^MT-") works for human gene names. You may need to adjust the pattern argument depending on your organism of interest. Additionally, if you weren't using gene names as the gene ID then this function wouldn't work as we have used it above as the pattern will not suffice. Since there are caveats to using this function, it is advisable to manually compute this metric. If you are interested, we have [code available to compute this metric on your own](https://github.com/hbctraining/scRNA-seq/blob/master/lessons/mitoRatio.md).
> 

### Additional metadata columns

We are a now all set with quality metrics required for assessing our data. However, we would like to **include some additional information** that would be useful to have in our metadata including **cell IDs and condition information**.

When we added columns of information to our metadata file above, we simply added it directly to the metadata slot in the Seurat object using the `$` operator. We could continue to do so for the next few columns of data, but instead we will extract the dataframe into a separate variable. In this way we can work with the metadata data frame as a seperate entity from the seurat object without the risk of affecting any other data stored inside the object.
 
Let's begin by creating the `metadata` dataframe by extracting the `meta.data` slot from the Seurat object: 

```R
#| label: create_metadata_dataframe
# Create metadata dataframe
metadata <- merged_seurat@meta.data
```

Next, we'll add a **new column for cell identifiers**. This information is currently located in the row names of our metadata dataframe. We will keep the rownames as is and duplicate it into a new column called `cells`:

```R
#| label: add_cell_IDs_to_metadata
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
```

Additionally we are going to create another column titled `sample` which is also stored in `orig.ident`. We do this because this because `sample` is a commonly used metadata field that others will look for when they download your data.

```R
#| label: add_sample_to_metadata
# Create sample column
metadata$sample <- metadata$orig.ident
```

And finally, we will **rename some of the existing columns** in our metadata dataframe to be more intuitive:



```R
#| label: rename_metadata_columns
# Rename columns
metadata <- metadata %>%
        dplyr::rename(nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
```

Now you are **all setup with the metrics you need to assess the quality of your data**! Your final metadata table will have rows that correspond to each cell, and columns with information about those cells.


### Saving the updated metadata to our Seurat object

Before we assess our metrics we are going to save all of the work we have done thus far back into our Seurat object. We can do this by simply assigning the dataframe into the `meta.data` slot:

```R
#| label: remerge_metadata
# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata
                           
# Create .RData object to load at any time
save(merged_seurat, file="data/merged_filtered_seurat.RData")
```

## Assessing the quality metrics

Now that we have generated the various metrics to assess, we can explore them with visualizations. We will assess various metrics and then decide on which cells are low quality and should be removed from the analysis:

- Cell counts
- UMI counts per cell
- Genes detected per cell
- Complexity (novelty score)
- Mitochondrial counts ratio


## What about doublets?

In single-cell RNA sequencing experiments, doublets are generated from two cells. They typically arise due to errors in cell sorting or capture, especially in droplet-based protocols involving thousands of cells. Doublets are obviously undesirable when the aim is to characterize populations at the single-cell level. In particular, they can incorrectly suggest the existence of intermediate populations or transitory states that do not actually exist. Thus, it is desirable to remove doublet libraries so that they do not compromise interpretation of the results.

Many workflows use maximum thresholds for UMIs or genes, with the idea that a much higher number of reads or genes detected indicate multiple cells. While this rationale seems to be intuitive, it is not accurate. Also, many of the tools used to detect doublets tend to get rid of cells with intermediate or continuous phenotypes, although they may work well on datasets with very discrete cell types. [Scrublet](https://github.com/AllonKleinLab/scrublet) is a popular tool for doublet detection. 

Currently, we recommend not including any thresholds at this point in time. When we have identified markers for each of the clusters, we suggest exploring the markers to determine whether the markers apply to more than one cell type.

### Cell counts

The cell counts are determined by the number of unique cellular barcodes detected. For this experiment, between 12,000 -13,000 cells are expected.

In an ideal world, you would expect the number of unique cellular barcodes to correpsond to the number of cells you loaded. However, this is not the case as capture rates of cells are only a proportion of what is loaded. For example, the inDrops cell **capture efficiency** is higher (70-80%) compared to 10X which is between 50-60%.

## Automated cell counting
The capture efficiency could appear much lower if the cell concentration used for library preparation was not accurate. Cell concentration should NOT be determined by FACS machine or Bioanalyzer (these tools are not accurate for concentration determination), instead use a hemocytometer or automated cell counter for calculation of cell concentration.


The cell numbers can also vary by protocol, **producing cell numbers that are much higher than what we loaded**. For example, during the inDrops protocol, the cellular barcodes are present in the hydrogels, which are encapsulated in the droplets with a single cell and lysis/reaction mixture. While each hydrogel should have a single cellular barcode associated with it, occasionally a hydrogel can have more than one cellular barcode. Similarly, with the 10X protocol there is a chance of obtaining only a barcoded bead in the emulsion droplet (GEM) and no actual cell.  Both of these, in addition to the presence of dying cells can lead to a higher number of cellular barcodes than cells.

```R
#| label: fig-cell_counts_per_sample
#| fig-cap: Barplot of the number of cells per sample.
# Visualize the number of cell counts per sample
metadata %>% 
  	ggplot(aes(x=sample, fill=sample)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells")
```


We see over 15,000 cells per sample, which is quite a bit more than the 12-13,000 expected. It is clear that we likely have some junk 'cells' present.

### UMI counts (transcripts) per cell

The UMI counts per cell should generally be above 500, that is the low end of what we expect. If UMI counts are between 500-1000 counts, it is usable but the cells probably should have been sequenced more deeply. 

```R
#| label: fig-UMIs_per_cell
#| fig-cap: Density plot showing the distribution of UMI counts per cell for each sample
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  	ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500)
```

We can see that majority of our cells in both samples have 1000 UMIs or greater, which is great. 

### Genes detected per cell

We have similar expectations for gene detection as for UMI detection, although it may be a bit lower than UMIs. For high quality data, the proportional histogram should contain **a single large peak that represents cells that were encapsulated**. If we see a **small shoulder** to the left of the major peak (not present in our data), or a bimodal distribution of the cells, that can indicate a couple of things. It might be that there are a set of **cells that failed** for some reason. It could also be that there are **biologically different types of cells** (i.e. quiescent cell populations, less complex cells of interest), and/or one type is much smaller than the other (i.e. cells with high counts may be cells that are larger in size). Therefore, this threshold should be assessed with other metrics that we describe in this lesson.
 
```R
#| label: fig-genes_detected
#| fig-cap: Density plot showing the distribution of genes detected per cell for each sample
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  	ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)
```

### Complexity

We can evaluate each cell in terms of how complex the RNA species are by using a measure called the novelty score. The novelty score is computed by taking the ratio of nGenes over nUMI. If there are many captured transcripts (high nUMI) and a low number of genes detected in a cell, this likely means that you only captured a low number of genes and simply sequenced transcripts from those lower number of genes over and over again. These low complexity (low novelty) cells could represent a specific cell type (i.e. red blood cells which lack a typical transcriptome), or could be due to an artifact or contamination. Generally, we expect the novelty score to be above 0.80 for good quality cells.

```R
#| label: fig-complexity_score_plot
#| fig-cap: Density plot showing the overall complexity of gene expression per cell for each sample

# Visualize the overall complexity of the gene expression by visualizing 
# the genes detected per UMI (novelty score)
metadata %>%
  	ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8)
```


### Mitochondrial counts ratio

This metric can identify whether there is a large amount of **mitochondrial contamination from dead or dying cells**. We define poor quality samples for mitochondrial counts as cells which surpass the 0.2 mitochondrial ratio mark, unless of course you are expecting this in your sample.

```R
#| label: fig-mito_ratio_plot
#| fig-cap: Density plot showing the distribution of mitochondrial gene expression detected per cell for each sample.

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  	ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)
```


## Reads per cell
This is another metric that can be useful to explore; however, the workflow used would need to save this information to assess. Generally, with this metric you hope to see all of the samples with peaks in relatively the same location between 10,000 and 100,000 reads per cell. 

### Joint filtering effects

Considering any of these QC metrics in isolation can lead to misinterpretation of cellular signals. For example, cells with a comparatively high fraction of mitochondrial counts may be involved in respiratory processes and may be cells that you would like to keep. Likewise, other metrics can have other biological interpretations.  A general rule of thumb when performing QC is to **set thresholds for individual metrics to be as permissive as possible, and always consider the joint effects** of these metrics. In this way, you reduce the risk of filtering out any viable cell populations. 


Two metrics that are often evaluated together are the number of UMIs and the number of genes detected per cell. Here, we have plotted the **number of genes versus the number of UMIs coloured by the fraction of mitochondrial reads**. Jointly visualizing the count and gene thresholds and additionally overlaying the mitochondrial fraction, gives a summarized persepective of the quality per cell.

```R
#| label: fig-joint_effects_plot
#| fig-cap: Scatterplot contrasting nUMI and nGenes while coloring each cell by mitochondrial ratio.
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
```


Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs (upper right quadrant of the plot). Cells that are **poor quality are likely to have low genes and UMIs per cell**, and correspond to the data points in the bottom left quadrant of the plot. With this plot we also evaluate the **slope of the line**, and any scatter of data points in the **bottom right hand quadrant** of the plot. These cells have a high number of UMIs but only a few number of genes. These could be dying cells, but also could represent a population of a low complexity celltype (i.e red blood cells).

**Mitochondrial read fractions are only high in particularly low count cells with few detected genes** (darker colored data points). This could be indicative of damaged/dying cells whose cytoplasmic mRNA has leaked out through a broken membrane, and thus, only mRNA located in the mitochondria is still conserved. We can see from the plot, that these cells are filtered out by our count and gene number thresholds. 


## Filtering

### Cell-level filtering

Now that we have visualized the various metrics, we can decide on the thresholds to apply which will result in the removal of low quality cells. Often the recommendations mentioned earlier are a rough guideline, and the specific experiment needs to inform the exact thresholds chosen. We will use the following thresholds:

- nUMI > 500
- nGene > 250
- log10GenesPerUMI > 0.8
- mitoRatio < 0.2

To filter, we wil go back to our Seurat object and use the `subset()` function:

```R
#| label: filtered_seurat_object
# Filter out low quality cells using selected 
# thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                         subset = (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (mitoRatio < 0.20))
filtered_seurat
```

### Gene-level filtering

Within our data we will have many genes with zero counts. These genes can dramatically reduce the average expression for a cell and so we will remove them from our data. We will start by identifying which genes have a zero count in each cell:

```R
#| label: extract_non_zero_counts
# Extract counts
counts <- GetAssayData(object = filtered_seurat, layer = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0
```

Now, we will perform some filtering by prevalence. If a gene is only expressed in a handful of cells, it is not particularly meaningful as it still brings down the averages for all other cells it is not expressed in. For our data we choose to **keep only genes which are expressed in 10 or more cells.** By using this filter, genes which have zero counts in all cells will effectively be removed.

```R
#| label: retain_expressed_genes
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
```

Finally, take those filtered counts and create a new Seurat object for downstream analysis.

```R
#| label: reassign_filtered counts
# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, 
                                      meta.data = filtered_seurat@meta.data)
filtered_seurat
```

### Re-assess QC metrics

After performing the filtering, it's recommended to look back over the metrics to make sure that your data matches your expectations and is good for downstream analysis. 

***

### [Exercise 1](05_quality_control-Answer_key.qmd#exercise-1)

Extract the new metadata from the filtered Seurat object using the code provided below:

```R
#| label: exercise_extract_metadata
# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data
```

1. Perform all of the same QC plots using the filtered data.

2. Report the number of cells left for each sample, and comment on whether the number of cells removed is high or low. Can you give reasons why this number is still not ~12K (which is how many cells were loaded for the experiment)?

3. After filtering for nGene per cell, you should still observe a small shoulder to the right of the main peak. What might this shoulder represent?

4. When plotting the nGene against nUMI do you observe any data points in the bottom right quadrant of the plot? What can you say about these cells that have been removed?

***

### Saving filtered cells

Based on these QC metrics we would identify any failed samples and move forward with our filtered cells. Often we iterate through the QC metrics using different filtering criteria; it is not necessarily a linear process. When satisfied with the filtering criteria, we would save our filtered cell object for clustering and marker identification.

```R
#| label: save_filtered_seurat_object
# Create .RData object to load at any time
save(filtered_seurat, file="data/seurat_filtered.RData")
```

> **Bad data**
> The data we are working with is pretty good quality. If you are interested in knowing what 'bad' data might look like when performing QC, we have some materials [linked here](Aside_QC_bad_data.qmd) where we explore similar QC metrics of a poor quality sample.




---

<br>
<br>

# Theory of PCA
> ## Theory of PCA
> 
> * Derive the covariance matrix used for Principal Components Analysis
> * Explain the roles of eigenvectors and eigenvalues within a Principal Components Analysis
> * Compare and contrast our Principal Components Analysis to the output of prcomp()
> 
> ### Setting up to calculate Eigenvalues and Eigenvectors
> 
> ```r
> library(tidyverse)
> 
> ## Creating an example data set
> 
> # Create a vector for Cell IDs
> cells <- c("Cell_1", "Cell_2", "Cell_3", "Cell_4")
> # Create a vector to hold expression values for Gene A across all of the cells
> Gene_A <- c(0, 12, 65, 23)
> # Create a vector to hold expression values for Gene B across all of the cells
> Gene_B <- c(4, 30, 57, 18)
> 
> # Create a tibble to hold the cell names and expression values
> expression_tibble <- tibble(cells, Gene_A, Gene_B)
> 
> # View the expression tibble
> expression_tibble
> 
> 
> # Let’s get an idea of what our data looks like by plotting it:
> 
> # Create a plot to view the raw expression data
> ggplot(expression_tibble, aes(x = Gene_A, y = Gene_B, label = cells)) +
>   geom_point(color = "cornflowerblue") +
>   geom_text(hjust = 0, vjust = -1) +
>   xlim(0, 80) +
>   ylim(0, 80) +
>   theme_bw() +
>   xlab("Gene A") +
>   ylab("Gene B") +
>   ggtitle("Example Expression Values from Four Cells") +
>   theme(plot.title = element_text(hjust = 0.5))
> 
> ## Re-centering the dataset
> 
> # Determine the center of the data by:
> # Finding the average expression of gene A
> Gene_A_mean <- mean(expression_tibble$Gene_A)
> # Finding the average expression of gene B
> Gene_B_mean <- mean(expression_tibble$Gene_B)
> 
> # Create a vector to hold the center of the data
> center_of_data <- c(Gene_A_mean, Gene_B_mean)
> # Assign names to the components of the vector
> names(center_of_data) <- c("Gene_A", "Gene_B")
> # Print out the center_of_data vector
> center_of_data
> 
> 
> # View where the center of the data is located
> ggplot(expression_tibble, aes(x = Gene_A, y = Gene_B, label = cells)) +
>   geom_point(color = "cornflowerblue") +
>     annotate("point", x = Gene_A_mean, y = Gene_B_mean, color = "red", size = 3) +
>   geom_text(hjust = 0, vjust = -1) +
>   annotate("text", x = Gene_A_mean, y = Gene_B_mean, color = "red", label="New Center", hjust = 0, vjust = -1) +
>   xlim(0, 80) +
>   ylim(0, 80) +
>   theme_bw() +
>   xlab("Gene A") +
>   ylab("Gene B") +
>   ggtitle("Example Expression Values from Four Cells") +
>   theme(plot.title = element_text(hjust = 0.5))
> 
> 
> ## 2. Translate the points from their raw expression coordinates to coordinates centered around the center of the data
> 
> # Shift the data points so that they data is centered on the origin
> recentered_expression_tibble <- expression_tibble %>%
>   mutate(
>     Gene_A = Gene_A - Gene_A_mean,
>     Gene_B = Gene_B - Gene_B_mean
>   )
> 
> # Print out the re-centered data
> recentered_expression_tibble
> 
> 
> # Plot the raw data after it has been shifted to have the center of the data align with the origin
> ggplot(recentered_expression_tibble, aes(x = Gene_A, y = Gene_B, label = cells)) +
>   geom_point(color = "cornflowerblue") +
>   annotate("point", x = 0, y = 0, color = "red", size = 3) +
>   geom_text(hjust = 0, vjust = -1) +
>   annotate("text", x = 0, y = 0, color = "red", label="New Center", hjust = 0, vjust = -1) +
>   xlim(-50, 50) +
>   ylim(-50, 50) +
>   theme_bw() +
>   xlab("Gene A") +
>   ylab("Gene B") +
>   ggtitle("Example Re-centered Expression Values from Four Cells") +
>   theme(plot.title = element_text(hjust = 0.5))
> ```
> 
> 





---

<br>
<br>

# Normalization and exploring data for unwanted variation

## Learning Objectives:

* Discuss why normalizing counts is necessary for accurate comparison between cells
* Describe different normalization approaches
* Evaluate the effects from any unwanted sources of variation and correct for them



## Methods for scRNA-seq normalization
ppt
## Explore sources of unwanted variation

Let load the data if not there


```R
#| label: load_libraries
#| echo: false
# Load libraries
library(Seurat)
library(tidyverse)
library(knitr)
load("data/seurat_filtered.RData")
```


The most common biological data correction (or source of "uninteresting" variation) in single cell RNA-seq is the effects of the cell cycle on the transcriptome. We need to explore the data and see if we observe any effects in our data.   

### Set-up

Let's start by creating a new script for the normalization and integration steps. Create a new script (File -> New File -> R script), and save it as `SCT_integration_analysis.R`.

For the remainder of the workflow we will be mainly using functions available in the Seurat package. Therefore, we need to load the Seurat library in addition to the tidyverse library and a few others listed below. 

```R
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

```R
#| label: normalize_data
# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)
```

Next, we take this normalized data and check to see if data correction methods are necessary. 


### Evaluating effects of cell cycle 

To assign each cell a score based on its expression of G2/M and S phase markers, we can use the Seuart function `CellCycleScoring()`. This function calculates cell cycle phase scores based on canonical markers that required as input.

We have provided a list of human cell cycle markers for you in the `data` folder as an Rdata file called `cycle.rda`. However, if you are not working with human data we have [additional materials](Aside_cell_cycle_scoring.qmd) detailing how to acquire cell cycle markers for other organisms of interest.

```R
#| label: cell_cycle
# Load cell cycle markers
load("data/cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)
```

Let's inspect our updated metadata:

```R
#| label: inspect_cell_cycle_scores
#| eval: false
# View cell cycle scores and phases assigned to cells
View(seurat_phase@meta.data)     
```

```R
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

PPT

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
 `FindVariableFeatures` arguments: 
For the `selection.method` and `nfeatures` arguments the values specified are the default settings. Therefore, you do not necessarily need to include these in your code. We have included it here for transparency and inform you what you are using.

Highly variable gene selection is extremely important since many downstream steps are computed only on these genes. Seurat allows us to access the ranked highly variable genes with the `VariableFeatures()` function. We can additionally visualize the dispersion of all genes using Seurat's `VariableFeaturePlot()`, which shows a gene's average expression across all cells on the x-axis and variance on the y-axis. Ideally we want to use genes that have high variance since this can indicate a change in expression depending on populations of cells. Adding labels using the `LabelPoints()` helps us understand which genes will be driving shape of our data.

```R
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

```R
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



```R

pbmc <- RunPCA(seurat_phase, features = VariableFeatures(object = seurat_phase))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") + NoLegend()

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(pbmc)

```





---

## Normalization and regressing out sources of unwanted variation using SCTransform

In the [Hafemeister and Satija, 2019 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1) the authors explored the issues with simple transformations. Specifically they evaluated the standard log normalization approach and found that genes with different abundances are affected differently and that **effective normalization (using the log transform) is only observed with low/medium abundance genes (Figure 1D, below)**. Additionally, **substantial imbalances in variance were observed with the log-normalized data (Figure 1E, below)**. In particular, cells with low total UMI counts exhibited disproportionately higher variance for high-abundance genes, dampening the variance contribution from other gene abundances.

The conclusion is, **we cannot treat all genes the same.**

The proposed solution was the use of **Pearson residuals for transformation**, as implemented in Seurat's `SCTransform` function. With this approach:

- Measurements are multiplied by a gene-specific weight
- Each gene is weighted based on how much evidence there is that it is non-uniformly expressed across cells
- More evidence == more of a weight; Genes that are expressed in only a small fraction of cells will be favored (useful for finding rare cell populations)
- Not just a consideration of the expression level is, but also the distribution of expression

**Why don't we just run SCTransform to normalize?**

While the functions `NormalizeData`, `VariableFeatures` and `ScaleData` can be replaced by the function `SCTransform`, the latter uses a more sophisticated way to perform the normalization and scaling. We suggest using log normalization because it is good to observe the data and any trends using a simple transformation, as methods like SCT can alter the data in a way that is not as intuitive to interpret. 


Now that we have established which effects are observed in our data, we can use the SCTransform method to regress out these effects. The **SCTransform** method was proposed as a better alternative to the log transform normalization method that we used for exploring sources of unwanted variation. The method not only **normalizes data, but it also performs a variance stabilization and allows for additional covariates to be regressed out**.

As described earlier, all genes cannot be treated the same. As such, the **SCTransform method constructs a generalized linear model (GLM) for each gene** with UMI counts as the response and sequencing depth as the explanatory variable. Information is pooled across genes with similar abundances, to regularize parameter estimates and **obtain residuals which represent effectively normalized data values** which are no longer correlated with sequencing depth.


### Regressing out covariates
Since the UMI counts are part of the GLM, the effects are automatically regressed out. The user can include any additional covariates (`vars.to.regress`) that may have an effect on expression and will be included in the model. 

To run the SCTransform we have the code below as an example. ****DO NOT RUN** this code**, as we prefer to run this for each sample separately in the next section below.


```R
#| label: SCTransform_example
#| eval: false
## **DO NOT RUN** CODE ##

# SCTransform
seurat_phase <- SCTransform(seurat_phase, 
                            vars.to.regress = c("mitoRatio"))
```

### Iterating over samples in a dataset

Since we have two samples in our dataset (from two conditions), we want to keep them as separate objects and transform them as that is what is required for integration. We will first split the cells in `seurat_phase` object into "Control" and "Stimulated":

```R
#| label: split_seurat_object
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")
split_seurat
```

**# Subsetting samples**:  If you only wanted to integrate on a subset of your samples (e.g. all Ctrl replicates only), you could select which ones you wanted from the `split_seurat` object as shown below and move forward with those.

```R
#| label: split_seurat_subset_example
#| eval: false
## **DO NOT RUN**
ctrl_reps <- split_seurat[c("ctrl_1", "ctrl_2")]
```


Now we will **use a 'for loop'** to run the `SCTransform()` on each sample, and regress out mitochondrial expression by specifying in the `vars.to.regress` argument of the `SCTransform()` function.

Before we run this `for loop`, we know that the output can generate large R objects/variables in terms of memory. If we have a large dataset, then we might need to **adjust the limit for allowable object sizes within R** (*Default is 500 x 1024 ^2 = 500 Mb*) using the following code:

```R
#| label: increase_memory
options(future.globals.maxSize = 4000 * 1024^2)
```

Now, we run the following loop to **perform the SCTransform on all samples**. This may take some time (~10 minutes):

```R
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

**SCTransform variable features:** By default, after normalizing, adjusting the variance, and regressing out uninteresting sources of variation, SCTransform will rank the genes by residual variance and output the 3000 most variant genes. If the dataset has larger cell numbers, then it may be beneficial to adjust this parameter higher using the `variable.features.n` argument.

Note, the last line of output specifies **"Set default assay to SCT"**. This specifies that moving forward we would like to use the data after SCT was implemented. We can view the different assays that we have stored in our seurat object.

```R
#| label: check_assays
# Check which assays are stored in objects
split_seurat$ctrl@assays
```

Now we can see that in addition to the raw RNA counts, we now have a SCT component in our `assays` slot. The most variable features will be the only genes stored inside the SCT assay. As we move through the scRNA-seq analysis, we will choose the most appropriate assay to use for the different steps in the analysis. 



### Save the object

Before finishing up, let's save this object to the `data/` folder. It can take a while to get back to this stage especially when working with large datasets, it is best practice to save the object as an easily loadable file locally.

```R
#| label: save_seurat_object
# Save the split seurat object
saveRDS(split_seurat, "data/split_seurat.rds")
```

To load the `.rds` file back into your environment you would use the following code:
```R
#| label: load_seurat_object_example
#| eval: false
## **DO NOT RUN**
# Load the split seurat object into the environment
split_seurat <- readRDS("data/split_seurat.rds")
```


---
<br>
<br>

# Performing Integration 

## Learning Objectives:

- Perform integration of cells across conditions to identify cells that are similar to each other
- Describe complex integration tasks and alternative tools for integration

## Introduction to integration
PPT

<br>

## Running CCA


First, we need to identify the shared variable genes for the integration. By default, this function only selects the top 2000 genes. In this step Seuart performs a more complex version of an intersect between the highly variable genes from each condition (based on SCTransform). We have specified 3000 genes for the size of the intersect set.
 
```R
#| label: integration_example
#| eval: false
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
```

***Loading `split_seurat`*** If you are missing the `split_seurat` object, you can first load it from your `data` folder:

```R
#| label: load_split_seurat
#| eval: false
# Load the split seurat object into the environment
split_seurat <- readRDS("data/split_seurat.rds")
```


If you do not have the `split_seurat.rds` file in your `data` folder, you can right-click [here](https://www.dropbox.com/scl/fi/7ion5yarjsko7rwfojzom/split_seurat.rds?rlkey=9x2b5t82y7hf805szneb6rnt2&st=knozf4r6&dl=1) to download it to the `data` folder (it may take a bit of time to download). 


Now, we need to **prepare the SCTransform object** for integration. This function basically prepares for integration analysis by ensuring all necessary data (specifically the SCTransform residuals) are present for the features chosen as anchors between datasets.

```R
#| label: prep_SCTIntegration
#| eval: false
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
```

Now, we are going to **perform CCA, find the best buddies or anchors and filter incorrect anchors**. For our dataset, this will take up to 15 minutes to run. *Also, note that the progress bar in your console will stay at 0%, but know that it is actually running.*

```R
#| label: integration_anchors
#| eval: false
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
```


Finally, we can **integrate across conditions**.

```R
#| label: integrate_data
#| eval: false
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# Rejoin the layers in the RNA assay that we split earlier
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])
```

### UMAP visualization

After integration, to visualize the integrated data we can use dimensionality reduction techniques, such as PCA and Uniform Manifold Approximation and Projection (UMAP). While PCA will determine all PCs, we can only plot two at a time. In contrast, UMAP will take the information from any number of top PCs to arrange the cells in this multidimensional space. It will take those distances in multidimensional space and plot them in two dimensions working to preserve local and global structure. In this way, the distances between cells represent similarity in expression. If you wish to explore UMAP in more detail, [this post](https://pair-code.github.io/understanding-umap/) is a nice introduction to UMAP theory.

To generate these visualizations we need to first run PCA and UMAP methods. Let's start with PCA.

```R
#| label: integrated_PCA
# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
PCAPlot(seurat_integrated,
        group.by = "sample",
        split.by = "sample")  
```


We can see with the PCA mapping that we have a good overlay of both conditions by PCA. 

Now, we can also **visualize with UMAP**. Let's run the method and plot. UMAP is a stochastic algorithm – this means that
it makes use of randomness both to speed up approximation steps, and to aid in solving hard optimization problems. Due to the stochastic nature, different runs of UMAP can produce different results. **We can set the seed to a specific (but random) number**, and this avoids the creation of a slightly different UMAP each time re-run our code.


**Reproducible randomness** Typically, the `set.seed()` would be placed at the beginning of your script. In this way the selected random number would be applied to any function that uses pseudorandom numbers in its algorithm.

```R
#| label: fig-integrated_UMAP
#| fig-cap: UMAP of dataset after integration with CCA.
# Set seed
set.seed(123456)

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
			                       reduction = "pca")

# Plot UMAP                             
DimPlot(seurat_integrated,
        group.by = "sample")                             
```

>> **Did integration improve our results?**

> When we compare the similarity between the ctrl and stim clusters in the above plot with what we see using the the unintegrated dataset, **it is clear that this dataset benefitted from the integration!**


### Side-by-side comparison of clusters

Sometimes it's easier to see whether all of the cells align well if we **split the plotting between conditions**, which we can do by adding the `split.by` argument to the `DimPlot()` function:

```R
#| label: fig-integrated_DimPlot
#| fig-cap: UMAP `split.by` sample to see if there are regions unique to one condition.
# Plot UMAP split by sample
DimPlot(seurat_integrated,
        group.by = "sample",
        split.by = "sample")  
```


### Save the "integrated" object!

Since it can take a while to integrate, it's often a good idea to **save the integrated seurat object**.

```R
#| label: save_integrated_object
#| eval: false
# Save integrated Seurat object
saveRDS(seurat_integrated, "data/integrated_seurat.rds")
```


---
<br>
<br>

# Clustering Analysis

  This part guides participants through the process of clustering analysis in single-cell RNA-seq data using Seurat. Participants learn how to evaluate principal components, construct K-nearest neighbor graphs, optimize clustering resolution and visualize clusters using dimensionality reduction techniques. The lesson emphasizes identifying meaningful biological clusters, recognizing technical artifacts and applying best practices for iterative refinement.

## Learning Objectives:

* Describe methods for evaluating the number of principal components used for clustering
* Perform clustering of cells based on significant principal components

<br>

Let us call our data if not already there

```R
#| label: load_libraries
#| echo: false
# Load libraries
library(Seurat)
load(bzfile("data/seurat_integrated.RData.bz2"))
```


## Clustering cells based on top PCs (metagenes)

### Set up

Before starting with this lesson, let's create a new script for the next few steps in the workflow called `clustering.R`. 

Next, let's load all the libraries that we need.

```R
#| label: script_set-up
# Single-cell RNA-seq - clustering

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
```


## Identify significant PCs

To overcome the extensive technical noise in the expression of any single gene for scRNA-seq data, **Seurat assigns cells to clusters based on their PCA scores derived from the expression of the integrated most variable genes**, with each PC essentially representing a "metagene" that combines information across a correlated gene set. **Determining how many PCs to include in the clustering step is therefore important to ensure that we are capturing the majority of the variation**, or cell types, present in our dataset. 

It is useful to explore the PCs prior to deciding which PCs to include for the downstream clustering analysis.

(a) One way of exploring the PCs is using a heatmap to visualize the most variant genes for select PCs with the **genes and cells ordered by PCA scores**. The idea here is to look at the PCs and determine whether the genes driving them make sense for differentiating the different cell types. 

The `cells` argument specifies the number of cells with the most negative or postive PCA scores to use for the plotting. The idea is that we are looking for a PC where the heatmap starts to look more "fuzzy", i.e. where the distinctions between the groups of genes is not so distinct.

```R
#| label: fig-integrated_DimHeatmap
#| fig-cap: Heatmap of the first 9 principal components, showing expression across 500 cells for the top genes (negative and positive) for the component.
# Explore heatmap of PCs
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
```

This method can be slow and hard to visualize individual genes if we would like to explore a large number of PCs. In the same vein and to explore a large number of PCs, we could print out the top 10 (or more) positive and negative genes by PCA scores driving the PCs.

```R
#| label: genes_driving_PCs
# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
```



(b) The **elbow plot** is another helpful way to determine how many PCs to use for clustering so that we are capturing majority of the variation in the data. The elbow plot visualizes the standard deviation of each PC, and we are looking for where the standard deviations begins to plateau. Essentially, **where the elbow appears is usually the threshold for identifying the majority of the variation**. However, this method can be quite subjective. 

Let's draw an elbow plot using the top 40 PCs:

```R
#| label: fig-PC_elbow_plot
#| fig-cap: Plot showing the standard deviation represented with each component to determine a cutoff value for the number of PCs.
# Plot the elbow plot
ElbowPlot(object = seurat_integrated, 
          ndims = 40)
```

Based on this plot, we could roughly determine the majority of the variation by where the elbow occurs around PC8 - PC10, or one could argue that it should be when the data points start to get close to the X-axis, PC30 or so. This gives us a very rough idea of the number of PCs needed to be included, we can extract the information visualized here in a [**more quantitative manner**](Aside_elbow_plot_metric.qmd), which may be a bit more reliable.  

While the above 2 methods were used a lot more with older methods from Seurat for normalization and identification of variable genes, they are no longer as important as they used to be. This is because the **SCTransform method is more accurate than older methods**.


>> **Why is selection of PCs more important for older methods?**
>
>> The older methods incorporated some technical sources of variation into some of the higher PCs, so selection of PCs was more important. SCTransform estimates the variance better and does not frequently include these sources of technical variation in the higher PCs. 
>
>> In theory, with SCTransform, the more PCs we choose the more variation is accounted for when performing the clustering, but it takes a lot longer to perform the clustering. Therefore for this analysis, we will use the **first 40 PCs** to generate the clusters. 





## Cluster the cells

Seurat uses a graph-based clustering approach using a K-nearest neighbor approach, and then attempts to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’ [[Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html)]. A nice in-depth description of clustering methods is provided in the [SVI Bioinformatics and Cellular Genomics Lab course](https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/clustering-and-cell-annotation.html).

### Find neighbors

The first step is to **construct a K-nearest neighbor (KNN) graph** based on the euclidean distance in PCA space. 

This is done in Seurat by using the `FindNeighbors()` function:

```R
#| label: FindNeighbors
# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                dims = 1:40)
```                               

### Find clusters

Next, Seurat will **iteratively group cells together with the goal of optimizing the standard modularity function**. 

We will use the `FindClusters()` function to perform the graph-based clustering. The `resolution` is an important argument that sets the "granularity" of the downstream clustering and will need to be optimized for every individual experiment.  For datasets of 3,000 - 5,000 cells, the `resolution` set between `0.4`-`1.4` generally yields good clustering. Increased resolution values lead to a greater number of clusters, which is often required for larger datasets. 

The `FindClusters()` function allows us to enter a series of resolutions and will calculate the "granularity" of the clustering. This is very helpful for testing which resolution works for moving forward without having to run the function for each resolution.

```R
#| label: FindClusters
# Determine the clusters for various resolutions                              
seurat_integrated <- FindClusters(object = seurat_integrated,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
```


## Visualize clusters of cells

To visualize the cell clusters, there are a few different dimensionality reduction techniques that can be helpful. The most popular methods include [t-distributed stochastic neighbor embedding (t-SNE)](https://kb.10xgenomics.com/hc/en-us/articles/217265066-What-is-t-Distributed-Stochastic-Neighbor-Embedding-t-SNE-) and [Uniform Manifold Approximation and Projection (UMAP)](https://umap-learn.readthedocs.io/en/latest/index.html) techniques. 

Both methods aim to place cells with similar local neighborhoods in high-dimensional space together in low-dimensional space. These methods will require you to input number of PCA dimentions to use for the visualization, we suggest using the same number of PCs as input to the clustering analysis. Here, we will proceed with the [UMAP method](https://umap-learn.readthedocs.io/en/latest/how_umap_works.html) for visualizing the clusters.


We can only visualize the results of one resolution setting at a time. If we look at the metadata of our Seurat object(`seurat_integrated@meta.data`), you should observe a separate column for each of the different resolutions calculated.

```R
#| label: evaluate_resolutions
#| eval: false
# Explore resolutions
seurat_integrated@meta.data %>% 
        View()
```

```R
#| label: tbl-evaluate_resolutions
#| tbl-cap: Table of `@meta.data` to see the results of clustering for each different resolution.
#| echo: false
# Explore resolutions
seurat_integrated@meta.data %>% 
  head(n = 5) %>% 
  remove_rownames() %>%
  relocate(cells) %>%
  DT::datatable() %>% 
  DT::formatStyle("cells", 
                  "white-space" = "nowrap")
```

To **choose a resolution to start with**, we often pick something in the middle of the range like 0.6 or 0.8. We will start with a resolution of 0.8 by assigning the identity of the clusters using the `Idents()` function.

```R
#| label: assign_ident
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"
```

Now, we can plot the UMAP to look at how cells cluster together at a resolution of 0.8:

> **Calculate UMAP code**
We calculate the UMAP coordinates in the last lesson with the following code:
>
> ```R
> #| label: UMAP_for_0.8_resolution
> #| eval: false
> ## Calculation of UMAP
> ## **DO NOT RUN** (calculated in the last lesson)
>
> seurat_integrated <- RunUMAP(seurat_integrated,
>                             reduction = "pca",
>                             dims = 1:40)
> ```


```R
#| label: fig-integrated_DimPlot_0.8
#| fig-cap: UMAP plot, representing each cell as a colored point corresponding with the cluster identified at resolution 0.8.
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
```


It can be useful to **explore other resolutions as well**. It will give you a quick idea about how the clusters would change based on the resolution parameter. For example, let's switch to a resolution of 0.4:

```R
#| label: fig-integrated_DimPlot_0.4
#| fig-cap: UMAP plot, representing each cell as a colored point corresponding with the cluster identified at resolution 0.4.

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
```


**How does your UMAP plot compare to the one above?**

It is possible that there is some variability in the way your clusters look compared to the image in this lesson. In particular **you may see a difference in the labeling of clusters**. This is an unfortunate consequence of slight variations in the versions of packages (mostly Seurat dependencies).


**If your clusters do look different from what we have in the lesson**, please follow the instructions provided below. 

Inside your `data` folder you will see a folder called `additional_data`. It contains the seurat_integrated object that we have created for the class. 

* **Load in the object to your R session and overwrite the existing one**: 

```R
#| label: load_data
#| eval: false
load(bzfile("data/seurat_integrated.RData.bz2"))
```


We will now continue with the 0.8 resolution to check the quality control metrics and known markers for the anticipated cell types. Plot the UMAP again to make sure your image now (or still) matches what you see in the lesson:

```R
#| label: fig-reassign_0.8_ident_and_DimPlot
#| fig-cap: UMAP plot, representing each cell as a colored point corresponding with the cluster identified at resolution 0.8 for the downloaded, processed data.

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
```

---
<br>
<br>


# Clustering Quality Control

## Learning Objectives:

- Evaluate whether clustering artifacts are present 
- Determine the quality of clustering with PCA and UMAP plots, and decide when to re-cluster
- Assess known cell type markers to hypothesize cell type identities of clusters


```R
#| label: load_libraries
#| echo: false
# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(cowplot)

load(bzfile("data/seurat_integrated.RData.bz2"))
```


Now that we have performed the integration, we want to know the different cell types present within our population of cells. 


## Exploration of quality control metrics

To determine whether our clusters might be due to artifacts such as cell cycle phase or mitochondrial expression, it can be useful to explore these metrics visually to see if any clusters exhibit enrichment or are different from the other clusters. However, if enrichment or differences are observed for particular clusters it may not be worrisome if it can be explained by the cell type. 

To explore and visualize the various quality metrics, we will use the versatile `DimPlot()` and `FeaturePlot()` functions from Seurat. 

### Segregation of clusters by sample

We can start by exploring the distribution of cells per cluster in each sample:

```R
#| label: fig-cells_per_cluster_barplot
#| fig-cap: Barplot representing the number of cells found in each cluster, split by sample identity.
#| fig.width: 15
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample")) %>%
           dplyr::count(ident, sample)

# Barplot of number of cells per cluster by sample
ggplot(n_cells, aes(x=ident, y=n, fill=sample)) +
    geom_bar(position=position_dodge(), stat="identity") +
    theme_classic() +
    geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))
```


We can visualize the cells per cluster for each sample using the UMAP:

```R
#| label: fig-cells_per_cluster_UMAP
#| fig-cap: UMAP `split.by` sample identity and colored by cluster identity.
#| fig.width: 10
# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()
```


Additionally, we can supply the metadata dataframe from our seurat object into ggplot to create more visuals. Looking at a UMAP is a great way to get a first pass look at your dataset, but we encourage you to look at your data in multiple different ways. For example, looking at the proportion of cells from a sample in each cluster.

```R
#| label: fig-cell_proportion_barplot
#| fig-cap: Barplot representing the proportion of cells from each samples for every cluster.
# Barplot of proportion of cells in each cluster by sample
ggplot(seurat_integrated@meta.data) +
    geom_bar(aes(x=integrated_snn_res.0.8, fill=sample), 
             position=position_fill())  +
    theme_classic()
```

Generally, we expect to see the majority of the cell type clusters to be present in all conditions; however, depending on the experiment we might expect to see some condition-specific cell types present. These clusters look pretty similar between conditions, which is good since we expected similar cell types to be present in both control and stimulated conditions.

### Segregation of clusters by cell cycle phase

Next, we can explore whether the **cells cluster by the different cell cycle phases**. We did not regress out variation due to cell cycle phase when we performed the SCTransform normalization and regression of uninteresting sources of variation. If our cell clusters showed large differences in cell cycle expression, this would be an indication we would want to re-run the SCTransform and add the `S.Score` and `G2M.Score` to our variables to regress, then re-run the rest of the steps.


```R
#| label: fig-clusters_DimPlot
#| fig-cap: UMAP `split.by` cell cycle phase and colored by cluster identity.
#| fig.width: 15
# Explore whether clusters segregate by cell cycle phase
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()
```

We do not see much clustering by cell cycle score, so we can proceed with the QC.

### Segregation of clusters by various sources of uninteresting variation

Next we will explore additional metrics, such as the number of UMIs and genes per cell, S-phase and G2M-phase markers, and mitochondrial gene expression by UMAP. Looking at the individual S and G2M scores can give us additional information to checking the phase as we did previously.

```R
#| label: fig-cell_phase_FeaturePlot
#| fig-cap: UMAP `FeaturePlot` showing a variety of QC metrics with cluster annotated on top.
#| fig.width: 9
#| fig.height: 9
# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
```


**`order` argument** : The `order` argument will plot the positive cells above the negative cells, while the `min.cutoff` argument will determine the threshold for shading. A `min.cutoff` of `q10` translates to the 10% of cells with the lowest expression of the gene will not exhibit any purple shading (completely gray).


The metrics seem to be relatively even across the clusters, with the exception of `nGene` exhibiting slightly higher values in clusters to the left of the plot. This can be more clearly seen when we look at the distribution as a boxplot. We can keep an eye on these clusters to see whether the cell types may explain the increase.

```R
#| label: fig-genes_per_cluster_boxplot
#| fig-cap: Boxplot distribution of nGenes for each cluster.
# Boxplot of nGene per cluster
ggplot(seurat_integrated@meta.data) +
    geom_boxplot(aes(x=integrated_snn_res.0.8, y=nGene, 
                     fill=integrated_snn_res.0.8)) +
    theme_classic() +
    NoLegend()
```

If we see differences corresponding to any of these metrics at this point in time, then we will often note them and then decide after identifying the cell type identities whether to take any further action.

### Exploration of the PCs driving the different clusters

We can also explore how well our clusters separate by the different PCs; we hope that the defined PCs separate the cell types well. To visualize this information, we need to extract the UMAP coordinate information for the cells along with their corresponding scores for each of the PCs to view by UMAP. 

First, we identify the information we would like to extract from the Seurat object, then, we can use the `FetchData()` function to extract it.

```R
#| label: FetchData_seurat_integrated
#| eval: false
# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
            "ident",
            "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)
head(pc_data)
```

```R
#| label: tbl-FetchData_seurat_integrated
#| tbl-cap: Principal component scores for the first several cells of the dataset.
#| echo: false
# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
            "ident",
            "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)
pc_data %>%
  head() %>%
  knitr::kable()
```

***# The `FetchData()` function***
How did we know in the `FetchData()` function to include `UMAP_1` to obtain the UMAP coordinates? The [Seurat cheatsheet](12_seurat_cheatsheet.qmd) describes the function as being able to pull any data from the expression matrices, cell embeddings, or metadata. 
 
For instance, if you explore the `seurat_integrated@reductions` list object, the first component is for PCA, and includes a slot for `cell.embeddings`. We can use the column names (`PC_1`, `PC_2`, `PC_3`, etc.) to pull out the coordinates or PC scores corresponding to each cell for each of the PCs. 

We could do the same thing for UMAP:

```R
#| label: extract_UMAP_coordinates
# Extract the UMAP coordinates for the first 10 cells
seurat_integrated@reductions$umap@cell.embeddings[1:10, 1:2]
```

The `FetchData()` function just allows us to extract the data more easily.


***# Seurat version < 5.0*** The pre-existing `seurat_integrated` loaded in previously was created using an older version of Seurat. As such the columns we `Fetch()` are in upper case (i.e `UMAP_1`). **If you are using your own seurat object using a newer version of Seurat you will need to change the column names as shown below.** Alternatively, explore your Seurat object to see how they have been stored.

```R
#| label: create_UMAP_columns
# Defining the information in the Seurat object of interest
columns <- c(paste0("PC_", 1:16),
         "ident",
         "umap_1", "umap_2")
columns
```


In the UMAP plots, the cells are colored by their PC score for each respective principal component. 

Let's take a quick look at the top 16 PCs:

```R
#| label: fig-UMAP_per_PC
#| fig-cap: UMAP scatterplot of top 16 PCs, coloring each cell by its score in the PC.
#| fig.height: 10
#| fig.width: 10
#| warning: false
# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  dplyr::summarise(x=mean(UMAP_1), y=mean(UMAP_2))
  
# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
        ggplot(pc_data, 
               aes(UMAP_1, UMAP_2)) +
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
```


We can see how the clusters are represented by the different PCs. For instance, the genes driving `PC_2` exhibit higher expression in clusters 8 and 12. We could look back at our genes driving this PC to get an idea of what the cell types might be:

```R
#| label: genes_driving_PCA
# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)
```

With the `GNLY` and `NKG7` genes as positive markers of `PC_2`, we can hypothesize that clusters 8 and 12 correspond to NK cells. This just hints at what the clusters identity could be, with the identities of the clusters being determined through a combination of the PCs. 

To truly determine the identity of the clusters and whether the `resolution` is appropriate, it is helpful to explore a handful of known gene markers for the cell types expected. 

## Exploring known cell type markers

With the cells clustered, we can explore the cell type identities by looking for known markers. The UMAP plot with clusters marked is shown, followed by the different cell types expected.

```R
#| label: fig-DimPlot_clusters
#| fig-cap: UMAP plot, representing each cell as a colored point corresponding with the cluster identified at resolution 0.8.
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE) + NoLegend()
```


| Cell Type | Marker |
|:---:|:---:|
| CD14+ monocytes | CD14, LYZ | 
| FCGR3A+ monocytes | FCGR3A, MS4A7 |
| Conventional dendritic cells | FCER1A, CST3 |
| Plasmacytoid dendritic cells | IL3RA, GZMB, SERPINF1, ITM2C |
| B cells | CD79A, MS4A1 |
| T cells | CD3D |
| CD4+ T cells | CD3D, IL7R, CCR7 |
| CD8+ T cells| CD3D, CD8A |
| NK cells | GNLY, NKG7 |
| Megakaryocytes | PPBP |
| Erythrocytes | HBB, HBA2 |

The `FeaturePlot()` function from seurat makes it easy to visualize a handful of genes using the gene IDs stored in the Seurat object. We can easily explore the expression of known gene markers on top of our UMAP visualizations. Let's go through and determine the identities of the clusters. To access the normalized expression levels of all genes, we can use the normalized count data stored in the `RNA` assay slot. 

---
**# SCTransform dimensions**
The SCTransform normalization and integration was performed only on the 3000 most variable genes, so many of our genes of interest may not be present in this data. 
```R
dim(seurat_integrated[["RNA"]])
dim(seurat_integrated[["integrated"]])
```
---

```R
#| label: NormalizeData_integrated
# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
seurat_integrated
```

---

***Assays***
Assay is a slot defined in the Seurat object, it has multiple slots within it. In a given assay, the `counts` slot stores non-normalized raw counts, and the `data` slot stores normalized expression data. Therefore, when we run the `NormalizeData()` function in the above code, the normalized data will be stored in the `data` slot of the RNA assay while the `counts` slot will remain unaltered.

---

Depending on our markers of interest, they could be positive or negative markers for a particular cell type. The combined expression of our chosen handful of markers should give us an idea on whether a cluster corresponds to that particular cell type. 

For the markers used here, we are looking for positive markers and consistency of expression of the markers across the clusters. For example, if there are two markers for a cell type and only one of them is expressed in a cluster - then we cannot reliably assign that cluster to the cell type.


**CD14+ monocyte markers**

```R
#| label: fig-CD4_FeaturePlot
#| fig-cap: UMAP `FeaturePlot()` of top CD14+ monocyte markers.
#| fig.width: 8
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD14", "LYZ"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
```

CD14+ monocytes appear to correspond to clusters 1, and 3. We wouldn't include clusters 14 and 10 because they do not highly express both of these markers.

**FCGR3A+ monocyte markers**

```R
#| label: fig-FCGR3A_monocyte_plot
#| fig-cap: UMAP `FeaturePlot()` of top FCGR3A+ monocyte markers.
#| fig.width: 8
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("FCGR3A", "MS4A7"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
```


FCGR3A+ monocytes markers distinctly highlight cluster 10, although we do see some decent expression in clusters 1 and 3 We would like to see additional markers for FCGR3A+ cells show up when we perform the marker identification.

**Macrophages**

```R
#| label: fig-Macrophages_plot
#| fig-cap: UMAP `FeaturePlot()` of top FCGR3A+ Macrophages markers.
#| fig.width: 8
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("MARCO", "ITGAM", "ADGRE1"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
```

We don't see much overlap of our markers, so no clusters appear to correspond to macrophages; perhaps cell culture conditions negatively selected for macrophages (more highly adherent).

To clearly see the expression levels of these genes, we can generate a `VlnPlot()` of these genes to more clearly see the distribution of expression.

```R
#| label: fig-Macrophages_plot_vln
#| fig-cap: Violin plot of top FCGR3A+ Macrophages markers.
#| fig.height: 6
VlnPlot(seurat_integrated,
        c("MARCO", "ITGAM", "ADGRE1"),
        ncol = 1)
```

**Conventional dendritic cell markers**

```R
#| label: fig-conventional_dendritic cell_plot
#| fig-cap: UMAP `FeaturePlot()` of top Conventional dendritic cell markers.
#| fig-width: 8
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("FCER1A", "CST3"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
```

The markers corresponding to conventional dendritic cells identify cluster 14 (both markers consistently show expression).

**Plasmacytoid dendritic cell markers**

```R
#| label: fig-plasmacytoid_dendritic_cell_plot
#| fig-cap: UMAP `FeaturePlot()` of top Plasmacytoid dendritic cell markers.
#| fig.width: 8
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
```

Plasmacytoid dendritic cells represent cluster 16. While there are a lot of differences in the expression of these markers, we see cluster 16 (though small) is consistently strongly expressed.

***

Seurat also has a built in visualization tool which allows us to view the average expression of genes across clusters called `DotPlot()`. This function additionally shows us how many cells within the cluster have expression of one gene. As input, we supply a list of genes - note that we cannot use the same gene twice or an error will be thrown. 

```R
#| label: fig-DotPlot_known_markers
#| fig-cap: DotPlot representing top marker genes for a variety of celltypes, with each circle reprsenting the average expression of that cluster and the size showing the percentage of cells that express that gene.
#| fig.width: 15
# List of known celltype markers
markers <- list()
markers[["CD14+ monocytes"]] <- c("CD14", "LYZ")
markers[["FCGR3A+ monocyte"]] <- c("FCGR3A", "MS4A7")
markers[["Macrophages"]] <- c("MARCO", "ITGAM", "ADGRE1")
markers[["Conventional dendritic"]] <- c("FCER1A", "CST3")
markers[["Plasmacytoid dendritic"]] <- c("IL3RA", "GZMB", "SERPINF1", "ITM2C")

# Create dotplot based on RNA expression
DotPlot(seurat_integrated, markers, assay="RNA")
```

***

[**Exercise 1**](11_clustering_quality_control-Answer_key.qmd#exercise-1)
Hypothesize the clusters corresponding to each of the different clusters in the table:

| Cell Type | Clusters |
|:---:|:---:|
| CD14+ monocytes | 1, 3 | 
| FCGR3A+ monocytes | 10 |
| Conventional dendritic cells | 14 |
| Plasmacytoid dendritic cells | 16 |
| Marcrophages | - |
| B cells | ? |
| T cells | ? |
| CD4+ T cells | ? |
| CD8+ T cells| ? |
| NK cells | ? |
| Megakaryocytes | ? |
| Erythrocytes | ? |
| Unknown | ? |

***

If any cluster appears to contain two separate cell types, it's helpful to increase our clustering resolution to properly subset the clusters. Alternatively, if we still can't separate out the clusters using increased resolution, then it's possible that we had used too few principal components such that we are just not separating out these cell types of interest. To inform our choice of PCs, we could look at our PC gene expression overlapping the UMAP plots and determine whether our cell populations are separating by the PCs included.


Now we have a decent idea as to the cell types corresponding to the majority of the clusters, but some questions remain:

1. *T cell markers appear to be highly expressed in may clusters. How can we differentiate and subset the larger group into smaller subset of cells?*
2. *Do the clusters corresponding to the same cell types have biologically meaningful differences? Are there subpopulations of these cell types?*
3. *Can we acquire higher confidence in these cell type identities by identifying other marker genes for these clusters?*

Marker identification analysis can help us address all of these questions!! 

The next step will be to perform marker identification analysis, which will output the genes that significantly differ in expression between clusters. Using these genes we can determine or improve confidence in the identities of the clusters/subclusters.



---
<br>
<br>

# Seurat Cheatsheet

This part provides participants with a Seurat cheatsheet that demonstrates commonly used functions for exploring, wrangling and visualizing single-cell RNA-seq data. Participants will learn how to access metadata, extract features and cell barcodes, work with assays and layers, retrieve dimensional reduction coordinates and generate a variety of Seurat-based plots. Each function is presented with examples using the workshop dataset so that participants can understand what each command does and how to apply it in their own analyses.


---
<br>
<br>


# Marker identification

## Learning Objectives:

- Describe how to determine markers of individual clusters
- Discuss the iterative processes of clustering and marker identification




```R
#| label: load_libraries
#| echo: false
# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(knitr)
library(DT)

load(bzfile("data/seurat_integrated.RData.bz2"))
```

Now that we have identified our desired clusters, we can move on to marker identification, which will allow us to verify the identity of certain clusters and help surmise the identity of any unknown clusters.


**Remember that we had the following questions from the clustering analysis**: 

1. _Do the clusters corresponding to the same cell types have biologically meaningful differences? Are there subpopulations of these cell types?_
2. _Can we acquire higher confidence in these cell type identities by identifying other marker genes for these clusters?_

There are a few different types of marker identification that we can explore using Seurat to get to the answer of these questions. Each with their own benefits and drawbacks:

1. **Identification of all markers for each cluster:** this analysis compares each cluster against all others and outputs the genes that are differentially expressed/present. 
	- *Useful for identifying unknown clusters and improving confidence in hypothesized cell types.*

2. **Identification of conserved markers for each cluster:** This analysis looks for genes that are differentially expressed/present within each condition first, and then reports those genes that are conserved in the cluster across all conditions. These genes can help to figure out the identity for the cluster. 
	- *Useful with more than one condition to identify cell type markers that are conserved across conditions.*  	
 
3. **Marker identification between specific clusters:** this analysis explores differentially expressed genes between specific clusters. 
	- *Useful for determining differences in gene expression between clusters that appear to be representing the same celltype (i.e with markers that are similar) from the above analyses.*



## Identification of all markers for each cluster
>> **## Dont Run this step**

This type of analysis is typically **recommended for when evaluating a single sample group/condition**. With the ` FindAllMarkers()` function we are comparing each cluster against all other clusters to identify potential marker genes. The cells in each cluster are treated as replicates, and essentially a **differential expression analysis is performed** with some statistical test. 

**`test.use` parameter** : The default is a Wilcoxon Rank Sum test, but there are other options available. 


Schematic of identifying markers for each cluster by comparing one cluster against all other cells (`FindAllMarkers`).


The `FindAllMarkers()` function has **three important arguments** which provide thresholds for determining whether a gene is a marker:

| Argument | Description | Cons |
|---|---|---|
| `logfc.threshold` | minimum log2 fold change for average expression of gene in cluster relative to the average expression in all other clusters combined. Default is 0.25. | - Could miss those cell markers that are expressed in a small fraction of cells within the cluster of interest, but not in the other clusters, if the average logfc doesn't meet the threshold<br>- Could return a lot of metabolic/ribosomal genes due to slight differences in metabolic output by different cell types, which are not as useful to distinguish cell type identities |
| `min.diff.pct` | minimum percent difference between the percent of cells expressing the gene in the cluster and the percent of cells expressing gene in all other clusters combined. | Could miss those cell markers that are expressed in all cells, but are highly up-regulated in this specific cell type |
| `min.pct` | only test genes that are detected in a minimum fraction of cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1. | If set to a very high value could incur many false negatives due to the fact that not all genes are detected in all cells (even if it is expressed) |
	
You could use any combination of these arguments depending on how stringent/lenient you want to be. Also, by default this function will return to you genes that exhibit both positive and negative expression changes. Typically, we add an argument `only.pos` to opt for keeping only the positive changes. The code to find markers for each cluster is shown below. **We will not run this code.**

```R
#| label: FindAllMarkers_example
#| eval: false
## **DO NOT RUN**

# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(object = seurat_integrated, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)                     
```

**# Why we are not running `FindAllMarkers()`** This command can take quite long to run, as it is processing each individual cluster against all other cells.


## Identification of conserved markers in all conditions
>> **## Run me**

Since we have samples representing different conditions in our dataset, **our best option is to find conserved markers**. This function internally separates out cells by sample group/condition, and then performs differential gene expression testing for a single specified cluster against all other clusters (or a second cluster, if specified). Gene-level p-values are computed for each condition and then combined across groups using meta-analysis methods from the MetaDE R package.


Before we start our marker identification we will explicitly set our default assay, we want to use the **normalized data, but not the integrated data**.

```R
#| label: set_default_assay
DefaultAssay(seurat_integrated) <- "RNA"
Idents(seurat_integrated) <- "integrated_snn_res.0.8"
```

The default assay should have already been `RNA`, because we set it up in the previous clustering quality control lesson. But we encourage you to run this line of code above to be absolutely sure in case the active slot was changed somewhere upstream in your analysis. 

>> ** Why don't we use SCT normalized data?**

>>Note that the raw and normalized counts are stored in the `counts` and `data` slots of `RNA` assay, respectively. By default, the functions for finding markers will use normalized data if RNA is the DefaultAssay. The number of features in the `RNA` assay corresponds to all genes in our dataset.
<br>
>> Now if we consider the `SCT` assay, functions for finding markers would use the `scale.data` slot which is the pearson residuals that come out of regularized NB regression. Differential expression on these values can be difficult interpret. Additionally, only the variable features are represented in this assay and so we may not have data for some of our marker genes.


The function `FindConservedMarkers()`, has the following structure:

**`FindConservedMarkers()` syntax:**

```R
#| label: FindConservedMarkers_example
#| eval: false
## **DO NOT RUN** ##
FindConservedMarkers(seurat_integrated,
                      ident.1 = cluster,
                      grouping.var = "sample",
                      only.pos = TRUE,
                      min.diff.pct = 0.25,
                      min.pct = 0.25,
                      logfc.threshold = 0.25)
```

You will recognize some of the arguments we described previously for the `FindAllMarkers()` function; this is because internally it is using that function to first find markers within each group. Here, we list **some additional arguments** which provide for when using `FindConservedMarkers()`:

- `ident.1`: this function only evaluates one cluster at a time; here you would specify the cluster of interest.
- `grouping.var`: the variable (column header) in your metadata which specifies the separation of cells into groups

For our analysis we will be fairly lenient and **use only the log fold change threshold greater than 0.25**. We will also specify to return only the positive markers for each cluster.


Let's **test it out on one cluster** to see how it works:

```R
#| label: FindConservedMarkers_cluster_0
#| eval: false
cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                    ident.1 = 0,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)

# Inspect the output of FindConservedMarkers
View(cluster0_conserved_markers)
```

```R
#| label: tbl-_FindConservedMarkers_cluster_0
#| tbl-cap: "Output from FindConservedMarkers"
#| echo: false
cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                    ident.1 = 0,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)

cluster0_conserved_markers %>% 
  head(n = 10) %>% 
  kable()
```

**The output from the `FindConservedMarkers()` function**, is a matrix containing a ranked list of putative markers listed by gene ID for the cluster we specified, and associated statistics. Note that the same set of statistics are computed for each group (in our case, Ctrl and Stim) and the last two columns correspond to the combined p-value across the two groups. We describe some of these columns below:

- **gene:** gene symbol
- **condition_p_val:** p-value not adjusted for multiple test correction for condition
- **condition_avg_logFC:** average log fold change for condition. Positive values indicate that the gene is more highly expressed in the cluster.	
- **condition_pct.1:** percentage of cells where the gene is detected in the cluster for condition		
- **condition_pct.2:** percentage of cells where the gene is detected on average in the other clusters for condition
- **condition_p_val_adj:** adjusted p-value for condition, based on bonferroni correction using all genes in the dataset, used to determine significance
- **max_pval:**	largest p value of p value calculated by each group/condition
- **minimump_p_val:** combined p value

**Inflated p-values** Since each cell is being treated as a replicate this will result in inflated p-values within each group! A gene may have an incredibly low p-value < 1e-50 but that doesn't translate as a highly reliable marker gene. 


When looking at the output, **we suggest looking for markers with large differences in expression between `pct.1` and `pct.2` and larger fold changes**. For instance if `pct.1` = 0.90 and `pct.2` = 0.80, it may not be as exciting of a marker. However, if `pct.2` = 0.1 instead, the bigger difference would be more convincing. Also, of interest is if the majority of cells expressing the marker is in my cluster of interest. If `pct.1` is low, such as 0.3, it may not be as interesting. Both of these are also possible parameters to include when running the function, as described above.


### Adding Gene Annotations

It can be helpful to add columns with gene annotation information. In order to do that we will load in an annotation file located in your `data` folder, using the code provided below: 

```R
#| label: read_in_annotations
annotations <- read.csv("data/annotation.csv")
```

**# How to create annotation file**
If you are interested in knowing how we obtained this annotation file, take a look at [the linked materials](Aside_fetching_annotations.qmd).

First, we will turn the row names with gene identifiers into its own columns. Then we will merge this annotation file with our results from the `FindConservedMarkers()`:

```R
#| label: wrangling_cluster_0_annotations
#| eval: false
# Combine markers with gene descriptions 
cluster0_ann_markers <- cluster0_conserved_markers %>% 
                rownames_to_column(var="gene") %>% 
                left_join(y = unique(annotations[, c("gene_name", "description")]),
                          by = c("gene" = "gene_name"))

head(cluster0_ann_markers)
```


```R
#| label: tbl-wrangling_cluster_0_annotations
#| tbl-cap: Output from `FindConservedMarkers` after merging with annotation information.
#| eval: false
# Combine markers with gene descriptions 
cluster0_ann_markers <- cluster0_conserved_markers %>% 
                rownames_to_column(var="gene") %>% 
                left_join(y = unique(annotations[, c("gene_name", "description")]),
                          by = c("gene" = "gene_name"))

cluster0_ann_markers %>%
  head() %>%
  knitr::kable()
```

**# Exercises**

1. In the previous lesson, we identified cluster 10 as FCGR3A+ monocytes by inspecting the expression of known cell markers FCGR3A and MS4A7. Use `FindConservedMarkers()` function to find conserved markers for cluster 10. What do you observe? Do you see FCGR3A and MS4A7 as highly expressed genes in cluster 10?


### Running on multiple samples

The function `FindConservedMarkers()` **accepts a single cluster at a time**, and we could run this function as many times as we have clusters. However, this is not very efficient. Instead we will first create a function to find the conserved markers including all the parameters we want to include. We will also **add a few lines of code to modify the output**. Our function will:

1. Run the `FindConservedMarkers()` function
2. Transfer row names to a column using `rownames_to_column()` function
3. Merge in annotations
3. Create the column of cluster IDs using the `cbind()` function


```R
#| label: create_get_cluster_function
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
```

Now that we have this function created  we can use it as an argument to the appropriate `map` function. We want the output of the `map` family of functions to be a **dataframe with each cluster output bound together by rows**, we will use the `map_dfr()` function.

**`map` family syntax:**

```R
#| label: map_dfr_example
#| eval: false
## **DO NOT RUN** ##
map_dfr(inputs_to_function, name_of_function)
```

Now, let's try this function to **find the conserved markers for the clusters that were identified as CD4+ T cells (4,0,6,2)** from our use of known marker genes. Let's see what genes we identify and of there are overlaps or obvious differences that can help us tease this apart a bit more.

```R
#| label: iterate_conserved_markers_across_select_clusters
#| eval: false
# Iterate function across desired clusters
conserved_markers <- map_dfr(c(4,0,6,2), get_conserved)
head(conserved_markers)
```

```R
#| label: tbl-iterate_conserved_markers_across_select_clusters
#| tbl-cap: Output from `FindConservedMarkers` for multiple clusters after wranging with the custom function `get_conserved`.
# Iterate function across desired clusters
conserved_markers <- map_dfr(c(4,0,6,2), get_conserved)
conserved_markers %>%
  head() %>%
  DT::datatable() %>% 
  DT::formatStyle("description", 
                  "white-space" = "nowrap")
```

**# Finding markers for all clusters**

For your data, you may want to run this function on all clusters, in which case you could input `0:20` instead of `c(4,0,6,2)`; however, it would take quite a while to run. Also, it is possible that when you run this function on all clusters, in **some cases you will have clusters that do not have enough cells for a particular group** - and  your function will fail. For these clusters you will need to use `FindAllMarkers()`.

### Evaluating marker genes

We would like to use these gene lists to see of we can **identify which celltypes these clusters identify with.** Let's take a look at the top genes for each of the clusters and see if that gives us any hints. We can view the top 10 markers by average fold change across the two groups, for each cluster for a quick perusal:

```R
#| label: extract_top_markers
# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (ctrl_avg_log2FC + stim_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)
```


```R
#| label: view_top_markers
#| eval: false
# Visualize top 10 markers per cluster
View(top10)
```

```R
#| label: tbl_top_markers
#| tbl-cap: "Inspecting the top 10 genes for each cluster."
#| echo: false
# Visualize top 10 markers per cluster
top10 %>%
  head() %>%
  DT::datatable() %>% 
  DT::formatStyle("description", 
                  "white-space" = "nowrap")
```

When we look at the entire list, we see clusters 0 and 6 have some overlapping genes, like CCR7 and SELL which correspond to **markers of memory T cells**. It is possible that these two clusters are more similar to one another and could be merged together as naive T cells. On the other hand, with cluster 2 we observe CREM as one of our top genes; a **marker gene of activation**. This suggests that perhaps cluster 2 represents activated T cells.

| Cell State | Marker |
|:---:|:---:|
| Naive T cells | CCR7, SELL | 
| Activated T cells | CREM, CD69 |

For cluster 4, we see a lot of heat shock and DNA damage genes appear in the top gene list. Based on these markers, it is likely that these are **stressed or dying cells**. However, if we explore the quality metrics for these cells in more detail (i.e. mitoRatio and nUMI overlayed on the cluster) we don't really support for this argument. There is a breadth of research supporting the association of heat shock proteins with reactive T cells in the induction of anti‐inflammatory cytokines in chronic inflammation. This is a cluster for which we would need a deeper understanding of immune cells to really tease apart the results and make a final conclusion.

### Visualizing marker genes

To get a better idea of cell type identity for **cluster 4** we can **explore the expression of different identified markers** by cluster using the `FeaturePlot()` function. We see that only a subset of cluster 4 are highly expressing these genes.

```R
#| label: fig-cluster_4_FeaturePlot
#| fig-cap: UMAP `FeaturePlot` showing expression across cells for the top genes in cluster 4.
#| fig.width: 8
# Plot interesting marker gene expression for cluster 4
FeaturePlot(object = seurat_integrated, 
            features = c("HSPH1", "HSPE1", "DNAJB1"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
```

We can also explore the range in expression of specific markers by using **violin plots**:

**# Violin plots**
Violin plots are similar to box plots, except that they also show the probability density of the data at different values, usually smoothed by a kernel density estimator. A violin plot is more informative than a plain box plot. While a box plot only shows summary statistics such as mean/median and interquartile ranges, the violin plot shows the full distribution of the data. The difference is particularly useful when the data distribution is multimodal (more than one peak). In this case a violin plot shows the presence of different peaks, their position and relative amplitude.

```R
#| label: fig-cluster_4_ViolinPlot
#| fig-cap: Violin plot showing expression across each cluster for the top genes in cluster 4.
#| fig.width: 15
# Vln plot - cluster 4
VlnPlot(object = seurat_integrated, 
        features = c("HSPH1", "HSPE1", "DNAJB1"))
```        

These results and plots can help us determine the identity of these clusters or verify what we hypothesize the identity to be after exploring the canonical markers of expected cell types previously.


## Identifying gene markers for each cluster

Sometimes the list of markers returned don't sufficiently separate some of the clusters. For instance, we had previously identified clusters 0, 4, 6 and 2 as CD4+ T cells, but when looking at marker gene lists we identfied markers to help us further subset cells. We were lucky and the signal observed from `FindAllMarkers()` helped us differentiate between naive and activated cells. Another option to identify biologically meaningful differences would be to use the **`FindMarkers()` function to determine the genes that are differentially expressed between two specific clusters**. 



We can try all combinations of comparisons, but we'll start with cluster 2 versus all other CD4+ T cell clusters:

```R
#| label: FindMarkers_T_cell
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
```

```R
#| label: show_t_cell_markers
#| eval: false
# View data
View(cd4_tcells)
```

```R
#| label: tbl-_t_cell_markers_table
#| tbl-cap: "Inspecting FindMarkers output"
#| echo: false
cd4_tcells %>% 
  head() %>%
  DT::datatable() %>% 
  DT::formatStyle("description", 
                  "white-space" = "nowrap")
```

Of these top genes the **CREM gene** stands out as a marker of activation with a positive fold change.  We also see markers of naive or memory cells include the SELL and CCR7 genes with negative fold changes, which is in line with previous results. 

As markers for the naive and activated states both showed up in the marker list, it is helpful to visualize expression. Based on these plots it seems as though clusters 0 and 2 are reliably the naive T cells. However, for the activated T cells it is hard to tell. We might say that clusters 4 and 18 are activated T cells, but the CD69 expression is not as apparent as CREM. We will label the naive cells and leave the remaining clusters labeled as CD4+ T cells.

Now taking all of this information, we can surmise the cell types of the different clusters and plot the cells with cell type labels.


| Cluster ID	| Cell Type |
|:-----:|:-----:|
|0	| Naive or memory CD4+ T cells|
|1	| CD14+ monocytes |
|2	| Activated T cells|
|3	| CD14+ monocytes|
|4	| Stressed cells / Unknown|
|5	| CD8+ T cells |
|6	| Naive or memory CD4+ T cells |
|7	| B cells |
|8	| NK cells |
|9	| CD8+ T cells |
|10	| FCGR3A+ monocytes |
|11	| B cells |
|12	| NK cells |
|13	| B cells|
|14	| Conventional dendritic cells |
|15| Megakaryocytes |
|16| Plasmacytoid dendritic cells |


We can then reassign the identity of the clusters to these cell types:

```R
#| label: fig-rename_clusters
#| fig-cap: UMAP visualization with each cell colored by celltype annotation.
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
```


If we wanted to remove the potentially stressed cells, we could use the `subset()` function:

```R
#| label: removed_dying_cells_DimPlot
#| fig-cap: UMAP visualization with each cell colored by celltype annotation after removing stressed/dying cells.
# Remove the stressed or dying cells
seurat_subset_labeled <- subset(seurat_integrated,
                               idents = "Stressed cells / Unknown", invert = TRUE)

# Re-visualize the clusters
DimPlot(object = seurat_subset_labeled, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
	repel = TRUE)
```


Now we would want to save our final labelled Seurat object and the output of `sessionInfo()`:

```R      
#| eval: false
#| label: save_output
# Save final R object
write_rds(seurat_integrated,
          file = "results/seurat_labelled.rds")

# Create and save a text file with sessionInfo
sink("results/sessionInfo_scrnaseq.txt")
sessionInfo()
sink()
```

***

# Downstream Analyses
Now that we have our clusters defined and the markers for each of our clusters, we have a few different questions we can answer:

- Determine if there is a shift in cell populations between `ctrl` and `stim`. Ideally this would be done with replicates to determine if the changes are significant.

```R
#| label: fig-cells_per_celltype_plot
#| fig-cap: Example of how to identify shifts in celltype proportion by condition.
#| fig.width: 18
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
```

- Perform differential expression analysis between conditions `ctrl` and `stim`. We can use the `FindMarkers()` function to do a simple wilcox test to see the difference in gene expression between conditions for the B cells

```R
#| label: FindMarkers_B_cells
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
```

  - For added visualization, we can used the `EnhancedVolcano()` function to see how the genes fall on a volcano plot.

```R
#| label: fig-B_cells_volcano_plot
#| fig-cap: Volcano plot showing the differences in p-value and logFoldChange between conditions `stim` and `ctrl` from the `FindMarkers` analysis.
#| fig.height: 8
#| fig.width: 10
library(EnhancedVolcano)
EnhancedVolcano(b_markers,
    row.names(b_markers),
    x="avg_log2FC",
    y="p_val_adj",
    title="B Cells",
    subtitle="Stim vs. Ctrl"
)
```

- Pseudobulk DESEq2 differential gene analysis. Biological replicates are **necessary** to proceed with this analysis, and we have [additional materials](https://hbctraining.github.io/Pseudobulk-for-scRNAseq/) to help walk through this analysis.
- [Pathway analysis](https://hbctraining.github.io/Pseudobulk-for-scRNAseq-Quarto/lessons/10_functional_analysis_pseudobulk.html) with GSEA and over-representation analysis (ORA)
- Experimentally validate intriguing markers for our identified cell types.
- Explore a subset of the cell types to discover subclusters of cells as described [here](Aside_seurat_subclustering.qmd)
- Trajectory analysis, or lineage tracing, could be performed if trying to determine the progression between cell types or cell states. For example, we could explore any of the following using this type of analysis:
	- Differentiation processes
	- Expression changes over time
	- Cell state changes in expression


---

