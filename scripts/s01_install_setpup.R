
# 1. Preparation of server or computer ----
## install one by one not all together ----

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

BiocManager::install("EnhancedVolcano",dependencies = TRUE, lib="/courses/software/R_libs")

# (2) Install the 4 packages listed below from Bioconductor using the the BiocManager::install() function.

library(BiocManager)
BiocManager::install("AnnotationHub",dependencies = TRUE, lib="/courses/software/R_libs")
BiocManager::install("ensembldb",dependencies = TRUE, lib="/courses/software/R_libs")
BiocManager::install("multtest",dependencies = TRUE, lib="/courses/software/R_libs")
BiocManager::install("glmGamPoi",dependencies = TRUE, lib="/courses/software/R_libs")

# (3) Install Presto from GitHub using the devtools::install_github() function:
library(devtools)
devtools::install_github("immunogenomics/presto")

## Library ----
# (4) Finally,check that all the packages were installed successfully by loading them one at a time using the library() function.

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

