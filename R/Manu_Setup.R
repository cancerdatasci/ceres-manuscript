
### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("ProjectTemplate",
           "doMC",
           "devtools",
           "rfigshare",
           "Matrix",
           "matrixStats",
           "scales",
           "magrittr",
           "tibble",
           "stringr",
           "readr",
           "plyr",
           "tidyr",
           "dplyr",
           "ggplot2",
           "grid",
           "cowplot",
           "squash",
           "broom",
           "ggrepel")

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}


bioc_pkgs <- c("Biostrings", "Rsamtools", "GenomeInfoDb",
               "GenomicRanges", "BSgenome", "BSgenome.Hsapiens.UCSC.hg19",
               "Gviz")

source("https://bioconductor.org/biocLite.R")
biocLite(pkgs = setdiff(bioc_pkgs, installed_pkgs),
         ask = F)

devtools::install_github("cancerdatasci/ceres")

### Download data from figshare

library(rfigshare)

ceres_figshare_id <- 5319388
data_dir <- "./data/raw"
dir.create(data_dir, showWarnings = F)

figshare_article <- fs_details(ceres_figshare_id)

for (figshare_file in figshare_article$files) {
    cat("downloading", figshare_file$name, "\n")
    file_path <- file.path(data_dir, figshare_file$name)
    download.file(figshare_file$download_url, file_path)
}


library(ProjectTemplate)
load.project()


