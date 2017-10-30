# CERES

Code repository for recreating figures in the publication describing CERES.

Meyers, Bryan, *et al*. **Computational correction of copy number effect improves specificity of CRISPR–Cas9 essentiality screens in cancer cells.** *Nature Genetics.* 2017. [Article](https://doi.org/10.1038/ng.3984)


## Installation

### Install the CERES R package and all dependencies

Follow all installation instructions for the [CERES R package](https://github.com/cancerdatasci/ceres)

### Clone this repository

From the command line:

```
git clone https://github.com/cancerdatasci/ceres-manuscript
```

### Open an R session from this directory

This can be done a variety of ways.

```
cd ceres-manuscript
R
```

or if you use RStudio

```
cd ceres-manuscript
open CERES_Manuscript.Rproj
```

From the R console, run the Manu_Setup.R script. This will install any missing packages and download data files from the [figshare record](https://doi.org/10.6084/m9.figshare.5319388). Follow the instructions when the program asks for figshare authentication.

```
source("./R/Manu_Setup.R")
```


## Run CERES

```
source("./R/Man_RunCERES.R")
```

## Prepare Data For Analysis

```
source("./R/Manu_CacheDatasets.R")
```

## Run Analysis Files

```
library(ProjectTemplate)
run.project()
```
