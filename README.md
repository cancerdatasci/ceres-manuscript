# CERES

Code repository for reproducing figures in the publication describing CERES.

Meyers, Bryan, *et al*. **Computational correction of copy number effect improves specificity of CRISPR–Cas9 essentiality screens in cancer cells.** *Nature Genetics.* 2017. [Article](https://doi.org/10.1038/ng.3984)

## Installation

We recommend running this code on a machine with at least 16GB RAM.

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

As you may have noticed, this R project uses a directory structure and functions from the R package, [ProjectTemplate](http://projecttemplate.net/). The following command runs each R file in the `src/` directory. Each script roughly corresponds to one type of analysis in the paper, and can be run individually using `source()`. Many of the scripts take advantage of parallelization using the package `doMC`. The number of threads/cores used may be edited towards the top of each script which uses them. If you are using a machine with 8GB of RAM, we recommend only using a single core for all analyses. With 16GB of RAM, we regularly use 4-6 cores.

```
library(ProjectTemplate)
run.project()
```
