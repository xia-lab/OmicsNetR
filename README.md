
# Description

**_OmicsNetR_** is the underlying R package synchronized with OmicsNet web server. It is designed for network-based multi-omics integration and systems-level interpretation. The R package is composed of R functions necessary for the web-server to perform network creation, trimming and analysis. 

Following installation and loading of _OmicsNetR_, users will be able to reproduce web server results from their local computers using the R command history downloaded from OmicsNet. Running the R functions will allow more flexibility and reproducibility.

Note - OmicsNetR is still under development - we cannot guarantee full functionality
# Installation

**Step 1. Install package dependencies**

To use OmcisNetR, make sure your R version is >4.0.3 and install all package dependencies. Ensure that you are able to download packages from Bioconductor. To install package dependencies, use the pacman R package. Note that some of these packages may require additional library dependencies that need to be installed prior to their own successful installation.

```
install.packages("pacman")

library(pacman)

pacman::p_load(igraph, RColorBrewer, qs, rjson, RSQLite)
```

**Step 2. Install the package**

OmicsNetR is freely available from GitHub. The package documentation, including the vignettes for each module and user manual is available within the downloaded R package file. If all package dependencies were installed, you will be able to install the OmicsNetR. 

Install the package directly from github using the _devtools_ package. Open R and enter:

```

# Step 1: Install devtools
install.packages(devtools)
library(devtools)

# Step 2: Install OmicsNetR WITHOUT documentation
devtools::install_github("xia-lab/OmicsNetR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"))

# Step 2: Install OmicsNetR WITH documentation
devtools::install_github("xia-lab/OmicsNetR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
```

# Tips for using the OmicsNetR package

1. The first function that you will use in every module is the `Init.Data` function, which initiates the _dataSet_ object that stores user's data for further processing and analysis.
2. The OmicsNetR package will output data files/tables/analysis/networks outputs in your current working directory.
3. Every function must be executed in sequence as it is shown on the R Command history, please do not skip any commands as this can result in errors downstream.
4. Each main function in OmicsNetR is documented. Use the _?Function_ format to open its documentation. For instance, use `?OmicsNetR::QueryNet` to find out more about this function.

# Examples

## Starting from a list of genes

```
library(OmicsNetR)

# Step 1. Initiate the dataSet object
dataSet<-Init.Data()

# Step 2. Map list of genes to the application
dataSet<-PrepareInputList(dataSet,"#Entrez  logFC
4495  61.12
4496  51.06
4499  23.79
6354  21.04
6369  19.76", "hsa", "gene", "entrez", "direct");

# Step 3. Build a PPI interaction network
dataSet<-QueryNet(dataSet, "gene", "innate")

# Step 4. Build a PPI interaction network
CreateGraph();

# Step 5. Prepare the network file, the output will be in JSON format.
dataSet<-PrepareNetwork(dataSet, "subnetwork1", "omicsnet_0.json")
```
