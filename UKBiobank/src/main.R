rm(list = ls())
main.dir = "my_folder/UKBiobank/"
setwd(main.dir)

# the following packages are required:
require(Biobase)
require(fgsea)
require(tableone)
require(ggplot2)
require(RColorBrewer)
require(svglite)
require(data.table)
require(tidyverse)
require(haven)

require(Hmisc)
require(reshape)
require(labelled)
require(plyr)
require(qvalue)
require(multtest)

################################################################################
### This code can be run with the provided example dataset 
### UKBiobank data can be requested here: 
### https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access 
################################################################################

# example dataset provided here
load("./data/poag_ukbb_es_example.RData")
# example phenotype data
pdata = pData(poag.ukb)
head(pdata)
# example metabolite annotation data
fdata = fData(poag.ukb)
head(fdata)
# example metabolite values
mdata = exprs(poag.ukb)
mdata[1:10,1:5]

#### create the dataset
source("./src/createStratifiedDatasets.R")

#### Table 1
source("./src/createTable1.R")

### calculate number of effective tests
source("./src/calculateNEF.R")

### run main analysis
source("./src/runCategoricalOutcomes.R")

# Figure 4: ./results/heatmap_poag_new.png
source("./src/plotCategoricalResults.R")
# Figure 5: ./results/heatmap_poag_MSEA.png
source("./src/runCategoricalOutcomes_MSEA.R")
