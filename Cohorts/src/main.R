rm(list = ls())
main.dir = "my_directory/Cohorts/"
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
### NHS, NHS2 and HPFS data can be requested from the study PI 
################################################################################

# example dataset provided here
load("./data/mets_eset_with_covars_example.RData")
# example phenotype data
pdata = pData(mets.glauc)
head(pdata)
# example metabolite annotation data
fdata = fData(mets.glauc)
head(fdata)
# example metabolite values
mdata = exprs(mets.glauc)
mdata[1:10,1:5]

### create the dataset
source("./src/createStratifiedDatasets.R")
source("./src/createResiduals.R")

# Table 1: ./results/table1_nhs.csv
source("./src/createTable1.R")

### calculate number of effective tests
# NEF = 132 based on the real data. This NEF value will be used in this analysis
# the following file will calculate NEF based on the example dataset
source("./src/runNEF.R")

### run main analysis
# Figure 1: ./results/heatmap_sig_metabolites_main_all_models.png
# Supplementary Table 1 and data for Figure 1: ./results/main_analysis_all_models.csv
source("./src/runMainAnalysisOnAll.R")

# Figure 2: ./results/heatmap_GSEA_broad_main_all_models.png & 
# Data for Figure 2: ./results/fgsea_main_all_models.csv
source("./src/runGSEAOnAll.R")

### analysis restricted to those with POAG with paracentral visual field loss 
source("./src/runMainAnalysisOnGlaucPa.R") 
# Figure 3a: heatmap_GSEA_broad_main_poag_pa_models.png
# Data for Figure 3a: ./results/fgsea_main_poag_pa_models.csv
source("./src/runGSEAOnPOAGPa.R")

### analysis restricted to those with POAG with peripheral visual field loss
source("./src/runMainAnalysisOnGlaucPe.R")
# Figure 3b: heatmap_GSEA_broad_main_poag_pe_models.png
# Data for Figure 3b: ./results/fgsea_main_poag_pe_models.csv
source("./src/runGSEAOnPOAGPe.R")

# This plot only works for the real data 
# because metabolite names need to include number of double bonds and Carbon atoms
# The exmaple dataset uses simpler metabolite names: Metabolite 1, Metabolite 2, etc.
# Supplementary Figure S1: ./results/POAG.png
# source("./src/createBubblePlots.R")


### stratified analyses
# Supplementary Figure S2: ./results/heatmap_sig_metabolites_by_Age_all_models.png
source("./src/runMainAnalysisByAge.R")

# Supplementary Figure S3: ./results/heatmap_sig_metabolites_by_sex_all_models.png
source("./src/runMainAnalysisBySex.R")

# Supplementary Figure S4: ./results/heatmap_sig_metabolites_by_TimeDX_all_models.png
source("./src/runMainAnalysisByTimeBloodDX.R")

# Supplementary Figure S5: ./results/heatmap_sig_metabolites_by_FamhHist_all_models.png
source("./src/runMainAnalysisByFamHist.R")

# Supplementary Figure S6: ./results/heatmap_sig_metabolites_by_BMI_all_models.png
source("./src/runMainAnalysisByBMI.R")
