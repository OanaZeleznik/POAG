## Plasma metabolite profile for primary open-angle glaucoma in three US cohorts and the UK Biobank

This code was used to generate the results of our [publication](https://www.medrxiv.org/content/10.1101/2022.02.24.22271483v1). 

### REFERENCE
If you use code or parts of this code in your work, please cite our publication:

Zeleznik, O.A., Kang, J.H., Lasky-Su, J.A., Eliassen, A.H., Frueh, L., Clish, C., Rosner, B.A., Elze, T., Hysi, P., Khawaja, A.P. and Wiggs, J.L., 2022. Plasma metabolomics of primary open-angle glaucoma in three prospective US cohorts and the UK Biobank. [medRxiv](https://www.medrxiv.org/content/10.1101/2022.02.24.22271483v1).

### ABSTRACT

Glaucoma is a progressive optic neuropathy that is a leading cause of irreversible blindness worldwide. Primary open-angle glaucoma is the most common form, and yet the etiology of this multifactorial disease is poorly understood. We aimed to identify plasma metabolites associated with risk of developing glaucoma in a case-control study (599 cases and 599 matched controls) nested within the Nurses’ Health Studies and Health Professionals Follow-Up Study. Plasma metabolites were measured with liquid chromatography mass spectrometry (Broad Institute, USA); 367 metabolites from 17 metabolite classes passed quality control analyses. For comparison, in a cross-sectional study in the UK Biobank, 168 metabolites were measured in plasma samples from 2,238 prevalent glaucoma cases and 44,723 controls using nuclear magnetic resonance spectroscopy (Nightingale, Finland). Here, we show higher levels of diglycerides and triglycerides are adversely associated with glaucoma in all cohorts, suggesting that they play an important role in glaucoma pathogenesis.

## Setup

The code is divided into two parts: The folder “Cohorts” includes all code needed 
to run the analyses among the health professional cohorts. The folder “UKBiobank” 
includes all code needed to run the analyses on the UK Biobank data. Each folder 
is further divided into the subfolders: “data”,”src”, “results”, and “expected_results”. 
The data folder includes an example dataset that can be used to run the code. 
The R code is stored in the “src” folder. Results will be stored in the “results” 
folder and expected results are provided in the “expected_results” folder.

### System requirements

The code was tested using R (version 4.0.3: https://cran.rstudio.com/) and 
RStudio (https://posit.co/downloads/) in CentOS Linux 7 and Windows 11. For Windows, 
RTools are required and can be downloaded at https://cran.rstudio.com/bin/windows/Rtools/.

The following R packages are required and the code was tested with the versions as stated: 
Biobase version 2.50, fgsea version 1.16, tableone version 0.13,  ggplot2 version 3.3.5, 
RColorBrewer version 1.1-2, svglite 2.1.0, data.table 1.14.2, tidyverse version 1.3.1, 
haven version 2.4.3, labelled version 2.10.0, plyr 1.8.8, qvalue version 2.22.0, 
multtest version 2.46.0

### Installation guide

All R code is available in the “Cohorts/src” and “UKBiobank/src” folders. The whole analysis workflow is provided in 
src/main.R. The code can be run by following the instructions in the src/main.R file. 
Make sure to set the path to your study folder in main.R (main.dir = “my_folder/Cohorts/” 
or main.dir = “my_folder/UKBiobank/”)

Typical installation time, including all required packages, is around 1 hour on a typical 
desktop computer.

### Demo

As data can not be shared due to privacy protocols in place, the authors provide example 
datasets in the Cohorts/data/ (mets_eset_with_covars_example.RData) and UKBiobank/data 
(poag_ukbb_es_example.RData) folders. The code will run on any other data formatted in 
the same way as the example data. Details on the format of the example data are provided 
in /src/main.R. The expected output is provided in the Cohorts/expected_results and 
UKBiobank/expected_results folders.

The expected run time for demo on a typical desktop computer is around 30 min for each study. 

### Instructions for use

The code will run with any data stored in the same format as the example data provided here 
(stored in the “data” folders: mets_eset_with_covars_example.RData and poag_ukb_es_example.RData). 
Each dataset is of type ExpressionSet (as defined in the Biobase package) and includes 3 
data frames: phenotype data (can be extracted with the function pData), metabolite 
annotation (can be extracted with the function fData) and metabolite values (can be extracted 
with the function exprs()).

All R code is available in the “Cohorts/src” and “UKBiobank/src” folders. The whole analysis workflow 
is provided in src/main.R. The code can be run by following the instructions in the src/main.R file.

## Disclaimer

The code and data of this repository are provided to promote reproducible research. They are not intended for clinical care or commercial use.

The software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the software or the use or other dealings in the software.
