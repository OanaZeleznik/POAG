rm(list = setdiff("main.dir", ls()))
print("Creating stratified datasets ...")

library(Biobase)
load(file = "./data/poag_ukbb_es_example.RData")

### create stratified datasets
# by age
pdata = pData(poag.ukb)
summary(pdata$age)

poag.young = poag.ukb[,which(pdata$age <= 58)]
dim(poag.young)

poag.old = poag.ukb[,which(pdata$age > 58)]
dim(poag.old)

# by BMI 
summary(pdata$BMI)
poag.lean = poag.ukb[,which(pdata$BMI <=25)]
dim(poag.lean)

poag.not.lean = poag.ukb[,which(pdata$BMI >25)]
dim(poag.not.lean)

# by sex
table(pdata$sex)

poag.women = poag.ukb[,which(pdata$sex == 0)]
dim(poag.women)

poag.men = poag.ukb[,which(pdata$sex == 1)]
dim(poag.men)

save(poag.ukb, poag.women, poag.men, poag.young, poag.old, poag.lean, poag.not.lean, file = "./data/poag_ukbb_es_strata.RData")

print("Done.")