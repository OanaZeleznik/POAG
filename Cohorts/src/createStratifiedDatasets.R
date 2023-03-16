rm(list = setdiff("main.dir", ls()))

library(Biobase)
library(Hmisc)

print("Create stratified datasets ...")

load("./data/mets_eset_with_covars_example.RData")

mets.glauc.pdata = pData(mets.glauc)

mets.glauclo = mets.glauc[,-1*which(mets.glauc$glaulo == "0")]
pdata = pData(mets.glauclo)
table(pdata$glaulo)

mets.glauchi = mets.glauc[,-1*which(mets.glauc$glauhi == "0")]
pdata = pData(mets.glauchi)
table(pdata$glauhi)

mets.glaucpa = mets.glauc[,-1*which(mets.glauc$glaupa == "0")]
pdata = pData(mets.glaucpa)
table(pdata$glaupa)

summary(pdata$ibmils)
mets.glaucpa.normal.bmi = mets.glaucpa[,which(pdata$ibmils<25)]
pdata = pData(mets.glaucpa.normal.bmi)
table(pdata$glaupa)

pdata = pData(mets.glaucpa)
mets.glaucpa.high.bmi = mets.glaucpa[,which(pdata$ibmils>=25)]
pdata = pData(mets.glaucpa.high.bmi)
table(pdata$glaupa)


mets.glaucpe = mets.glauc[,-1*which(mets.glauc$glaupe == "0")]
pdata = pData(mets.glaucpe)
table(pdata$glaupe)

save(mets.glauc, mets.glauchi, mets.glauclo, mets.glaucpa, 
     mets.glaucpa.normal.bmi,mets.glaucpa.high.bmi, mets.glaucpe, 
     file = "data/mets_eset_with_covars_by_subtype.RData")

# create additional stratified datasets
mets.glauc.pdata = pData(mets.glauc)

# mean years between blood collection and DX by sex
tapply(mets.glauc.pdata$blood_to_dx, mets.glauc.pdata$sex, mean, na.rm=T)


# mean age at DX, by sex
mets.glauc.pdata$ageDX = mets.glauc.pdata$ageyr+mets.glauc.pdata$blood_to_dx
summary(mets.glauc.pdata$ageDX)

# age at blood draw, DX and time between these
mean(mets.glauc.pdata$ageyr[which(mets.glauc.pdata$caco == "1")],na.rm=T) 
mean(mets.glauc.pdata$ageDX[which(mets.glauc.pdata$caco == "1")],na.rm=T) 
mean(mets.glauc.pdata$blood_to_dx[which(mets.glauc.pdata$caco == "1")],na.rm=T) 

sd(mets.glauc.pdata$ageyr[which(mets.glauc.pdata$caco == "1")],na.rm=T) 
sd(mets.glauc.pdata$ageDX[which(mets.glauc.pdata$caco == "1")],na.rm=T) 
sd(mets.glauc.pdata$blood_to_dx[which(mets.glauc.pdata$caco == "1")],na.rm=T) 

# by cohort 
table(mets.glauc.pdata$cohort)

nhs.glauc = mets.glauc[,which(mets.glauc.pdata$cohort == "NHS")]
dim(nhs.glauc)

nhs2.glauc = mets.glauc[,which(mets.glauc.pdata$cohort == "NHS2")]
dim(nhs2.glauc)

hpfs.glauc = mets.glauc[,which(mets.glauc.pdata$cohort == "HPFS")]
dim(hpfs.glauc)

# by sex
table(mets.glauc.pdata$sex)
women.glauc = mets.glauc[,which(mets.glauc.pdata$sex == "F")]
dim(women.glauc)

men.glauc = mets.glauc[,which(mets.glauc.pdata$sex == "M")]
dim(men.glauc)

# by age
summary(mets.glauc.pdata$ageyr)
age.low.glauc = mets.glauc[,which(mets.glauc.pdata$ageyr<58.5)]
dim(age.low.glauc)

age.high.glauc = mets.glauc[,which(mets.glauc.pdata$ageyr>=58.5)]
dim(age.high.glauc)

# by family history
table(mets.glauc.pdata$hgla)
fh.glauc = mets.glauc[,which(mets.glauc.pdata$hgla == " : yes")]
dim(fh.glauc)

nfh.glauc = mets.glauc[,which(mets.glauc.pdata$hgla == " : no")]
dim(nfh.glauc)

# by BMI
summary(mets.glauc.pdata$ibmils)
normal.bmi = mets.glauc[,which(mets.glauc.pdata$ibmils<25)]
dim(normal.bmi)
pdata = pData(normal.bmi)
table(pdata$caco)

high.bmi = mets.glauc[,which(mets.glauc.pdata$ibmils>=25)]
dim(high.bmi)
pdata = pData(high.bmi)
table(pdata$caco)

# by BMI=22
summary(mets.glauc.pdata$ibmils)
normal.bmi22 = mets.glauc[,which(mets.glauc.pdata$ibmils<22)]
dim(normal.bmi22)
pdata = pData(normal.bmi22)
table(pdata$caco)

high.bmi22 = mets.glauc[,which(mets.glauc.pdata$ibmils>=22)]
dim(high.bmi22)
pdata = pData(high.bmi22)
table(pdata$caco)

# by time between sample collection and DX
summary(mets.glauc.pdata$blood_to_dx[which(mets.glauc.pdata$caco == "1")])
closeDX.glauc = mets.glauc[,which(mets.glauc.pdata$blood_to_dx<10.17)]
dim(closeDX.glauc)

farDX.glauc = mets.glauc[,which(mets.glauc.pdata$blood_to_dx>=10.17)]
dim(farDX.glauc)

save(nhs.glauc, nhs2.glauc, hpfs.glauc,
     women.glauc, men.glauc,
     age.low.glauc, age.high.glauc,
     fh.glauc, nfh.glauc,
     normal.bmi, high.bmi,
     normal.bmi22, high.bmi22,
     closeDX.glauc, farDX.glauc,
     file="data/mets_eset_with_covars_for_stratified_analyses.RData")

### created additional stratified datasets
mets.glauc.pdata = pData(mets.glaucpa)
dim(mets.glaucpa)

# by BMI
summary(mets.glauc.pdata$ibmils)
normal.bmi.glauc.pa = mets.glaucpa[,which(mets.glauc.pdata$ibmils<25)]
dim(normal.bmi.glauc.pa)
pdata = pData(normal.bmi.glauc.pa)
table(pdata$caco)

high.bmi.glauc.pa = mets.glaucpa[,which(mets.glauc.pdata$ibmils>=25)]
dim(high.bmi.glauc.pa)
pdata = pData(high.bmi.glauc.pa)
table(pdata$caco)

save(normal.bmi.glauc.pa, high.bmi.glauc.pa, file = "data/mets_eset_with_covars_glauc_paracentral_byBMI.RData")
print("Done")
