rm(list = setdiff("main.dir", ls()))

library(Biobase)
library(ggplot2)
library(tableone)

print("creating Table 1 ...")

load("data/mets_eset_with_covars_example.RData")

pdata = pData(mets.glauc)
pdata$bldtime = factor(pdata$bldtime)
levels(pdata$bldtime) = c(rep("12:00am-7:59am",4),"8:00am-9:59am","10:00am-11:59am", rep("12:00pm-11:59pm",6))

pdata$menobld = as.character(pdata$menobld)
pdata$menodiag = as.character(pdata$menodiag)

pdata$menobld[which(pdata$menobld == "missing" & pdata$sex == "M")] = "men"
pdata$menodiag[which(pdata$menodiag == "missing" & pdata$sex == "M")] = "men"

all.vars = c("ageyr", "sex", "hrace", "bldtime", "bldmon",
              "smokec", "ibmils", "acts",  
              "hgla", "agemeno", "nses", 
              "nitras", "alcos", "caffs", "aheis" ,"calors",
              "hbph","cholh","sdbh")
cat.vars = c("sex", "hrace", "bldtime", "bldmon",
             "smokec", 
             "hbph","cholh","sdbh")

all.t1 = CreateTableOne(vars = all.vars, strata = "caco", data = pdata, factorVars = cat.vars)
all.t1.mat = print(all.t1, catDigits = 1, contDigits = 1, pDigits = 1, test = TRUE)
write.csv(all.t1.mat, file = file.path(main.dir,"results/table1_nhs.csv"), quote = F)
print("Done")