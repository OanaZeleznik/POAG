rm(list = setdiff("main.dir", ls()))
print("Creating Table 1 ...")

library(tableone)
library(Biobase)

load("./data/poag_ukbb_es_example.RData")

pdata.t1 = pData(poag.ukb)

pdata.t1$sex = factor(pdata.t1$sex)
levels(pdata.t1$sex) = c("female", "male")
table(pdata.t1$sex)

pdata.t1$smoking = factor(pdata.t1$smoking)
levels(pdata.t1$smoking) = c("non-smoker", "past smoker", "current smoker")
table(pdata.t1$smoking)

pdata.t1$alcohol = factor(pdata.t1$alcohol)
levels(pdata.t1$alcohol) = c("daily or almost daily", "3-4 times a week", "1-2 times a week",
                             "1-3 times a month","special occasions only or never or prefer not to answer ")
table(pdata.t1$alcohol)

pdata.t1$ethnicity = factor(pdata.t1$ethnicity)
levels(pdata.t1$ethnicity) = c("white", "asian","black", "other")
table(pdata.t1$ethnicity)

all.vars = c("age","sex","smoking","sumMETperWK_0_0","BMI","ethnicity", 
               "coffee_0_new6", "tea_0_new6","alcohol","sbpa_0_mean","diabetes", 
               "cad", "cholesterol")
               
               
cat.vars = c("sex","smoking","ethnicity", 
             "alcohol","diabetes", 
             "cad")

all.t1 = CreateTableOne(vars = all.vars, strata = "any_glaucoma_Q_cohort", data = pdata.t1, factorVars = cat.vars)
all.t1.mat = print(all.t1, catDigits = 1, contDigits = 1, pDigits = 1, test = TRUE)
write.csv(all.t1.mat, file = "./results/table1.csv", quote = F)

print("Done.")
