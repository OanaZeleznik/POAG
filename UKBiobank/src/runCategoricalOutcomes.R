rm(list = setdiff("main.dir", ls()))
print("Running multivariable unconditional logistic regression for each metabolite - this can take up to 15 min, depending on your system ...")

# functions to run logistic regressions
source("./src/myRegressions.R")

library(Biobase)
load("./data/poag_ukbb_es_strata.RData")
load("./results/ukb_nef.RData")


addFDRandNEF = function(raw.res, or.res, nef, fdata){
  
  raw.res = raw.res[rownames(or.res),]
  
  if(!all.equal(raw.res$METABOLITE, or.res$METABOLITE))
    print("STOP - order does not match!")
  or.res$p.value = raw.res$p.value
  
  res.m1 = merge(x = or.res, by.x = "HMDB_ID", y = fdata[,c("HMDB_ID","is.indicator")], by.y = "HMDB_ID") 
  res.m1 = res.m1[order(res.m1$p.value),]
  res.m1$FDR = p.adjust(p = res.m1$p.value,method = "fdr")
  res.m1$NEF = ifelse(res.m1$p.value*nef>1,1,res.m1$p.value*nef) # nef = 133 see file "./QCAnalysis.R"
  
  print(paste("There are ", length(which(res.m1$p.value<0.05)),"mets with p<0.05"))
  print(paste("There are ", length(which(res.m1$FDR<0.2)),"mets with FDR<0.2"))
  print(paste("There are ", length(which(res.m1$NEF<0.2)),"mets with NEF<0.2"))
  
  return(res.m1)
}

fdata = fData(poag.ukb)
mdata = exprs(poag.ukb)
pdata = pData(poag.ukb)

uclog.reg.res.m2 = runRegression(regression.type = "unconditional.logistic.regression",
                                 my.eset = poag.ukb,
                                 outcome = "any_glaucoma_Q_cohort", 
                                 covariates = c("age","age.squared","sex","deprivation", 
                                                "SE.y", "smoking", "no_cigarettes_current_0_0",  
                                                "sumMETperWK_0_0", "BMI", "ethnicity", 
                                                "coffee_0_new6", "tea_0_new6",  "alcohol",
                                                "sbpa_0_mean", "diabetes", "cad","bblocker2", "statins"), 
                                 strata = NULL, allow.na.values.in.exposure = T,
                                 file.with.result.prefix = file.path(main.dir, 
                                                                     "/results/uclog_new"),
                                 store.complete.result = F)

or.table.m2 = calculateORCI(result.data = uclog.reg.res.m2, or.scale = 1, file.with.result.prefix =
                              file.path(main.dir, "/results/uclog_new"),
                            estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = T, estimate.CI = T)

res.m2 = addFDRandNEF(raw.res = uclog.reg.res.m2, or.res = or.table.m2, nef = ukb.nef, fdata = fdata)
save(uclog.reg.res.m2, res.m2, file = "results/uclog_all_new.RData")

all.res = res.m2[,c(1:4,8)]
all.res$Classification = fdata$Group[match(x = res.m2$HMDB_ID, table = fdata$HMDB_ID)]
all.res$Subclassification = fdata$Subgroup[match(x = res.m2$HMDB_ID, table = fdata$HMDB_ID)]
all.res$Unit = "mmol/l"

all.res = all.res[,c(2,8,6,7,3,4,5)]
write.table(x = all.res, file = "./results/ukbb_all_results.csv",quote = F, row.names = F, col.names = T, sep = "\t")

print("File with all results: ./results/ukbb_all_results.csv")

print("Done")
