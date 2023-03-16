rm(list = setdiff("main.dir", ls()))

library(Biobase)

source("./src/myRegressions.R")

print("Running analyses stratified by BMI ...")

load("data/mets_eset_with_covars_for_stratified_analyses_and_residuals.RData")
load("data/mets_eset_with_covars_example.RData")

fdata = fData(mets.glauc)
colnames(fdata) = gsub(pattern = "class_broad",replacement = "super_class",x = colnames(fdata),fixed = T)
fdata$sub_class = fdata$super_class
colnames(fdata)[4] = "mean_cv"

pdata = pData(mets.glauc)

#### NHS 1 M1
nhs.clog.reg.res = runRegression(regression.type = "unconditional.logistic.regression",
                                 my.eset = normal.bmi.m1,
                                 outcome = "caco", 
                                 covariates = c("ageyr"),                                                
                                 strata = NULL, allow.na.values.in.exposure = T,
                                 file.with.result.prefix = file.path(main.dir, 
                                                                     "/results/log_nbmi_m1"),
                                 store.complete.result = T)

nhs.or.table = calculateORCI(result.data = nhs.clog.reg.res, or.scale = 1, file.with.result.prefix =
                               file.path(main.dir, "/results/log_nbmi_m1"),
                             estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = T, estimate.CI = F)

nhs.clog.reg.res = nhs.clog.reg.res[rownames(nhs.or.table),]
all.equal(nhs.clog.reg.res$METABOLITE, nhs.or.table$METABOLITE)
nhs.or.table$p.value = nhs.clog.reg.res$p.value

res1.m1 = merge(x = nhs.or.table, by.x = "METABOLITE", y = fdata[,c("METABOLITE","super_class","sub_class","missing","mean_cv")], by.y = "METABOLITE") 
res1.m1 = res1.m1[order(res1.m1$p.value),]
res1.m1$FDR = p.adjust(p = res1.m1$p.value,method = "fdr")
res1.m1$NEF = ifelse(res1.m1$p.value*132>1,1,res1.m1$p.value*132) # NEF = 132 see file "./runNEFF.R"

#### HPFS M1
hpfs.clog.reg.res = runRegression(regression.type = "unconditional.logistic.regression",
                                  my.eset = high.bmi.m1,
                                  outcome = "caco", 
                                  covariates = c("ageyr"),
                                  strata = NULL, allow.na.values.in.exposure = T,
                                  file.with.result.prefix = file.path(main.dir, 
                                                                      "/results/log_hbmi_m1"),
                                  store.complete.result = T)

hpfs.or.table = calculateORCI(result.data = hpfs.clog.reg.res, or.scale = 1, file.with.result.prefix =
                               file.path(main.dir, "/results/log_hbmi_m1"),
                             estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = T, estimate.CI = F)

hpfs.clog.reg.res = hpfs.clog.reg.res[rownames(hpfs.or.table),]
all.equal(hpfs.clog.reg.res$METABOLITE, hpfs.or.table$METABOLITE)
hpfs.or.table$p.value = hpfs.clog.reg.res$p.value

res2.m1 = merge(x = hpfs.or.table, by.x = "METABOLITE", y = fdata[,c("METABOLITE","super_class","sub_class","missing","mean_cv")], by.y = "METABOLITE") 
res2.m1 = res2.m1[order(res2.m1$p.value),]
res2.m1$FDR = p.adjust(p = res2.m1$p.value,method = "fdr")
res2.m1$NEF = ifelse(res2.m1$p.value*132>1,1,res2.m1$p.value*132) 

# overall result
load("./results/image_runMainAnalysisOnAll.RData")

m1.res.h1 = merge(x = res.m1[,c(1,2,3,4,8,9)], by.x = "METABOLITE", y = res1.m1[,c(1,3,4,10,11)], by.y = "METABOLITE", suffixes = c(".m1.all",""))
m1.res = merge(x = m1.res.h1, by.x = "METABOLITE", y = res2.m1[,c(1,3,4,10,11)], by.y = "METABOLITE", suffixes = c(".m1.NHS",".m1.HPFS"))

fdata = fData(mets.glauc)
colnames(fdata) = gsub(pattern = "class_broad",replacement = "super_class",x = colnames(fdata),fixed = T)
fdata$sub_class = fdata$super_class
colnames(fdata)[4] = "mean_cv"


#### NHS 1 M4 
nhs.clog.reg.res = runRegression(regression.type = "unconditional.logistic.regression",
                                 my.eset = normal.bmi.m5,
                                 outcome = "caco", 
                                 covariates = c("ageyr",
                                               "ibmils", "acts", 
                                               "nses" , "agemeno", # "nsesm"
                                               "nitras", "alcos", "caffs", "aheis" ,"calors"),                                               
                                 strata = NULL, allow.na.values.in.exposure = T,
                                 file.with.result.prefix = file.path(main.dir, 
                                                                     "/results/log_nbmi_m5"),
                                 store.complete.result = T)

nhs.or.table = calculateORCI(result.data = nhs.clog.reg.res, or.scale = 1, file.with.result.prefix =
                               file.path(main.dir, "/results/log_nbmi_m5"),
                             estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = T, estimate.CI = T)

nhs.clog.reg.res = nhs.clog.reg.res[rownames(nhs.or.table),]
all.equal(nhs.clog.reg.res$METABOLITE, nhs.or.table$METABOLITE)
nhs.or.table$p.value = nhs.clog.reg.res$p.value

res1.m4 = merge(x = nhs.or.table, by.x = "METABOLITE", y = fdata[,c("METABOLITE","super_class","sub_class","missing","mean_cv")], by.y = "METABOLITE") 
res1.m4 = res1.m4[order(res1.m4$p.value),]
res1.m4$FDR = p.adjust(p = res1.m4$p.value,method = "fdr")
res1.m4$NEF = ifelse(res1.m4$p.value*132>1,1,res1.m4$p.value*132) # NEF = 132 see file "./runNEFF.R"

hpfs.clog.reg.res = runRegression(regression.type = "unconditional.logistic.regression",
                                  my.eset = high.bmi.m5,
                                  outcome = "caco", 
                                  covariates =  c("ageyr",
                                                  "ibmils", "acts", 
                                                  "nses", "agemeno",# "nsesm","hrace"
                                                  "nitras", "alcos", "caffs", "aheis" ,"calors"),                                                  
                                  strata = NULL, allow.na.values.in.exposure = T,
                                  file.with.result.prefix = file.path(main.dir, 
                                                                      "/results/log_hbmi_m5"),
                                  store.complete.result = T)

hpfs.or.table = calculateORCI(result.data = hpfs.clog.reg.res, or.scale = 1, file.with.result.prefix =
                                file.path(main.dir, "/results/log_hbmi_m5"),
                              estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = T, estimate.CI = F)
hpfs.clog.reg.res = hpfs.clog.reg.res[rownames(hpfs.or.table),]
all.equal(hpfs.clog.reg.res$METABOLITE, hpfs.or.table$METABOLITE)
hpfs.or.table$p.value = hpfs.clog.reg.res$p.value

res2.m4 = merge(x = hpfs.or.table, by.x = "METABOLITE", y = fdata[,c("METABOLITE","super_class","sub_class","missing","mean_cv")], by.y = "METABOLITE") 
res2.m4 = res2.m4[order(res2.m4$p.value),]
res2.m4$FDR = p.adjust(p = res2.m4$p.value,method = "fdr")
res2.m4$NEF = ifelse(res2.m4$p.value*132>1,1,res2.m4$p.value*132)

save.image("./results/image_runMainAnalysisBMI.RData")
##############################################################
library(ggplot2)

res.m5$strata = "Overall"
res.m5$model = "Model 5"

res.m1$strata = "Overall"
res.m1$model = "Model 1"

res1.m1$strata = "BMI<25"
res1.m1$model = "Model 1"

res2.m1$strata = "BMI>=25"
res2.m1$model = "Model 1"

res1.m4$strata = "BMI<25"
res1.m4$model = "Model 5"

res2.m4$strata = "BMI>=25"
res2.m4$model = "Model 5"

res.all = rbind(res.m1[,c(1,3,4,8:11)], res1.m1[,c(1,3,4,10,11:13)], res2.m1[,c(1,3,4,10,11:13)],
                res.m5[,c(1,3,4,8:11)], res1.m4[,c(1,3,4,10,11:13)], res2.m4[,c(1,3,4,10,11:13)])

res.all$OR = sapply(X = res.all$OR.CI,FUN = function(x){as.numeric(strsplit(x = x,split = " ")[[1]][1])})

met.sig = unique(c(res.m1$METABOLITE[which(res.m1$p.value<0.05)],res.m4$METABOLITE[which(res.m4$p.value<0.05)]))

res.sig = res.all[which(res.all$METABOLITE %in% met.sig),]
res.sig$strata = factor(res.sig$strata)
res.sig$strata = relevel(factor(res.sig$strata), ref = "Overall")

# plot 
res.sig$stars = ""
res.sig$stars[which(res.sig$p.value<0.05)] = "*"
res.sig$stars[which(res.sig$NEF<0.2)] = "**"
res.sig$stars[which(res.sig$NEF<0.05)] = "***"

res.sig$METABOLITE = factor(res.sig$METABOLITE, levels = met.sig)
range(res.sig$OR, na.rm =T)
library(RColorBrewer)
hm.palette <- rev(brewer.pal(11, "RdBu"))

ggplot(data = res.sig, aes(x = strata, y = METABOLITE, fill = OR)) + 
  geom_tile(color = "white", size = 0.1) + 
  geom_text(aes(y = METABOLITE, x = strata, label=stars), color="black", size=5,
           nudge_x = 0, nudge_y = -0.3, fontface = "bold") +
  facet_wrap(.~model)+
  scale_fill_gradientn(colours = hm.palette, limits = c(0.3,1.7), breaks = c(0.5, 1, 1.5), na.value = "black")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5,
                                   size = 9, face = "bold", colour = "black"),
        strip.text = element_text(face = "bold", colour = "black",size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", colour = "black"),
        axis.text.y = element_text(face = "bold", size = 9, color = "black"),
        legend.title = element_text(face = "bold", size = 11, color = "black"),
        legend.text = element_text(face = "bold", size = 10, color = "black"),
        plot.caption = element_text(hjust = 1,size = 10, face = "bold"))+
  labs(caption = "* p<0.05; ** NEF-p<0.2; *** NEF-p<0.05 ")

ggsave(filename=file.path("./results/heatmap_sig_metabolites_by_BMI_all_models.png"), 
       device = "png", dpi = 600, width = 10, height = 9)  

print("Figure: ./results/heatmap_sig_metabolites_by_BMI_all_models.png")

# results table
rownames(res.m1) = res.m1$METABOLITE
colnames(res.m1) = paste(colnames(res.m1),".all.m1", sep ="")
rownames(res.m5) = res.m5$METABOLITE
colnames(res.m5) = paste(colnames(res.m5),".all.m5", sep ="")

rownames(res1.m1) = res1.m1$METABOLITE
colnames(res1.m1) = paste(colnames(res1.m1),".normalBMI.m1", sep ="")
rownames(res1.m4) = res1.m4$METABOLITE
colnames(res1.m4) = paste(colnames(res1.m4),".normalBMI.m5", sep ="")

rownames(res2.m1) = res2.m1$METABOLITE
colnames(res2.m1) = paste(colnames(res2.m1),".highBMI.m1", sep ="")
rownames(res2.m4) = res2.m4$METABOLITE
colnames(res2.m4) = paste(colnames(res2.m4),".highBMI.m5", sep ="")


res.all = cbind(res.m1[rownames(res.m1),c(1,2,6,7,3,4,10)], res1.m1[rownames(res.m1),c(3,4,10)], res2.m1[rownames(res.m1),c(3,4,10)],
                res.m5[rownames(res.m1),c(1,2,6,7,3,4,10)], res1.m4[rownames(res.m1),c(3,4,10)], res2.m4[rownames(res.m1),c(3,4,10)])
write.table(res.all, file = "./results/analysis_by_BMI_all_models.csv", sep = "\t", row.names = F, col.names = T, quote = F)

print("Done")
