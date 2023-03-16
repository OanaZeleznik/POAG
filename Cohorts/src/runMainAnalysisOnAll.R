rm(list = setdiff("main.dir", ls()))

library(Biobase)
library(ggplot2)

print("Running individual metabolites analyses - this will take a few minutes ...")

source("./src/myRegressions.R")
load("./data/mets_eset_with_covars_and_residuals.RData")

fdata = fData(mets.glauc)
pdata = pData(mets.glauc)

# categorical variables are adjusted using residuals

# m1 - CLR without adjustments
clog.reg.res.m1 = runRegression(regression.type = "conditional.logistic.regression",
                             my.eset = mets.glauc.m1,
                             outcome = "caco", 
                             covariates = NULL, 
                             strata = "matchid", allow.na.values.in.exposure = T,
                             file.with.result.prefix = file.path(main.dir, 
                                                                 "/results/clog_all_m1"),
                             store.complete.result = T)

or.table.m1 = calculateORCI(result.data = clog.reg.res.m1, or.scale = 1, file.with.result.prefix =
                                file.path(main.dir, "/results/clog_all_m1"),
                              estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = F, estimate.CI = F)

clog.reg.res.m1 = clog.reg.res.m1[rownames(or.table.m1),]
all.equal(clog.reg.res.m1$METABOLITE, or.table.m1$METABOLITE)
or.table.m1$p.value = clog.reg.res.m1$p.value

res.m1 = merge(x = or.table.m1, by.x = "METABOLITE", y = fdata[,c("METABOLITE","missing","meanCV")], by.y = "METABOLITE") 
res.m1 = res.m1[order(res.m1$p.value),]
res.m1$FDR = p.adjust(p = res.m1$p.value,method = "fdr")
res.m1$NEF = ifelse(res.m1$p.value*132>1,1,res.m1$p.value*132) # NEF = 132 see file "./runNEFF.R"

res.m1[which(res.m1$FDR<0.2),]

# model 2: factors that influence metabolite levels
clog.reg.res.m2 = runRegression(regression.type = "conditional.logistic.regression",
                             my.eset = mets.glauc.m2,
                             outcome = "caco", 
                             covariates = c("ageyr","ibmils", "acts"), 
                             strata = "matchid", allow.na.values.in.exposure = T,
                             file.with.result.prefix = file.path(main.dir, 
                                                                 "/results/clog_all_m2"),
                             store.complete.result = T)

or.table.m2 = calculateORCI(result.data = clog.reg.res.m2, or.scale = 1, file.with.result.prefix =
                           file.path(main.dir, "/results/clog_all_m2"),
                           estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = F, estimate.CI = F)

clog.reg.res.m2 = clog.reg.res.m2[rownames(or.table.m2),]
all.equal(clog.reg.res.m2$METABOLITE, or.table.m2$METABOLITE)
or.table.m2$p.value = clog.reg.res.m2$p.value

res.m2 = merge(x = or.table.m2, by.x = "METABOLITE", y = fdata[,c("METABOLITE","missing","meanCV")], by.y = "METABOLITE") 
res.m2 = res.m2[order(res.m2$p.value),]
res.m2$FDR = p.adjust(p = res.m2$p.value,method = "fdr")
res.m2$NEF = ifelse(res.m2$p.value*132>1,1,res.m2$p.value*132) # NEF = 132 see file "./runNEFF.R"

res.m2[which(res.m2$FDR<0.2),]

# model 3: set 1 of established risk factors for glaucoma
clog.reg.res.m3 = runRegression(regression.type = "conditional.logistic.regression",
                                my.eset = mets.glauc.m3,
                                outcome = "caco", 
                                covariates = c("ageyr", "ibmils", "acts", "nses" , "agemeno"),
                                strata = "matchid", allow.na.values.in.exposure = T,
                                file.with.result.prefix = file.path(main.dir, 
                                                                    "/results/clog_all_m3"),
                                store.complete.result = T)

or.table.m3 = calculateORCI(result.data = clog.reg.res.m3, or.scale = 1, file.with.result.prefix =
                              file.path(main.dir, "/results/clog_all_m3"),
                            estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = F, estimate.CI = F)

clog.reg.res.m3 = clog.reg.res.m3[rownames(or.table.m3),]
all.equal(clog.reg.res.m3$METABOLITE, or.table.m3$METABOLITE)
or.table.m3$p.value = clog.reg.res.m3$p.value

res.m3 = merge(x = or.table.m3, by.x = "METABOLITE", y = fdata[,c("METABOLITE","missing","meanCV")], by.y = "METABOLITE") 
res.m3 = res.m3[order(res.m3$p.value),]
res.m3$FDR = p.adjust(p = res.m3$p.value,method = "fdr")
res.m3$NEF = ifelse(res.m3$p.value*132>1,1,res.m3$p.value*132) # NEF = 132 see file "./runNEFF.R"

res.m3[which(res.m3$p.value<0.2),]

# model 4: set 2 of established risk factors for glaucoma
clog.reg.res.m4 = runRegression(regression.type = "conditional.logistic.regression",
                                my.eset = mets.glauc.m4,
                                outcome = "caco", 
                                covariates = c("ageyr", "ibmils", "acts", "nses" , "agemeno", 
                                               "nitras", "alcos", "caffs", "aheis" ,"calors"),
                                strata = "matchid", allow.na.values.in.exposure = T,
                                file.with.result.prefix = file.path(main.dir, 
                                                                    "/results/clog_all_m4"),
                                store.complete.result = T)

or.table.m4 = calculateORCI(result.data = clog.reg.res.m4, or.scale = 1, file.with.result.prefix =
                              file.path(main.dir, "/results/clog_all_m4"),
                            estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = F, estimate.CI = F)

clog.reg.res.m4 = clog.reg.res.m4[rownames(or.table.m4),]
all.equal(clog.reg.res.m4$METABOLITE, or.table.m4$METABOLITE)
or.table.m4$p.value = clog.reg.res.m4$p.value

res.m4 = merge(x = or.table.m4, by.x = "METABOLITE", y = fdata[,c("METABOLITE","missing","meanCV")], by.y = "METABOLITE") 
res.m4 = res.m4[order(res.m4$p.value),]
res.m4$FDR = p.adjust(p = res.m4$p.value,method = "fdr")
res.m4$NEF = ifelse(res.m4$p.value*132>1,1,res.m4$p.value*132) # NEF = 132 see file "./runNEFF.R"

res.m4[which(res.m4$p.value<0.05),]

# model 5 : co-morbidities
clog.reg.res.m5 = runRegression(regression.type = "conditional.logistic.regression",
                                my.eset = mets.glauc.m5,
                                outcome = "caco", 
                                covariates = c("ageyr", "menobld", "menodiag", "ibmils", "acts", 
                                               "nses" , "agemeno", "nitras", "alcos", "caffs", "aheis" ,"calors"),
                                strata = "matchid", allow.na.values.in.exposure = T,
                                file.with.result.prefix = file.path(main.dir, 
                                                                    "/results/clog_all_m5"),
                                store.complete.result = T)

or.table.m5= calculateORCI(result.data = clog.reg.res.m5, or.scale = 1, file.with.result.prefix =
                              file.path(main.dir, "/results/clog_all_m5"),
                            estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = F, estimate.CI = F)

clog.reg.res.m5 = clog.reg.res.m5[rownames(or.table.m5),]
all.equal(clog.reg.res.m5$METABOLITE, or.table.m5$METABOLITE)
or.table.m5$p.value = clog.reg.res.m5$p.value

res.m5 = merge(x = or.table.m5, by.x = "METABOLITE", y = fdata[,c("METABOLITE","missing","meanCV")], by.y = "METABOLITE") 
res.m5 = res.m5[order(res.m5$p.value),]
res.m5$FDR = p.adjust(p = res.m5$p.value,method = "fdr")
res.m5$NEF = ifelse(res.m5$p.value*132>1,1,res.m5$p.value*132) # NEF = 132 see file "./runNEFF.R"

res.m5[which(res.m5$p.value<0.05),]
save.image(file = "./results/image_runMainAnalysisOnAll.RData")
# ################################################################################
library(ggplot2)
# load("./results/image_runMainAnalysisOnAll.RData")

res.all.h1 = merge(x = res.m1, by.x = "METABOLITE", y = res.m2[,c(-2,-6,-7)], by.y = "METABOLITE", suffixes = c(".m1",""))
res.all.h23 = merge(x = res.all.h1, by.x = "METABOLITE", y = res.m3[,c(-2,-6,-7)], by.y = "METABOLITE", suffixes = c(".m2",""))
res.all.h3 = merge(x = res.all.h23, by.x = "METABOLITE", y = res.m4[,c(-2,-6,-7)], by.y = "METABOLITE", suffixes = c(".m3",""))
res.all.h4 = merge(x = res.all.h3, by.x = "METABOLITE", y = res.m5[,c(-2,-6,-7)], by.y = "METABOLITE", suffixes = c(".m4",".m5"))

res.all = res.all.h4[,c(1,2,7, 3,4,9, 10,11,14, 15,16,19, 20,21,24, 25,26,29)]
res.all = res.all[order(res.all$METABOLITE),]
save(res.all, file ="./results/main_analysis_all_models.RData")
write.table(res.all, file = "./results/main_analysis_all_models.csv", sep = "\t", row.names = F, col.names = T, quote = F)

res.all.sig = res.all[which(res.all$NEF.m1<0.2 | res.all$NEF.m2<0.2 | res.all$NEF.m3<0.2 | res.all$NEF.m4<0.2 | res.all$NEF.m5<0.2),]
# no metabolites are significant at NEF<0.2 in the example dataset:
res.all.sig = res.all[which(res.all$p.value.m1<0.01 | res.all$p.value.m1<0.01 | res.all$p.value.m1<0.01 | res.all$p.value.m1<0.01 | res.all$p.value.m1<0.01),]

write.table(res.all.sig, file = "./results/main_analysis_all_models_sig_mets.csv", sep = "\t", row.names = F, col.names = T, quote = F)

# nominal p-value
data.padj = res.all.sig[,c("METABOLITE",colnames(res.all.sig)[grep(pattern = "p.value",x = colnames(res.all.sig),fixed = T)])]

# NEF
data.nef= res.all.sig[,c("METABOLITE",colnames(res.all.sig)[grep(pattern = "NEF",x = colnames(res.all.sig),fixed = T)])]

data.nes =  res.all.sig[,c("METABOLITE",colnames(res.all.sig)[grep(pattern = "OR.CI",x = colnames(res.all.sig),fixed = T)])]
data.nes[,-1] = apply(data.nes[,-1],MARGIN = c(1,2), function(x) {return(as.numeric(strsplit(x = as.character(x), split = " ",fixed = T)[[1]][1]))})
data.help =  res.all.sig[,c("METABOLITE",colnames(res.all.sig)[grep(pattern = "OR.CI",x = colnames(res.all.sig),fixed = T)])]
data.help[,-1] = apply(data.help[,-1],MARGIN = c(1,2), function(x) {return(strsplit(x = as.character(x), split = " ",fixed = T)[[1]][2])})
data.help[,-1] = apply(data.help[,-1],MARGIN = c(1,2), function(x) {return(gsub(pattern = "(", replacement = "",x = x, fixed = T))})
data.help[,-1] = apply(data.help[,-1],MARGIN = c(1,2), function(x) {return(gsub(pattern = ")", replacement = "",x = x, fixed = T))})

data.lci = data.help
data.lci[,-1] = apply(data.lci[,-1],MARGIN = c(1,2), function(x) {return(as.numeric(strsplit(x = as.character(x), split = "-",fixed = T)[[1]][1]))})

data.uci = data.help
data.uci[,-1] = apply(data.uci[,-1],MARGIN = c(1,2), function(x) {return(as.numeric(strsplit(x = as.character(x), split = "-",fixed = T)[[1]][2]))})

colnames(data.nes) = c("Metabolite","Model 1","Model 2","Model 3","Model 4","Model 5")
colnames(data.padj) = c("Metabolite","Model 1","Model 2","Model 3","Model 4","Model 5")
colnames(data.nef) = c("Metabolite","Model 1","Model 2","Model 3","Model 4","Model 5")
colnames(data.lci) = c("Metabolite","Model 1","Model 2","Model 3","Model 4","Model 5")
colnames(data.uci) = c("Metabolite","Model 1","Model 2","Model 3","Model 4","Model 5")

my.order = order(data.nes$`Model 1`)
data.nes = data.nes[my.order,]
data.padj = data.padj[my.order,]
data.nef = data.nef[my.order,]
data.lci = data.lci[my.order,]
data.uci = data.uci[my.order,]

library(reshape)
nes = melt(data.nes)
padj = melt(data.padj)
pnef = melt(data.nef)
lci = melt(data.lci)
uci = melt(data.uci)

nes$stars = ""
nes$stars[which(padj$value<0.05)] = "p<0.05"
nes$stars[which(pnef$value<=0.2)] = "NEF<0.2"
nes$stars[which(pnef$value<=0.05)] = "NEF<0.05"

nes$lci = lci$value
nes$uci = uci$value

nes$Metabolite = factor(nes$Metabolite, levels = unique(data.nes$Metabolite))
range(nes$value, na.rm = T)
library(RColorBrewer)
library(svglite)
hm.palette <- rev(brewer.pal(11, "RdBu"))
range(nes$value)

gg <- ggplot(data = nes, aes(x = Metabolite, y = value, ymin=lci, ymax = uci)) + 
  geom_pointrange(aes(color = stars)) +
  coord_flip()+
  geom_hline(yintercept=1, lty=2) +
  ylab("Odds ratio and 95% confidence interval")+
  theme()+ # minimal theme
  facet_wrap(.~variable,ncol = 5)+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, 
                                   size = 9, face = "bold", colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(face = "bold", size = 10, color = "black"),
        axis.text.y = element_text(face = "bold", size = 12, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 10, color = "black"),
        strip.text = element_text(face = "bold", size = 12, color = "black"),
        plot.caption = element_text(hjust = 1,size = 10, face = "bold"))+
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold", 
                                    size = 10), legend.position="bottom")
gg
ggsave(filename=file.path("./results/heatmap_sig_metabolites_main_all_models.png"), 
       device = "png", dpi = 600, height = 4, width = 9) 

ggsave(filename=file.path("./results/heatmap_sig_metabolites_main_all_models.svg"), 
       device = "svg", dpi = 600, height = 4, width = 9) 

print("Figure: ./results/heatmap_sig_metabolites_main_all_models.png")
print("Done")
