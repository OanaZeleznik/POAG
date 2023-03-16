rm(list = setdiff("main.dir", ls()))

library(Biobase)
library(ggplot2)

print("Running analysis restricted to those with POAG with paracentral visual field loss - this will take a few minutes ...")

load("data/mets_eset_with_covars_by_subtype_and_residuals.RData")
source("./src/myRegressions.R")

fdata = fData(mets.glaucpa.m1)
pdata = pData(mets.glaucpa.m1)

# categorical variables are adjusted using residuals

# m1
uclog.reg.res.m1 = runRegression(regression.type = "unconditional.logistic.regression",
                             my.eset = mets.glaucpa.m1,
                             outcome = "glaupa", 
                             covariates = c("ageyr"),
                             strata = NULL, allow.na.values.in.exposure = T,
                             file.with.result.prefix = file.path(main.dir, 
                                                                 "/results/uclog_glauc_pa_m1"),
                             store.complete.result = T)

or.table.m1 = calculateORCI(result.data = uclog.reg.res.m1, or.scale = 1, file.with.result.prefix =
                                file.path(main.dir, "/results/uclog_glauc_pa_m1"),
                              estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = T, estimate.CI = F)

uclog.reg.res.m1 = uclog.reg.res.m1[rownames(or.table.m1),]
all.equal(uclog.reg.res.m1$METABOLITE, or.table.m1$METABOLITE)
or.table.m1$p.value = uclog.reg.res.m1$p.value

res.m1 = merge(x = or.table.m1, by.x = "METABOLITE", y = fdata[,c("METABOLITE","missing","meanCV")], by.y = "METABOLITE") 
res.m1 = res.m1[order(res.m1$p.value),]
res.m1$FDR = p.adjust(p = res.m1$p.value,method = "fdr")
res.m1$NEF = ifelse(res.m1$p.value*132>1,1,res.m1$p.value*132) # NEF = 132 see file "./runNEFF.R"

res.m1[which(res.m1$NEF<0.3),]

# model 2: factors that iunfluence metabolite levels
uclog.reg.res.m2 = runRegression(regression.type = "unconditional.logistic.regression",
                             my.eset = mets.glaucpa.m2,
                             outcome = "glaupa", 
                             covariates = c("ageyr",
                                            "ibmils", "acts"), 
                             strata = NULL, allow.na.values.in.exposure = T,
                             file.with.result.prefix = file.path(main.dir, 
                                                                 "/results/uclog_glau_pa_m2"),
                             store.complete.result = T)

or.table.m2 = calculateORCI(result.data = uclog.reg.res.m2, or.scale = 1, file.with.result.prefix =
                           file.path(main.dir, "/results/uclog_glau_pa_m2"),
                           estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = T, estimate.CI = F)

uclog.reg.res.m2 = uclog.reg.res.m2[rownames(or.table.m2),]
all.equal(uclog.reg.res.m2$METABOLITE, or.table.m2$METABOLITE)
or.table.m2$p.value = uclog.reg.res.m2$p.value

res.m2 = merge(x = or.table.m2, by.x = "METABOLITE", y = fdata[,c("METABOLITE","missing","meanCV")], by.y = "METABOLITE") 
res.m2 = res.m2[order(res.m2$p.value),]
res.m2$FDR = p.adjust(p = res.m2$p.value,method = "fdr")
res.m2$NEF = ifelse(res.m2$p.value*132>1,1,res.m2$p.value*132) # NEF = 132 see file "./runNEFF.R"

res.m2[which(res.m2$NEF<0.2),]

# model 3: set 1 of established risk factors for glaucoma
uclog.reg.res.m3 = runRegression(regression.type = "unconditional.logistic.regression",
                                my.eset = mets.glaucpa.m3,
                                outcome = "glaupa", 
                                covariates =c("ageyr",
                                              "ibmils", "acts",
                                              "nses","agemeno"),
                                strata = NULL, allow.na.values.in.exposure = T,
                                file.with.result.prefix = file.path(main.dir, 
                                                                    "/results/uclog_glau_pa_m3"),
                                store.complete.result = T)

or.table.m3 = calculateORCI(result.data = uclog.reg.res.m3, or.scale = 1, file.with.result.prefix =
                              file.path(main.dir, "/results/uclog_glau_pa_m3"),
                            estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = , estimate.CI = T)

uclog.reg.res.m3 = uclog.reg.res.m3[rownames(or.table.m3),]
all.equal(uclog.reg.res.m3$METABOLITE, or.table.m3$METABOLITE)
or.table.m3$p.value = uclog.reg.res.m3$p.value

res.m3 = merge(x = or.table.m3, by.x = "METABOLITE", y = fdata[,c("METABOLITE","missing","meanCV")], by.y = "METABOLITE") 
res.m3 = res.m3[order(res.m3$p.value),]
res.m3$FDR = p.adjust(p = res.m3$p.value,method = "fdr")
res.m3$NEF = ifelse(res.m3$p.value*132>1,1,res.m3$p.value*132) # NEF = 132 see file "./runNEFF.R"

res.m3[which(res.m3$NEF<0.3),]

# model 4: set 2 of established risk factors for exf glaucoma
uclog.reg.res.m4 = runRegression(regression.type = "unconditional.logistic.regression",
                                my.eset = mets.glaucpa.m4,
                                outcome = "glaupa", 
                                covariates =c("ageyr",
                                              "ibmils", "acts",
                                              "nses","agemeno",
                                              "nitras", "alcos", "caffs", "aheis" ,"calors"),
                                strata = NULL, allow.na.values.in.exposure = T,
                                file.with.result.prefix = file.path(main.dir, 
                                                                    "/results/uclog_glau_pa_m4"),
                                store.complete.result = T)

or.table.m4 = calculateORCI(result.data = uclog.reg.res.m4, or.scale = 1, file.with.result.prefix =
                              file.path(main.dir, "/results/uclog_glau_pa_m4"),
                            estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = T, estimate.CI = T)

uclog.reg.res.m4 = uclog.reg.res.m4[rownames(or.table.m4),]
all.equal(uclog.reg.res.m4$METABOLITE, or.table.m4$METABOLITE)
or.table.m4$p.value = uclog.reg.res.m4$p.value

res.m4 = merge(x = or.table.m4, by.x = "METABOLITE", y = fdata[,c("METABOLITE","missing","meanCV")], by.y = "METABOLITE") 
res.m4 = res.m4[order(res.m4$p.value),]
res.m4$FDR = p.adjust(p = res.m4$p.value,method = "fdr")
res.m4$NEF = ifelse(res.m4$p.value*132>1,1,res.m4$p.value*132) # NEF = 132 see file "./runNEFF.R"

uclog.reg.res.m5 = runRegression(regression.type = "unconditional.logistic.regression",
                                my.eset = mets.glaucpa.m5,
                                outcome = "glaupa", 
                                covariates = c("ageyr",
                                               "ibmils", "acts",
                                               "nses"  ,"agemeno",
                                               "nitras", "alcos", "caffs", "aheis" ,"calors"),                                               
                                strata = NULL, allow.na.values.in.exposure = T,
                                file.with.result.prefix = file.path(main.dir, 
                                                                    "/results/uclog_glau_pa_m5"),
                                store.complete.result = T)

or.table.m5= calculateORCI(result.data = uclog.reg.res.m5, or.scale = 1, file.with.result.prefix =
                              file.path(main.dir, "/results/uclog_glau_pa_m5"),
                            estimate = "estimate", l.ci = "l.ci", u.ci = "u.ci", exp.ci = T, estimate.CI = T)

uclog.reg.res.m5 = uclog.reg.res.m5[rownames(or.table.m5),]
all.equal(uclog.reg.res.m5$METABOLITE, or.table.m5$METABOLITE)
or.table.m5$p.value = uclog.reg.res.m5$p.value

res.m5 = merge(x = or.table.m5, by.x = "METABOLITE", y = fdata[,c("METABOLITE","missing","meanCV")], by.y = "METABOLITE") 
res.m5 = res.m5[order(res.m5$p.value),]
res.m5$FDR = p.adjust(p = res.m5$p.value,method = "fdr")
res.m5$NEF = ifelse(res.m5$p.value*132>1,1,res.m5$p.value*132) # NEF = 132 see file "./runNEFF.R"

save.image(file = "./results/image_runMainAnalysisOnPOAG-Pa.RData")
#load("./results/image_runMainAnalysisOnPOAG-Pa.RData")

res.all.h1 = merge(x = res.m1, by.x = "METABOLITE", y = res.m2[,c(-2,-6,-7)], by.y = "METABOLITE", suffixes = c(".m1",""))
res.all.h2 = merge(x = res.all.h1, by.x = "METABOLITE", y = res.m3[,c(-2,-6,-7)], by.y = "METABOLITE", suffixes = c(".m2",""))
res.all.h3 = merge(x = res.all.h2, by.x = "METABOLITE", y = res.m4[,c(-2,-6,-7)], by.y = "METABOLITE", suffixes = c(".m3",""))
res.all.h4 = merge(x = res.all.h3, by.x = "METABOLITE", y = res.m5[,c(-2,-6,-7)], by.y = "METABOLITE", suffixes = c(".m4",".m5"))

res.all = res.all.h4

res.all[which(res.all$NEF.m5<0.2),]
write.table(res.all, file = "./results/main_analysis_all_models_poag_pa.csv", sep = "\t", row.names = F, col.names = T, quote = F)

# plot
res.all.sig = res.all[which(res.all$p.value.m1<0.05 | res.all$p.value.m2<0.05 | res.all$p.value.m3<0.05 | res.all$p.value.m4<0.05 | res.all$p.value.m5<0.05),]
#res.all.sig = res.all[which(res.all$FDR.m1<0.2 | res.all$FDR.m2<0.2 | res.all$FDR.m3<0.2 | res.all$FDR.m4<0.2 | res.all$FDR.m5<0.2),]

# nominal p-value
data.padj = res.all.sig[,c("METABOLITE",colnames(res.all.sig)[grep(pattern = "p.value",x = colnames(res.all.sig),fixed = T)])]
# FDR
data.fdr = res.all.sig[,c("METABOLITE",colnames(res.all.sig)[grep(pattern = "FDR",x = colnames(res.all.sig),fixed = T)])]
# NEF
data.nef= res.all.sig[,c("METABOLITE",colnames(res.all.sig)[grep(pattern = "NEF",x = colnames(res.all.sig),fixed = T)])]

data.nes =  res.all.sig[,c("METABOLITE",colnames(res.all.sig)[grep(pattern = "OR.CI",x = colnames(res.all.sig),fixed = T)])]
data.nes[,-1] = apply(data.nes[,-1],MARGIN = c(1,2), function(x) {return(as.numeric(strsplit(x = as.character(x), split = " ",fixed = T)[[1]][1]))})

colnames(data.nes) = c("Metabolite","Model 1","Model 2","Model 3","Model 4","Model 5")
colnames(data.padj) = c("Metabolite","Model 1","Model 2","Model 3","Model 4","Model 5")
colnames(data.fdr) = c("Metabolite","Model 1","Model 2","Model 3","Model 4","Model 5")
colnames(data.nef) = c("Metabolite","Model 1","Model 2","Model 3","Model 4","Model 5")

my.order = order(data.nes$`Model 1`)
data.nes = data.nes[my.order,]
data.padj = data.padj[my.order,]
data.fdr = data.fdr[my.order,]
data.nef = data.nef[my.order,]

library(reshape)
nes = melt(data.nes)
padj = melt(data.padj)
pfdr = melt(data.fdr)
pnef = melt(data.nef)

nes$stars = ""
nes$stars[which(padj$value<0.05)] = "."
nes$stars[which(pfdr$value<=0.2)] = "*"
nes$stars[which(pnef$value<=0.2)] = "+"

nes$Metabolite = factor(nes$Metabolite, levels = unique(data.nes$Metabolite))
range(nes$value)
library(RColorBrewer)
library(ggplot2)
hm.palette <- rev(brewer.pal(11, "RdBu"))
range(nes$value)
ggplot(data = nes, aes(x = variable, y = Metabolite, fill = value)) + 
  geom_tile(color = "white", size = 0.1) +
  scale_fill_gradientn(colours = hm.palette, limits = c(0,2.1))+
  geom_text(aes(x = variable, y = Metabolite, label=stars), color="black", size=4, 
            nudge_x = 0, nudge_y = 0, fontface = "bold") +
  theme_minimal()+ # minimal theme
  ggtitle("POAG with paracentral visual field loss")+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, 
                                   size = 9, face = "bold", colour = "black"),
        axis.title.x = element_text(face = "bold", colour = "white"),
        axis.title.y = element_text(face = "bold", colour = "black"),
        axis.text.y = element_text(face = "bold", size = 9, color = "black"),
        legend.title = element_text(face = "bold", size = 11, color = "black"),
        legend.text = element_text(face = "bold", size = 10, color = "black"),
        plot.caption = element_text(hjust = 1,size = 10, face = "bold"))+
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold", 
                                    size = 10), legend.position="bottom")+
  guides(fill=guide_colorbar(title = "Odds Ratio",title.position = "top",nrow = 1, label.position = "bottom"))+
  labs(caption = ". p<0.05; * FDR<0.2; + NEF<0.2 ")

ggsave(filename=file.path("./results/heatmap_sig_metabolites_main_all_models_poag_pa.png"), device = "png", dpi = 600, height = 15, width = 10)  
print("Done")
