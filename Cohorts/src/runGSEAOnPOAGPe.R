rm(list = setdiff("main.dir", ls()))

library(Biobase)
library(fgsea)
library(ggplot2)

print("Running MSEA on results from analysis restricted to those with POAG with paracentral visual field loss ...")

load("./results/main_analysis_all_models.RData")
res.m1 = res.all

# MSEA data and function
load("./data/mets_eset_with_covars_example.RData")
hmdb.groups.raw = fData(mets.glauc)
hmdb.groups.raw$class = hmdb.groups.raw$class_broad

res.m1.annot.new = merge(x = res.m1, 
                         y = hmdb.groups.raw[,c(2,7)], 
                         by.x = "METABOLITE", by.y = "METABOLITE", all.x = T)

length(which(is.na(res.m1.annot.new$class)))

# class annotation is not complete:
res.m1.annot.new$class[grep(pattern = "^CAR", x = res.m1.annot.new$METABOLITE)] = "Carnitines"
length(which(is.na(res.m1.annot.new$class)))

met.classes = vector("list",length = length(unique(res.m1.annot.new$class)))
names(met.classes) = unique(res.m1.annot.new$class)
met.classes[["NA"]] = NULL

for(i in 1:length(met.classes)){
  met.classes[[i]] = res.m1.annot.new$METABOLITE[which(res.m1.annot.new$class == names(met.classes)[i])]
}

runFGSEA = function(pathways = met.classes, stats = ovca.mets, minSize=3, maxSize=500, file.result){

  fgsea.ovca <- fgseaMultilevel(pathways = pathways, stats = stats, minSize=minSize, maxSize=maxSize, eps=0)
  fgsea.ovca = fgsea.ovca[order(fgsea.ovca$padj),]
  to.print = fgsea.ovca[,-8]
  write.table(x = to.print, file = file.result ,col.names = T, row.names = F, sep = "\t", quote = F)
  return(to.print)
  
}

set.seed(1429)

# model 1
load("./results/uclog_glauc_pe_m1_table_raw.RData")
res.m1 = result.data
m1.mets = res.m1$estimate
names(m1.mets) = res.m1$METABOLITE
m1.mets = m1.mets[which(!is.na(m1.mets))]
m1.fgsea = runFGSEA(pathways = met.classes, stats = m1.mets, file.result = "./results/gsea_poag_pe_m1.csv")

# model 2
load("./results/uclog_glau_pe_m2_table_raw.RData")
res.m2 = result.data 
m2.mets = res.m2$estimate
names(m2.mets) = res.m2$METABOLITE
m2.mets = m2.mets[which(!is.na(m2.mets))]
m2.fgsea = runFGSEA(pathways = met.classes, stats = m2.mets, file.result = "./results/gsea_poag_pe_m2.csv")

# model 3
load("./results/uclog_glau_pe_m3_table_raw.RData")
res.m3 = result.data
m3.mets = res.m3$estimate
names(m3.mets) = res.m3$METABOLITE
m3.mets = m3.mets[which(!is.na(m3.mets))]
m3.fgsea = runFGSEA(pathways = met.classes, stats = m3.mets, file.result = "./results/gsea_poag_pe_m3.csv")

# model 4
load("./results/uclog_glau_pe_m4_table_raw.RData")
res.m4 = result.data
m4.mets = res.m4$estimate
names(m4.mets) = res.m4$METABOLITE
m4.mets = m4.mets[which(!is.na(m4.mets))]
m4.fgsea = runFGSEA(pathways = met.classes, stats = m4.mets, file.result = "./results/gsea_poag_pe_m4.csv")

# model 5
load("./results/uclog_glau_pe_m5_table_raw.RData")
res.m5 = result.data
m5.mets = res.m5$estimate
names(m5.mets) = res.m5$METABOLITE
m5.mets = m5.mets[which(!is.na(m5.mets))]
m5.fgsea = runFGSEA(pathways = met.classes, stats = m5.mets, file.result = "./results/gsea_poag_pe_m5.csv")

fgsea.res.h1 = merge(x = m1.fgsea[,c(1,2,3,6)], by.x = "pathway", y = m2.fgsea[,c(1,2,3,6)], by.y = "pathway", suffixes = c(".m1",""))
fgsea.res.h2 = merge(x = fgsea.res.h1, by.x = "pathway", y = m3.fgsea[,c(1,2,3,6)], by.y = "pathway", suffixes = c(".m2",""))
fgsea.res.h3 = merge(x = fgsea.res.h2, by.x = "pathway", y = m4.fgsea[,c(1,2,3,6)], by.y = "pathway", suffixes = c(".m3",""))
fgsea.res = merge(x = fgsea.res.h3, by.x = "pathway", y = m5.fgsea[,c(1,2,3,6)], by.y = "pathway", suffixes = c(".m4",".m5"))

write.table(fgsea.res, file = "./results/fgsea_main_poag_pe_models.csv", row.names = F, col.names = T, quote = F, sep = "\t")

# plot
data.padj = fgsea.res[,c(1,3,6,9,12,15)]
data.nes = fgsea.res[,c(1,4,7,10,13,16)]
colnames(data.nes) = c("Pathway","Model 1","Model 2","Model 3","Model 4","Model 5")
colnames(data.padj) = c("Pathway","Model 1","Model 2","Model 3","Model 4","Model 5")
load(file = "./results/metabolite_classes_plot_order.RData")
my.order = order(data.nes$`Model 1`)
data.nes = data.nes[my.order,]
data.nes$Pathway[which(data.nes$Pathway %in% my.order.sig.m4)] = paste("+",data.nes$Pathway[which(data.nes$Pathway %in% my.order.sig.m4)])
data.padj = data.padj[my.order,]

library(reshape)
nes = melt(data.nes)
padj = melt(data.padj)

nes$stars = ""
nes$stars[which(padj$value<=0.2)] = "*"
nes$stars[which(padj$value<=0.05)] = "**"
nes$stars[which(padj$value<=0.001)] = "***"

nes$Pathway = factor(nes$Pathway, levels = data.nes$Pathway)
range(nes$value)
library(RColorBrewer)
library(svglite)
hm.palette <- rev(brewer.pal(11, "RdBu"))
range(nes$value)

gg <- ggplot(data = nes, aes(x = variable, y = Pathway, fill = value)) + #data and lables
  geom_tile(color = "white", size = 0.1) + # what to plot
  scale_fill_gradientn(colours = hm.palette, limits = c(-3,3))+
  geom_text(aes(x = variable, y = Pathway, label=stars), color="black", size=5, 
            nudge_x = 0, nudge_y = -0.3, fontface = "bold") +
  theme_minimal()+ # minimal theme
  ggtitle("B.     POAG with peripheral visual field loss")+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, 
                                   size = 12, face = "bold", colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", colour = "black"),
        title = element_text(face = "bold", colour = "black"),
        axis.text.y = element_text(face = "bold", size = 12, color = "black"),
        legend.title = element_text(face = "bold", size = 11, color = "black"),
        legend.text = element_text(face = "bold", size = 10, color = "black"),
        plot.caption = element_text(hjust = 1,size = 10, face = "bold"))+
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold", 
                                    size = 10), legend.position="bottom")+
  guides(fill=guide_colorbar(title = "MSEA Enrichment Score",title.position = "top",nrow = 1, label.position = "bottom"))+
  labs(caption = "* FDR<0.2; ** FDR<0.05; *** FDR<0.001\n + Metabolite classes identified in the main analysis")

ggsave(filename=file.path("./results/heatmap_GSEA_broad_main_poag_pe_models.png"), 
       device = "png", dpi = 600, width = 10, height = 6, plot = gg)  

ggsave(filename=file.path("./results/heatmap_GSEA_broad_main_poag_pe_models.svg"), 
       device = "svg", dpi = 600, width = 10, height = 6, plot = gg)  

print("Figure: ./results/heatmap_GSEA_broad_main_poag_pe_models.png")
print("Done")
