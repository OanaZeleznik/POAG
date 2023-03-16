rm(list = setdiff("main.dir", ls()))
print("Running Metabolite Set Enrichment Analysis among all participants ...")

library(Biobase)
library(fgsea)
load("./data/poag_ukbb_es_strata.RData")

fdata = fData(poag.ukb)
mdata = exprs(poag.ukb)
pdata = pData(poag.ukb)

met.classes = vector("list",length = length(unique(fdata$Group)))
names(met.classes) = unique(fdata$Group)

for(i in 1:length(met.classes)){
  met.classes[[i]] = fdata$title[which(fdata$Group == names(met.classes)[i])]
}

runFGSEA = function(pathways = met.classes, stats = ovca.mets, minSize=3, maxSize=500, file.result){
  fgsea.ovca <- fgseaMultilevel(pathways = pathways, stats = stats, minSize=minSize, maxSize=maxSize, eps=0)
  fgsea.ovca = fgsea.ovca[order(fgsea.ovca$padj),]
  to.print = fgsea.ovca[,-8]
  write.table(x = to.print, file = file.result ,col.names = T, row.names = F, sep = "\t", quote = F)
  return(to.print)
  
}

set.seed(1348)
# model 1
load("./results/uclog_all_new.RData")
res.m1 = uclog.reg.res.m2
res.help = merge(x = res.m1, y = fdata)
m1.mets = res.m1$estimate
names(m1.mets) = res.m1$METABOLITE

m1.fgsea = runFGSEA(pathways = met.classes, stats = m1.mets, file.result = "./results/gsea_poag_all.csv")
m1.fgsea$Strata = "All"

m1.fgsea$stars = ""
m1.fgsea$stars[which(m1.fgsea$padj<0.2)] = "*"
m1.fgsea$stars[which(m1.fgsea$padj<0.05)] = "**"
m1.fgsea$stars[which(m1.fgsea$padj<0.001)] = "***"

m1.fgsea = m1.fgsea[order(m1.fgsea$NES),]
m1.fgsea$pathway = factor(m1.fgsea$pathway, levels = m1.fgsea$pathway)

write.table(x = m1.fgsea, file = "results/poag_all_msea.csv", quote = F, sep = ",",row.names = T, col.names = T)
print("File with all results: ./results/poag_all_msea.csv")

library(ggplot2)
library(RColorBrewer)
hm.palette <- rev(brewer.pal(11, "RdBu"))

gg <- ggplot(data = m1.fgsea, aes(x = Strata, y = pathway, fill = NES)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_gradientn(colours = hm.palette, limits = c(-3.5,3.5))+
  geom_text(aes(x = Strata, y = pathway, label=stars), color="black", size=5, 
            nudge_x = 0, nudge_y = -0.3, fontface = "bold") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, 
                                   size = 9, face = "bold", colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "bold", size = 12, color = "black"),
        legend.title = element_text(face = "bold", size = 11, color = "black"),
        legend.text = element_text(face = "bold", size = 10, color = "black"),
        plot.caption = element_text(hjust = 1,size = 10, face = "bold"))+
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold", 
                                    size = 10), legend.position="bottom")+
  guides(fill=guide_colorbar(title = "MESA Enrichment Score",title.position = "top",nrow = 1, label.position = "bottom"))+
  labs(caption = "* FDR<0.2; ** FDR<0.05; *** FDR<0.001")

ggsave(filename=file.path("./results/heatmap_poag_MSEA.png"), device = "png", dpi = 600, 
       width = 8, height = 6, plot = gg)  
print(paste("Figure: ", file.path("./results/heatmap_poag_MSEA.png"), sep = ""))
print("Done")
