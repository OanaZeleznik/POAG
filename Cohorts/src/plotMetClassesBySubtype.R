rm(list = ls())
main.dir = "my_directory"
main.dir = "/udd/nhotz/projects/POAG_NatComm/Cohorts/"
setwd(main.dir)

library(Biobase)
library(chanmetab)
library(ggplot2)


all.fgsea = read.csv(file = "./results/fgsea_main_all_models.csv", header = T, sep = "\t")
m5.all = all.fgsea[,c(1,14:16)]
m5.all$model = "Model 5"
m5.all$strata = "POAG"
colnames(m5.all) = gsub(pattern = ".m5",replacement = "",x = colnames(m5.all))

pa.fgsea = read.csv(file = "./results/fgsea_main_poag_pa_models.csv", header = T, sep = "\t")
m5.pa = pa.fgsea[,c(1,14:16)]
m5.pa$model = "Model 5"
m5.pa$strata = "POAG with \n paracentral VF loss"
colnames(m5.pa) = gsub(pattern = ".m5",replacement = "",x = colnames(m5.pa))


pe.fgsea = read.csv(file = "./results/fgsea_main_poag_pe_models.csv", header = T, sep = "\t")
m5.pe = pe.fgsea[,c(1,14:16)]
m5.pe$model = "Model 5"
m5.pe$strata = "POAG with \n peripheral VF loss"
colnames(m5.pe) = gsub(pattern = ".m5",replacement = "",x = colnames(m5.pe))

mets.file1 = merge(x = m5.all, y = m5.pa, by = "pathway", suffixes = c(".POAG",""))
mets.file2 = merge(x = mets.file1, y = m5.pe, by = "pathway", suffixes = c(".paracentral.VF.los",".peripheral.VF.loss"))


write.table(mets.file2, file = "./results/main_analysis_m5_byVFlos_fgsea.csv", 
            sep = "\t", row.names = F, col.names = T, quote = F)

to.plot = rbind(m5.all, m5.pa, m5.pe)

help = m5.pa
help = help[order(help$NES),]

to.plot$pathway = factor(to.plot$pathway, levels = help$pathway)
to.plot$strata = factor(to.plot$strata, levels = c("POAG with \n paracentral VF loss", "POAG with \n peripheral VF loss", "POAG"))

to.plot$stars = ""
to.plot$stars[which(to.plot$pval<0.05)] = "*"
to.plot$stars[which(to.plot$padj<0.2)] = "**"
to.plot$stars[which(to.plot$padj<0.05)] = "***"

range(to.plot$NES, na.rm = T)
library(RColorBrewer)
hm.palette <- rev(brewer.pal(11, "RdBu"))

gg <- ggplot(data = to.plot, aes(x = strata, y = pathway, fill = NES)) + #data and lables
  geom_tile(color = "white", size = 0.1) + # what to plot
  scale_fill_gradientn(colours = hm.palette, limits = c(-3,3))+#, values = rescale(c(0,1,3.1)))+
  geom_text(aes(x = strata, y = pathway, label=stars), color="black", size=5, 
            nudge_x = 0, nudge_y = -0.3, fontface = "bold") +
  theme_minimal()+ # minimal theme
  #facet_wrap(model~.)+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, 
                                   size = 9, face = "bold", colour = "black"),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "bold", size = 12, color = "black"), 
        strip.text = element_text(face = "bold", size = 12, color = "black"),
        legend.title = element_text(face = "bold", size = 11, color = "black"),
        legend.text = element_text(face = "bold", size = 11, color = "black"),
        plot.caption = element_text(hjust = 1,size = 11, face = "bold"))+
  theme(strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold", 
                                    size = 12), legend.position="bottom")+
  guides(fill=guide_colorbar(title = "MSEA Enrichment Score",title.position = "top",nrow = 1, label.position = "bottom"))+
  labs(caption = "* p<0.05; ** FDR<0.2; *** FDR<0.05\n VF: visual field")

ggsave(filename=file.path("./results/heatmap_met_classes_by_VFlos.png"), 
       device = "png", dpi = 600, width = 9, height = 7, plot = gg)  

