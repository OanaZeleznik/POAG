################################################################################
# this file only runs if you have access to the real data
################################################################################
rm(list = setdiff("main.dir", ls()))

library(Biobase)
library(ggplot2)
library(RColorBrewer)

print("Creating bubble plot ...")

createBubblePlot = function(my.data, title.plot){
  
  colnames(my.data) = c("metabolite","OR.CI","pval")
  # subset to metabolites with info on # of Carbon atoms and # double bonds
  
  my.data$class = sapply(X = my.data$metabolite, FUN = function(x) {(strsplit(x = x,split = "(", fixed = T)[[1]][1])})
  my.data$help = sapply(X = my.data$metabolite, FUN = function(x) {(strsplit(x = x,split = "(", fixed = T)[[1]][2])})
  
  my.data = my.data[complete.cases(my.data),]
  help.classes = table(my.data$class)
  
  classes.to.keep = names(help.classes)[which(help.classes>3)]
  my.data = my.data[which(my.data$class %in% classes.to.keep),] 
  
  my.data$nC = as.numeric(sapply(X = my.data$help, FUN = function(x) {(strsplit(x = x,split = ":", fixed = T)[[1]][1])}))
  my.data$nDB.help = sapply(X = my.data$help, FUN = function(x) {(strsplit(x = x,split = ":", fixed = T)[[1]][2])})
  my.data$nDB = as.numeric(sapply(X = my.data$nDB.help, FUN = function(x) {gsub(pattern = ")",replacement = "",x = x,fixed = T)}))
  my.data$OR = as.numeric(sapply(X = my.data$OR.CI, FUN = function(x) {(strsplit(x = x,split = "(", fixed = T)[[1]][1])}))
  my.data$beta = log(my.data$OR, base = exp(1))
  
  
  to.plot = my.data[c("metabolite","class","nC","nDB","beta","pval","OR")]
  to.plot$log10p = -log10(to.plot$pval)
  
  table(to.plot$class)
  
  curr.data = to.plot
  hm.palette <- rev(brewer.pal(11, "RdBu"))
  curr.data = curr.data[which(curr.data$class == "TG"),]
  range(curr.data$OR)
  
  my.plots <- ggplot(curr.data, aes(x = nC, y = nDB, size = log10p, fill = OR))+
    geom_point(shape = 21)+
    facet_wrap(~class, nrow = 5, ncol = 3,scales = "free")+
    scale_fill_gradientn(colours = hm.palette, limits = c(0.6,1.4), breaks = c(0.7,1,1.3))+
    #scale_fill_gradient2(midpoint = 1, low = "#053061",high = "#67001f",mid = "#f7f7f7")+#, limits = c(0,5), breaks = c(0,1,2,3,4,5))+
    scale_size_area(max_size = 5, limits = c(0,3), breaks = 0:3)+  
    ggtitle(title.plot)+
    labs(x = "Number of Carbon atoms", y = "Number of double bonds")+
    theme(axis.text.x = element_text(angle = 0, vjust = 0, size = 10, hjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold", size = 10),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(face = "bold", size = 10),
          legend.key.size = unit(0.75,"cm"),
          title = element_text(face = "bold", size = 12),
          strip.text = element_text(face = "bold", size = 12, hjust = 0.05), legend.position = "bottom", legend.direction = "horizontal")+
    guides(fill = guide_colorbar(order = 1, title = "OR", title.vjust = +0.5), size = guide_legend(order = 2,title = "-log10(p)"))
   
  my.file = paste(main.dir,"results/",gsub(pattern = " ",replacement = "_",x = title.plot, fixed = T),".png", sep = "")
  ggsave(filename = my.file, plot = my.plots, device = "png",width = 7, height = 5)
  
}

load("./results/clog_all_m5_table_raw.RData")
res.m5 = result.data
m5.mets = res.m5$estimate
res.m5$OR = as.character(exp(res.m5$estimate))
names(m5.mets) = res.m5$METABOLITE

poag.all = res.m5[,c(1,12,8)]

createBubblePlot(my.data = poag.all,
                 title.plot = "POAG")

print("Figure: ./results/POAG.png")
print("Done")

