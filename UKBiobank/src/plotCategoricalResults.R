rm(list = setdiff("main.dir", ls()))
print("Ploting results ...")

library(Biobase)
library(ggplot2)
library("svglite")

processResultsToPlot = function(my.res, my.strata){
    
  my.res$strata = my.strata
  my.res$OR = as.numeric(sapply(X = my.res$OR.CI, FUN = function(x){strsplit(x = x, split = " ", fixed = T)[[1]][1]}))
  my.res$help = sapply(X = my.res$OR.CI, FUN = function(x){strsplit(x = x, split = " ", fixed = T)[[1]][2]})
  my.res$help = sapply(X = my.res$help, FUN = function(x){gsub(pattern = "(", replacement = "", x = x, fixed = T)})
  my.res$help = sapply(X = my.res$help, FUN = function(x){gsub(pattern = ")", replacement = "", x = x, fixed = T)})
  my.res$lCI = as.numeric(sapply(X = my.res$help, FUN = function(x){strsplit(x = x, split = "-", fixed = T)[[1]][1]}))
  my.res$uCI = as.numeric(sapply(X = my.res$help, FUN = function(x){strsplit(x = x, split = "-", fixed = T)[[1]][2]}))
    
  my.help = my.res[,c("METABOLITE","OR","lCI","uCI","p.value","FDR","NEF","strata")]  
  
  return(my.help)
}

load("./results/uclog_all_new.RData")
all.res = processResultsToPlot(my.res = res.m2, my.strata = "All")

to.plot = rbind(all.res)

to.plot$stars = ""
to.plot$stars[which(to.plot$p.value<0.05)] = "p<0.05"
to.plot$stars[which(to.plot$NEF<0.2)] = "NEF<0.2"
to.plot$stars[which(to.plot$NEF<0.05)] = "NEF<0.05"
table(to.plot$stars)

to.plot$stars = factor(to.plot$stars, levels = c("p<0.05", "NEF<0.2", "NEF<0.05"))

help = all.res[order(all.res$OR),]
help.v2 = help[which(help$p.value <0.05 | help$NEF < 0.2 | help$FDR < 0.2),]
mets.to.plot = help.v2$METABOLITE
to.plot = to.plot[which(to.plot$METABOLITE %in% mets.to.plot),]
to.plot$METABOLITE = factor(to.plot$METABOLITE, levels = help.v2$METABOLITE)

to.plot$strata = factor(to.plot$strata, levels = c("All","Men","Women","BMI <=25", "BMI >25","Age <=58","Age >58"))

range(to.plot$OR)
library(RColorBrewer)
hm.palette <- rev(brewer.pal(11, "RdBu"))

gg <- ggplot(data = to.plot, aes(x = METABOLITE, y = OR, ymin=lCI, ymax = uCI)) + #data and lables
  geom_pointrange(aes(color = stars)) + # what to plot
  coord_flip()+
  geom_hline(yintercept=1, lty=2) +
  ylab("Odds ratio and 95% confidence interval")+
  theme()+ # minimal theme
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

ggsave(filename=file.path("./results/heatmap_poag_new.svg"), device = "svg", dpi = 600, 
       width = 8, height = 6, plot = gg)  

ggsave(filename=file.path("./results/heatmap_poag_new.png"), device = "png", dpi = 600,
       width = 8, height = 6, plot = gg)

print("Figure: ./results/heatmap_poag_new.png")

cat("Done")
