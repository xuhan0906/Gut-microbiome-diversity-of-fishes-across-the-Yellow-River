suppressWarnings(suppressMessages(library(amplicon)))
library(ggplot2)
library(RColorBrewer)
library(randomcoloR)
color<-c("#002e6d",
         "#ffc747",
         "#539cff",
         "#beba27",
         "#5a005f",
         "#01e294",
         "#ec50aa",
         "#01942e",
         "#e43a8a",
         "#8ec844",
         "#ff84d7",
         "#02ceaa",
         "#a30042",
         "#dfd86e",
         "#ffa9f9",
         "#716b00",
         "#ff719b",
         "#a5914a",
         "#64000d",
         "#ffb86e",
         "#883042",
         "#f75a53",
         "#924100",
         "#ff7b84")
#
distance_mat <- read.delim('bray.tsv', row.names = 1, sep = '\t', stringsAsFactors = F, check.names = F)
metadata=read.table("group.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors=F)
#
beta_pcoa(distance_mat,metadata,groupID="species",ellipse=T)+geom_point(alpha = 0.0001)+
  stat_ellipse(level = 0.68,alpha=1,size=0.2)+
  scale_color_manual(values = color)+theme_bw()
ggsave(paste0("p8.PCoA.unifrac.pdf"), p, width=12, height=8)



