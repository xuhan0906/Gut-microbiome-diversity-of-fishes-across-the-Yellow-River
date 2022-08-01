library(vegan)
library(picante)
library(geosphere)
library(tidyverse)
library(ecodist)
####read table
geo <- read.table('site.txt', header= TRUE,row.names = 1)
otu<-read.table("site_otu_new.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors=F)
otu <- data.frame(t(otu))
dist.abund <- vegdist(otu, method = 'bray')
####match sample name
geo<-geo[otu$sample.id,]
#####
d.geo <- distm(geo, fun = distHaversine)       #library(geosphere)
dist.geo <- as.dist(d.geo)
#####
MRM(dsit.geo~dist.abund,nperm = 10000,mrank = T)

aa <- as.vector(dist.abund)

gg <- as.vector(dist.geo)
mat <- data.frame(aa, gg)
mm <- ggplot(mat, aes(y = aa, x = gg/1000)) + 
  geom_point(size = 3, alpha = 0.5, colour="black") +
  geom_smooth(method = "lm", colour = "blue",alpha=0.2) + 
  labs(x = "Physical separation (km)", y = "Bray-Curtis Dissimilarity") + 
  theme( axis.text.x = element_text(colour = "black", size = 12), 
         axis.text.y = element_text(size = 11, colour = "black"), 
         axis.title= element_text(size = 14, colour = "black"),
         panel.border = element_rect(fill = NA, colour = "black"))+theme_bw()
mm
