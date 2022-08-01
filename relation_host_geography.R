library(vegan)
library(picante)
library(geosphere)
library(tidyverse)
library(ecodist)
####read table
site <- read.table('site.txt', header= TRUE)
tree<-read.tree('COI.contree')
tree =multi2di(tree)
metadata<-read.delim('sample-metadata.txt')
geo<-inner_join(site,metadata,by=c('sample_ID'='sample_ID')) %>% select(species,longitude,latitude) %>% column_to_rownames(var='species')
####match sample name
geo<-geo[tree$tip.label,]
#####
d.geo <- distm(geo, fun = distHaversine)       #library(geosphere)
dist.geo <- as.dist(d.geo)
####host phylogenetic distance
rescale_dist_mtx = function(m){
  m = m %>% as.matrix
  labs = m %>% colnames
  n_row = m %>% nrow
  n_col = m %>% ncol
  x = m %>% as.vector 
  x = scales::rescale(x) 
  m = matrix(x, nrow=n_row, ncol=n_col)
  colnames(m) = labs
  rownames(m) = labs
  m = m %>% as.dist
  return(m)
}
host_D = tree %>% cophenetic %>% as.dist

MRM(dist.geo~host_D,nperm = 10000,mrank = T)

aa <- as.vector(host_D)

gg <- as.vector(dist.geo)
mat <- data.frame(aa, gg)

 ggplot(mat, aes(y = aa, x = gg/1000)) + 
  geom_point(size = 1, alpha = 0.5, colour="purple") +
  geom_smooth(method = "lm", colour = "blue",alpha=0.2) + theme_bw()


