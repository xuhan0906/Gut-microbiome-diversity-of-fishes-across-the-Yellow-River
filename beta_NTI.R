library(NST)
library(picante)
library(ape)
 
#
otu <- read.csv('community.csv', row.names = 1)
otu <- as.data.frame(t(otu))
 
group <- read.csv('treatment.csv', row.names = 1)
 
#
tree <- read.tree('tree.nwk')
 
#
tree <- prune.sample(otu, tree)
 

set.seed(123)
pnst <- pNST(comm = otu, tree = tree, group = group, phylo.shuffle = TRUE, taxo.null.model = NULL, 
    pd.wd = tempdir(), abundance.weighted = TRUE, rand = 1000, nworker = 16, SES = TRUE, RC = FALSE)
 
#
#names(pnst)
betaMNTD <- pnst$index.pair
head(betaMNTD)
 
#
write.csv(betaMNTD, 'betaMNTD.csv', quote = FALSE, row.names= FALSE)
