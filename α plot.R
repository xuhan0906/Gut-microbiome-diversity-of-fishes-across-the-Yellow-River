library(picante)  
alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson Ö¸Êý
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}

#
otu <- read.delim('otu_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- t(otu)
tree <- read.tree('otu_tree.tre')

#not include phylogenetic diversity
alpha_all <- alpha(otu, base = 2)

#rarecurve
rarecurve(otu, step = 2000,label=F,col = '#01e294')
#
write.csv(alpha_all, 'alpha.csv', quote = FALSE)

group=read.delim('group.txt')
alpha_all$sample_ID=row.names(alpha_all)
Data<-merge(x=alpha_all,y=group,by.x = "sample_ID")
##plot
#Data <- read.delim('alpha.txt', row.names = 1, sep = '\t')
library(ggsci)
library(ggpubr)
color=c("#baffa1",
        "#160090",
        "#92d000",
        "#af3aec",
        "#4eff91",
        "#8465ff",
        "#00ce65",
        "#cc0092",
        "#087e00",
        "#7c8bff",
        "#f18100",
        "#001e6c",
        "#ffdf74",
        "#67004b",
        "#cefff6",
        "#b90006",
        "#5ccfff",
        "#6f001f",
        "#009075",
        "#ffb6fe",
        "#4a2b00",
        "#ffdcf4",
        "#004944",
        "#004262")
my_comparisons <- list(c("Fgut", "Hgut"), c("Fgut", "Mgut"), c("Mgut", "Hgut"))
ggboxplot(Data, x="species", y="Shannon", fill = 'species')+
 stat_compare_means(method = "kruskal.test",label.y = 10)+scale_fill_manual(values = color)+
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggboxplot(Data, x="geography", y="Shannon", fill = 'geography')+
  stat_compare_means(method = "kruskal.test",label.y = 10)+scale_fill_manual(values = color)+
  theme(axis.text.x = element_text(angle=45, hjust=1))

