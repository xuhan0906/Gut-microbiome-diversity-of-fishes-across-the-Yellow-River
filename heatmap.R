library(pheatmap)
library(RColorBrewer)
ht=read.delim("lipa_sort.txt",row.names = 1)
group=read.delim("group.txt",row.names = 1)
col_color<-c("#6e6ec7",
             "#bfa13d",
             "#a34e97",
             "#6ca24d",
             "#bc4862",
             "#46c19a",
             "#b65c36")
names(col_color)<-unique(group$group)
anColor<-list(group=col_color)
pheatmap(ht,na_col = "white",cluster_rows=F,cluster_cols=F,annotation_col = group,
         annotation_colors = anColor,border=F,color = colorRampPalette(c("#fdbb2d", "#1E9600"))(10))
