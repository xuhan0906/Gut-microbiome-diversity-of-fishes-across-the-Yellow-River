library(tidyverse)
library(ggthemes)

#read table
data1<-read.delim('contribution.txt')

#

data2<-data1 %>% pivot_longer(-name,names_to='type',values_to='values')

#barplot
data2 %>% ggplot(aes(type, values, fill = name)) +geom_col(position = 'fill', width = 0.6,colour='black') +
  scale_fill_manual(values=rev(c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3')))+theme_classic()
 
