#######The following codes reference published literature
#####Youngblut, Nicholas D., et al. "Host diet and evolutionary history explain different aspects of gut microbiome diversity among vertebrate clades." Nature communications 10.1 (2019): 1-15.
#work_dir="/media/project/2.R/LIPA"
#setwd(work_dir)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ape)
library(paco)
library(picante)
library(future)
library(future.batchtools)
library(doParallel)
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
## load files
host_tree = read.tree("COI.contree")
host_tree =multi2di(host_tree)   #convert to rooted tree
micro_tree=read.tree("micro_tree.nwk")
otu=read.delim("species_group.txt",sep = '\t',row.names = 1)
###abundance
otu<-otu/69000
#filter
otu <- otu[which(rowSums(otu) >= 0.005), ]
####
phy.tree = prune.sample(t(otu), micro_tree)
####
##16S dist
micro_D = phy.tree %>% cophenetic %>% rescale_dist_mtx %>% as.matrix
micro_D %>% dim
##host dist
host_D = host_tree %>% cophenetic %>% rescale_dist_mtx %>% as.matrix
host_D %>% dim
###Preparing data
otu = otu %>% t %>% apply(2, function(x) ifelse(x > 0, 1, 0)) %>% as.matrix
D = prepare_paco_data(H=host_D, P=micro_D, HP=otu)
D %>% names
D = add_pcoord(D, correction='cailliez')
D %>% names
#########
###PACO: diffused model
PACo_file = file.path('physeq_IndD_PACo-Con.RDS')
D %<-% { PACo(D, nperm=1000, seed=3874, method='quasiswap', symmetric=TRUE) } %packages% c("paco")
# save results
saveRDS(D, PACo_file)
cat('File written:', PACo_file, '\n')
D = readRDS(PACo_file)
# goodness of fit
D$gof
#######Individual contributions
# for loading results instead of re-running
PACo_links_file = file.path('physeq_IndD_PACo-Con-links.RDS')
# wrapper
paco_links_run = function(D, threads=16){
  doParallel::registerDoParallel(cores=threads)
  paco::paco_links(D, .parallel=TRUE)
}
# run on cluster
D_links %<-% { paco_links_run(D, threads) } %packages% c("paco", "doParallel")
# saving object
saveRDS(D_links, PACo_links_file)
cat('File written:', PACo_links_file, '\n')
# loading object
D_links = readRDS(PACo_links_file)
# residuals
res = residuals_paco(D_links$proc) %>% as.data.frame 
colnames(res) = 'residuals'
res = res %>%
  mutate(comparison = rownames(.)) %>%
  separate(comparison, c('host', 'microbe'), sep='-') 
#res %>% status
res$microbe %>% unique %>% length
###Formatting output
D_links_jk = do.call(rbind, D_links["jackknife"]) %>%
  t %>% as.data.frame %>%
  mutate(comparison = rownames(.)) %>%
  separate(comparison, c('host', 'microbe'), sep='-') %>%
  inner_join(res, c('host'='host', 'microbe'='microbe'))
#Summarizing results
# adding taxonomy
tax = read.delim("taxonomy.txt")
tax$microbe=tax$OTU
D_links_jk = D_links_jk %>%
  inner_join(tax, c('microbe'))
# summarizing by taxonomy
tmp = D_links_jk %>%
  filter(!is.na(residuals),
         Genus != '') %>%
  group_by(Genus) %>%
  mutate(mean_res = mean(residuals, na.rm=TRUE),
         n = n() %>% as.numeric) %>%
  ungroup() %>%
  mutate(Genus = Genus %>% reorder(mean_res),
         Genus = Genus %>% reorder(Phylum %>% as.factor %>% as.numeric))
ggplot(tmp, aes(Genus, n, fill=Phylum)) +
  geom_bar(stat='identity') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)
  )
##Summarizing residuals
host_tax=read.delim("host_class.txt")
D_links_l = D_links_jk %>%
  inner_join(host_tax, c('host'='species'))
tmp = D_links_l %>% 
  group_by(order) %>%
  mutate(median_resid = median(residuals)) %>%
  ungroup() %>%
  mutate(class = class %>% reorder(-median_resid))
p=ggplot(D_links_l, aes(order, residuals, color=order)) +
  geom_boxplot() +
  scale_color_discrete('order') +
  labs(x='Order', y='Residuals') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )
plot_write = function(p, file=NULL, path=NULL, width=NA, height=NA, ...){
  # file path
  if(is.null(path)){
    path = file.path(getwd(), '.figures')
    if(! dir.exists(path)){
      dir.create(path, showWarnings=FALSE)
    }
  }
  # file name
  if(is.null(file)){
    file = paste0(fig_uuid(), '.pdf')
  } 
  file = file.path(path, file)
  
  # width & height
  if(is.na(width)){
    width = options()$repr.plot.width
  }
  if(is.na(height)){
    height = options()$repr.plot.height
  }
  # writting figure
  if(length(class(p)) >= 2 & class(p)[2] == 'ggplot2'){
    ggplot2::ggsave(filename=file, plot=p, width=width, height=height, ...)
  } else {
    pdf(file=file, width=width, height=height)
    plot(p, ...)
    dev.off()
  }
  cat('File written:', file, '\n')
  # plotting
  plot(p)
}
options(repr.plot.width=5, repr.plot.height=2.5)
plot_write(p, file='PAco.pdf')