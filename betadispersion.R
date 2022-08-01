#######The following codes reference published literature
#####Youngblut, Nicholas D., et al. "Host diet and evolutionary history explain different aspects of gut microbiome diversity among vertebrate clades." Nature communications 10.1 (2019): 1-15.
library(cluster)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(phyloseq)
library(vegan)
library(doParallel)
######Convert four beta diversity distances into RDS files
files = list.files(pattern='*.tsv', full.names=TRUE)
for(f in files){
  data=read.delim(f, sep='\t', row.names=1)
  final_file = gsub('.tsv', '.RDS', f)
  saveRDS(data, final_file)
}
beta_div = list(
  'bray_curtis' = 'bray_curtis.RDS',
  'jaccard' = 'jaccard.RDS',
  'unweighted_unifrac' = 'unweighted_unifrac.RDS',
  'weighted_unifrac' = 'weighted_unifrac.RDS'
)
beta_div_rds = list()
for(n in names(beta_div)){
  x = readRDS(beta_div[[n]])
  beta_div_rds[[n]] = x
}

beta_div_rds_file = file.path('beta_div.RDS')
saveRDS(beta_div_rds, beta_div_rds_file)
######
# load beta-diversity dist. matrices
beta_d = readRDS("beta_div.RDS")
metadata<- read.delim('metadata.tsv', sep = '\t',row.names = 1)

dist_mtx_order = function(d, x){
  # Ordering distance matrixes
  # d = distance matrix (dist class)
  # x = vector to order dist. obj. by
  m = d %>% as.matrix
  d = as.dist(m[x,x])
  return(d)
}

beta_disp = function(d, fac, samples){
  # checking beta dispersion
  d = dist_mtx_order(d, samples)
  ret = vegan::betadisper(d,fac)
  ret = data.frame(group = ret$group,         # note: group is in wrong order; but can just match host_subj_id based on sampleID
                   distances = ret$distances)
  return(ret)
}
format_data = function(x){
  # format output
  x = do.call(rbind, x) 
  x$metric = gsub('\\..+', '', rownames(x))
  x$sample_id = gsub('^.*?\\.', '', rownames(x))
  return(x)
}
.beta_disp_all = function(group, d, metadata, threads=2){
  #doParallel::registerDoParallel(cores=threads)
  ret = plyr::llply(d, beta_disp, fac=metadata[,group], 
                    samples=rownames(metadata), .parallel=TRUE)
  ret = format_data(ret)
  ret$grouping = group %>% as.character
  return(ret)
}
beta_disp_all = function(list_of_dists, metadata, groups, threads=2){
  groups = as.list(groups)
  ret = lapply(groups, .beta_disp_all, d=list_of_dists, metadata=metadata, threads=threads)
  ret = do.call(rbind, ret)
  return(ret)
}
groups = c('host_subject_id', 'scientific_name', 'genus', 'family', 'order', 'class') #host_subject_id¿ÉÒÔÉ¾³ý
beta_disp_ret = beta_disp_all(beta_d, metadata, groups=groups, threads=threads)

# formatting
## renaming grouping
old_names = c('host_subject_id', 'scientific_name','genus', 'family', 'order', 'class')
new_names = c('Individual', 'Species','Genus', 'Family', 'Order', 'Class')
df = data.frame(old_names, new_names)
beta_disp_ret_j = beta_disp_ret %>%
  inner_join(df, c('grouping'='old_names')) %>%
  dplyr::select(-grouping) %>%
  rename('grouping' = new_names)
beta_disp_ret_j$grouping = factor(beta_disp_ret_j$grouping, levels=new_names)
## Distance metrics
old_names = c('jaccard', 'bray_curtis', 'unweighted_unifrac', 'weighted_unifrac')
new_names = c('Jaccard', 'Bray Curtis', 'Unweighted\nUnifrac', 'Weighted\nUnifrac')
df = data.frame(old_names, new_names)
beta_disp_ret_j = beta_disp_ret_j %>%
  inner_join(df, c('metric'='old_names')) %>%
  dplyr::select(-metric) %>%
  rename('metric' = new_names)
## filtering out any with just one in group
beta_disp_ret_j = beta_disp_ret_j %>%
  group_by(group, metric, grouping) %>%
  mutate(N_SAMP = n()) %>%
  ungroup() %>%
  filter(N_SAMP > 1) %>%
  mutate(metric = factor(metric, levels=new_names))
########delete individual
beta_disp_ret_j<-beta_disp_ret_j[-c(1:10),]
########plot
p = ggplot(beta_disp_ret_j, aes(grouping, distances)) +
  geom_boxplot() +
  labs(x='Taxonomic level', y='Beta dispersion') +
  facet_grid(metric ~ .) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )
##########
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
plot_write(p, file='betadisper_diet_all-dist_tax-lev.pdf')


