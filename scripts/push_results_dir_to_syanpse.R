library(synapser)
library(tidyverse)

synLogin()

parent_id <- 'syn51117833'

d <- list.dirs(paste0(here::here(),'/results'), recursive = F) %>% 
  str_subset('input_gene_lists|wKDA', negate=T)

for( j in 1:length(d) ){
  foo <- synStore( Folder(d[j] %>% str_remove_all('^.*results/'), parent = parent_id) )
  f <- list.files(d[j]) %>% str_subset('kda',negate=T)
  for(i in 1:length(f)){
    foo2 <- synStore( File(
      paste0(d[j],'/',f[i]),
      parent=foo$properties$id
    ))
  }
  foo3 <- synStore( Folder('kda', parent = foo$properties$id) )
  f <- list.files( paste0(d[j],'/','kda')) 
  for(i in 1:length(f)){
    foo4 <- synStore( File(
      paste0(d[j],'/kda/',f[i]),
      parent=foo3$properties$id
    ))
  }
}

