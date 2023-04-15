library(synapser)
library(tidyverse)

synLogin()

# enumerate results directory on synapse
parent_id <- 'syn51117833'

## enumerate results directory on synapse
  # syn.dir <- synapser::synGetChildren(parent_id)$asList() %>% 
  #   tibble(f = .) %>% unnest_wider(f) 
  
  # kda.dir <- map_dfr(
  #   syn.dir %>% filter(grepl('Folder', type)) %>% pull(id),
  #   ~ {
  #     id = synGetChildren(.x)$asList() %>% tibble(f = .) %>% unnest_wider(f) %>%
  #       filter(grepl('Folder', type)) %>% pull(id)
  #     synGetChildren(id)$asList() %>% tibble(f=.) %>% unnest_wider(f)
  #     }
  # )

# local dirs
d <- list.dirs(paste0(here::here(),'/results'), recursive = F) %>% 
  str_subset('input_gene_lists|wKDA', negate=T)

for( j in 1:length(d) ){
  
  # top-level directory
  foo <- synStore( Folder(d[j] %>% str_remove_all('^.*results/'), parent = parent_id) )
  
  # list graph files
  f <- list.files(d[j]) %>% str_subset('kda',negate=T)

  # store graph files  
  for(i in 1:length(f)){
    foo2 <- synStore( File(
      paste0(d[j],'/',f[i]),
      parent=foo$properties$id
    ))
  }
  
  # kda sub-directory
  foo3 <- synStore( Folder('kda', parent = foo$properties$id) )

  # list kda files 
  f <- list.files( paste0(d[j],'/','kda')) 
  
  # store kda files
  for(i in 1:length(f)){
    foo4 <- synStore( File(
      paste0(d[j],'/kda/',f[i]),
      parent=foo3$properties$id
    ))
  }
  
}