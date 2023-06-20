
# setup -------------------------------------------------------------------

library(synapser)
library(igraph)
library(tidyverse)

synLogin()

base_net <- igraph::read_graph(synGet('syn51080932')$path, format = 'graphml') 

nw.edges <- tibble( ea = edge.attributes(base_net) ) %>% 
  t() %>% as_tibble(rownames = NA, .name_repair = 'unique') %>% unnest(everything()) %>% 
  rename_with(., ~names(edge.attributes(base_net)), everything())

# ROSMAP ------------------------------------------------------------------

st = Sys.time()

# download and pivot pcor table
synLogin()
rm.pcor <- read_csv(synGet('syn51061790')$path) %>% 
  rename(gene1 = X1) %>% 
  pivot_longer(cols = -gene1, names_to = 'gene2', values_to = 'pcor') %>% 
  filter( gene1 != gene2 ) 

cat('download: \n')
Sys.time() - st

# get gene symbols corresponding to Ensembl IDs
genes <- union(rm.pcor$gene1, rm.pcor$gene2)
sym <- gprofiler2::gconvert(
  genes, 
  organism = 'hsapiens',
  target = 'HGNC'
) %>% 
  select( gene1 = input, gene2 = input, symbol1 = target, symbol2 = target)

cat('id conversion: \n')
Sys.time() - st

# add symbols to table
# remove genes without symbols (therefore not in NW)
# & put genes in consistent order to reduce dataset size
rm.pcor <- 
  left_join( rm.pcor, sym %>% select(ends_with('1')), by = 'gene1') %>% 
  left_join( .,       sym %>% select(ends_with('2')), by = 'gene2') %>% 
  filter(!is.na(symbol1), !is.na(symbol2)) 

rm.pcor$edge = NA_character_
rm.pcor$in_PC_KB = NA_character_

for(i in 1:nrow(rm.pcor)){
  
  # sort the gene pair
  if( rm.pcor$symbol1[i] > rm.pcor$symbol2[i] ){
    s1 = rm.pcor$symbol2[i]
    s2 = rm.pcor$symbol1[i]
    rm.pcor$symbol1[i] = s1
    rm.pcor$symbol2[i] = s2
  }
  
  # define the edge
  rm.pcor$edge[i] = paste0(rm.pcor$symbol1[i],':',rm.pcor$symbol2[i])
  
  # annotate if the edge is in the pathway commons network
  if( rm.pcor$edge[i] %in% nw.edges$edge){
    rm.pcor$in_PC_KB[i] = 'Y'
  } else {
    rm.pcor$in_PC_KB[i] = 'N'
  }
}

rm.pcor <- rm.pcor %>% 
  select(symbol1 = T1, symbol2 = T2, edge, in_PC_KB, pcor) %>%
  distinct()

saveRDS(rm.pcor, paste0(here::here(), '/data/rosmap_pcor.rds'))
rm(rm.pcor)

# # Mayo --------------------------------------------------------------------
# 
# st = Sys.time()
# 
# # download and pivot pcor table
# synLogin()
# my.pcor <- read_csv(synGet('syn51061614')$path) %>% 
#   rename(gene1 = X1) %>% 
#   pivot_longer(cols = -gene1, names_to = 'gene2', values_to = 'pcor') %>% 
#   filter( gene1 != gene2 ) 
# 
# cat('download: \n')
# Sys.time() - st
# 
# # get gene symbols corresponding to Ensembl IDs
# genes <- union(my.pcor$gene1, my.pcor$gene2)
# sym <- gprofiler2::gconvert(
#   genes, 
#   organism = 'hsapiens',
#   target = 'HGNC'
# ) %>% 
#   select( gene1 = input, gene2 = input, symbol1 = target, symbol2 = target)
# 
# cat('id conversion: \n')
# Sys.time() - st
# 
# # add symbols to table
# # remove genes without symbols (therefore not in NW)
# # & put genes in consistent order to reduce dataset size
# my.pcor <- 
#   left_join( my.pcor, sym %>% select(ends_with('1')), by = 'gene1') %>% 
#   left_join( .,       sym %>% select(ends_with('2')), by = 'gene2') %>% 
#   filter(!is.na(symbol1), !is.na(symbol2)) %>% 
#   mutate(
#     T1=if_else(symbol1 < symbol2, symbol1, symbol2),
#     T2=if_else(symbol1 < symbol2, symbol2, symbol1),
#     edge = paste0(T1,':',T2), 
#     in_PC_KB = if_else(edge %in% nw.edges$edge, 'Y','N')
#   ) %>% 
#   select(symbol1 = T1, symbol2 = T2, edge, in_PC_KB, pcor) %>%
#   distinct()
# 
# saveRDS(my.pcor, paste0(here::here(), '/data/mayo_pcor.rds'))
# rm(my.pcor)
# 
# # MSSM --------------------------------------------------------------------
# 
# st = Sys.time()
# 
# # download and pivot pcor table
# synLogin()
# ms.pcor <- read_csv(synGet('syn51061324')$path) %>% 
#   rename(gene1 = X1) %>% 
#   pivot_longer(cols = -gene1, names_to = 'gene2', values_to = 'pcor') %>% 
#   filter( gene1 != gene2 ) 
# 
# cat('download: \n')
# Sys.time() - st
# 
# # get gene symbols corresponding to Ensembl IDs
# genes <- union(ms.pcor$gene1, ms.pcor$gene2)
# sym <- gprofiler2::gconvert(
#   genes, 
#   organism = 'hsapiens',
#   target = 'HGNC'
# ) %>% 
#   select( gene1 = input, gene2 = input, symbol1 = target, symbol2 = target)
# 
# cat('id conversion: \n')
# Sys.time() - st
# 
# # add symbols to table
# # remove genes without symbols (therefore not in NW)
# # & put genes in consistent order to reduce dataset size
# ms.pcor <- 
#   left_join( ms.pcor, sym %>% select(ends_with('1')), by = 'gene1') %>% 
#   left_join( .,       sym %>% select(ends_with('2')), by = 'gene2') %>% 
#   filter(!is.na(symbol1), !is.na(symbol2)) %>% 
#   mutate(
#     T1=if_else(symbol1 < symbol2, symbol1, symbol2),
#     T2=if_else(symbol1 < symbol2, symbol2, symbol1),
#     edge = paste0(T1,':',T2), 
#     in_PC_KB = if_else(edge %in% nw.edges$edge, 'Y','N')
#   ) %>% 
#   select(symbol1 = T1, symbol2 = T2, edge, in_PC_KB, pcor) %>%
#   distinct()
# 
# saveRDS(ms.pcor, paste0(here::here(), '/data/mssm_pcor.rds'))
# rm(ms.pcor)
# 
