library(synapser)
synLogin()
library(tidyverse)

set.seed(113)

#### use this script to make input random gene lists for KDA
## not derived from path trace app, just pure random network
download_biodomain <- readRDS(synGet("syn25428992")$path) %>%
  dplyr::select(GOterm_Name, Biodomain, hgnc_symbol) 
biodomains <- separate_rows(download_biodomain, hgnc_symbol, sep = '","')
#118,467 obs

# make random network by GO Term
num_rows <- 1000

new_network_head <- biodomains %>% slice_sample(n = num_rows, replace = TRUE) %>%
  rename(HEAD = hgnc_symbol) %>%
  dplyr::select(HEAD)
new_network_tail <- biodomains %>% slice_sample(n = num_rows, replace = TRUE) %>%
  rename(TAIL = hgnc_symbol) %>%
  dplyr::select(TAIL)

KDA_network <- bind_cols(new_network_head, new_network_tail)

filename <- paste0('/projects/carter-lab/keegas/network_tracing/random_KDA_',num_rows,'.txt')
write.table(KDA_network, filename,
            append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


#### make random networks for pathway tracing
download_biodomain <- readRDS(synGet("syn25428992")$path) %>%
  dplyr::select(GOterm_Name, Biodomain, hgnc_symbol) 
biodomains <- separate_rows(download_biodomain, hgnc_symbol, sep = '","')

num_rows <- 1000

path_network <- biodomains %>% slice_sample(n = num_rows, replace = TRUE) %>%
  rename(HEAD = hgnc_symbol) %>%
  dplyr::select(HEAD)

filename <- paste0('/projects/carter-lab/keegas/network_tracing/random_path_input_',num_rows,'.txt')
write.table(path_network, filename,
            append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)