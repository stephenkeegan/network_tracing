library(synapser)
synLogin()
library(tidyverse)

set.seed(113)

#### use this script to make input random gene lists for KDA
## not derived from path trace app, just pure random network
download_biodomain <- readRDS(synGet("syn25428992")$path) %>%
  dplyr::select(GOterm_Name, Biodomain, hgnc_symbol) 
biodomains <- separate_rows(download_biodomain, hgnc_symbol, sep = '","')

# make random network by GO Term
num_rows <- 10

new_network_head <- biodomains %>% slice_sample(n = num_rows, replace = TRUE) %>%
  rename(HEAD = hgnc_symbol) %>%
  dplyr::select(HEAD)
new_network_tail <- biodomains %>% slice_sample(n = num_rows, replace = TRUE) %>%
  rename(TAIL = hgnc_symbol) %>%
  dplyr::select(TAIL)

new_network <- bind_cols(new_network_head, new_network_tail)

