# initialize biodomain NW trace queries and directories

# Package names
packages <- c('synapser','tidyverse')

# Load packages
suppressPackageStartupMessages(
  invisible(lapply(packages, library, character.only = TRUE))
)

synLogin()

# biological domain annotations
biodom <- full_join(
  # biodomains
  readRDS(synGet('syn25428992')$path),
  # domain labels
  read_csv(synGet('syn26856828')$path,
           col_types = cols()),
  by = c('Biodomain'='domain')
) %>%
  mutate(Biodomain = case_when(Biodomain == 'none' ~ NA_character_, T ~ Biodomain))

domains <- biodom %>% pull(Biodomain) %>% unique() %>% sort() %>% .[!is.na(.)]

# enriched biodomain terms
trs.enr <- read_csv(synGet('syn45824995')$path, col_types = cols()) %>% 
  mutate(leadingEdge_genes = str_split(leadingEdge_genes, '\\|'))
gen.enr <- read_csv(synGet('syn45824969')$path, col_types = cols()) %>% 
  mutate(leadingEdge_genes = str_split(leadingEdge_genes, '\\|'))
omic.enr <- read_csv(synGet('syn45824835')$path, col_types = cols()) %>% 
  mutate(leadingEdge_genes = str_split(leadingEdge_genes, '\\|'))

enr.bd <-  trs.enr

# initialize directories and leading edge gene list files
setwd('/projects/carter-lab/caryg/network_tracing/results/biodomain')
for( bd in domains ){
  
  # pull gene list from leading edge genes
  gl <- enr.bd %>% 
    filter(
      Biodomain == bd,
      padj < 0.01,
      NES > 1.7
    ) %>% 
    pull(leadingEdge_genes) %>% 
    unlist() %>% 
    unique()

  cat(paste0(bd, ': ', length(gl), ' leading edge genes\n'))
  
  if(length(gl) == 0){next}
  
  bd_filename <- bd %>% str_replace_all(.,' ','_')
  if(!dir.exists(bd_filename)){ dir.create(bd_filename) }
  
  # write query list to file
  write_tsv(tibble(x = gl), paste0(bd_filename, '/queryList_', bd_filename, '.tsv'), col_names = F)
  
}
