# library(paxtoolsr)
# library(RJSONIO)
# library(httr)
# library(DOSE)
# library(org.Hs.eg.db)
# library(clusterProfiler)
# library(parallel)
# library(doParallel)
# library(RColorBrewer)
#devtools::install_github("jhrcook/HotNetvieweR")
# library(HotNetvieweR)
# source('igraphhack/igraphplot2.R')

library(synapser)
library(igraph)
library(tidyverse)

synLogin()

pc.sources <- read_tsv('https://www.pathwaycommons.org/archives/PC2/v12/datasources.txt', skip = 2,
                       col_names = c("ID","DESCRIPTION","TYPE","HOMEPAGE","PATHWAYS","INTERACTIONS","PARTICIPANTS")) %>% 
  mutate(source = str_remove_all(ID, 'http://pathwaycommons.org/pc12/')) %>% 
  relocate(source)

pc.pathways <- read_tsv( synGet('syn47046495')$path , guess_max = )

pc.data.dir <- synGetChildren(
  # 'syn46947105' # v13 sif versions
  # 'syn47045122' # v13 txt versions
  'syn21907759' # v12 txt versions
  )$asList() %>% 
  tibble(f = .) %>% unnest_wider(f) %>% 
  # mutate(source = str_extract(name, '(?<=PathwayCommons13\\.).+(?=\\.hgnc)')) %>% 
  mutate(source = str_extract(name, '(?<=PathwayCommons12\\.).+(?=\\.hgnc)')) %>% 
  relocate(source) %>% 
  filter( 
    !is.na(source)
    # (source %in% c('All', 'Detailed'))==F 
    )

# all <- read_tsv(synGet(pc.data.dir$id[1])$path)
# detailed <- read_tsv(synGet(pc.data.dir$id[2])$path)
 
pc.data <- map_dfr(
  3:nrow(pc.data.dir),
  ~ read_tsv(synGet(pc.data.dir$id[.x])$path, col_types = 'ccccccc') 
)

# pc.data.sif <- map_dfr(
#   3:nrow(pc.data.dir),
#   ~ read_tsv(synGet(pc.data.dir$id[.x])$path, col_names = c('from','interaction','to')) %>%
#     mutate(source = pc.data.dir$source[.x]) %>% relocate(source)
# )

pc.filt = pc.data %>% 
  filter( 
    !grepl('CHEBI:', PARTICIPANT_A), !grepl('CHEBI:', PARTICIPANT_B),
    !grepl('DnaReference|DnaRegionReference|ProteinReference|RnaReference|RnaRegionReference|PARTICIPANT_TYPE', INTERACTION_TYPE)
    ) %>% 
  rename(from = PARTICIPANT_A, to = PARTICIPANT_B, interaction = INTERACTION_TYPE) %>%
  mutate(
    T1=ifelse(from < to, from, to),
    T2=ifelse(from < to, to, from),
    edge = paste0(T1,':',T2)
  ) %>% 
  select(-T1,-T2) %>% 
  group_by(edge) %>% 
  mutate(
    occurrance = paste0(from,'-',to,':', interaction),
    n_edge = n(),
    n_edge_types = length(unique(interaction)),
    n_edge_evidence = str_split(INTERACTION_PUBMED_ID[!is.na(INTERACTION_PUBMED_ID)], ';') %>% unlist() %>% unique() %>% length(),
    ) %>% 
  ungroup() %>% 
  group_by(from, to, interaction
           , edge, occurrance
           , n_edge, n_edge_types, n_edge_evidence
           ) %>% 
  summarise(
  # mutate(
    # n_edge = length(edge), 
    n_source = length(unique(INTERACTION_DATA_SOURCE)),
    sources = paste0(unique(INTERACTION_DATA_SOURCE), collapse = ','),
    n_evidence = str_split(INTERACTION_PUBMED_ID[!is.na(INTERACTION_PUBMED_ID)], ';') %>% unlist() %>% unique() %>% length(),
    evidence_pmid = paste0(unique(INTERACTION_PUBMED_ID[!is.na(INTERACTION_PUBMED_ID)]), collapse = '|'),
    n_pathways = str_split(PATHWAY_NAMES[!is.na(PATHWAY_NAMES)], ';') %>% unlist() %>% unique() %>% length(),
    pathway_names = paste0( unique(PATHWAY_NAMES[!is.na(PATHWAY_NAMES)]), collapse = '|' )
  ) 

# # PC interaction type specifications:
# # https://github.com/BioPAX/Paxtools/blob/2f93afa94426bf8b5afc2e0e61cd4b269a83288d/pattern/src/main/resources/org/biopax/paxtools/pattern/miner/sif-guide.html
# directed_edge_types = c("catalysis-precedes",
#                         "controls-expression-of",
#                         "controls-phosphorylation-of",
#                         "controls-state-change-of", 
#                         "controls-transport-of")

# make graphs
undirected = graph_from_data_frame(d = pc.filt, directed = F)
# directed = graph_from_data_frame(d = pc.filt %>% filter(interaction %in% directed_edge_types), directed = T)

# output base network graphs
write_graph(undirected, paste0(here::here(),'/data/pathway_commons_v12_undirected.graphml'), format = 'graphml')
# write_graph(undirected, paste0(here::here(),'/data/pathway_commons_v13_undirected.graphml'), format = 'graphml')
# write_graph(directed, paste0(here::here(),'/data/pathway_commons_v13_directed.graphml'), format = 'graphml')

# upload to synapse
foo <- synStore(
  File(
    paste0(here::here(),'/data/pathway_commons_v12_undirected.graphml')
    # paste0(here::here(),'/data/pathway_commons_v13_undirected.graphml')
    , parent = 'syn51080283'
  )
)
