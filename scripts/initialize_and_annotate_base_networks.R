
# setup -------------------------------------------------------------------


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


# pathway commons data ----------------------------------------------------

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


# summarise & filter PC data ----------------------------------------------

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


# >> make base NW graph ------------------------------------------------------

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
synLogin()
foo <- synStore(
  File(
    paste0(here::here(),'/data/pathway_commons_v12_undirected.graphml')
    # paste0(here::here(),'/data/pathway_commons_v13_undirected.graphml')
    , parent = 'syn51080283'
  )
)

# data for annotation -----------------------------------------------------

synLogin()

###
# TREAT-AD scores

scores <- read_csv(synTableQuery('select * from syn25575156', includeRowIdAndRowVersion = F)$filepath)
omics <- read_csv(synTableQuery('select * from syn22758536',  includeRowIdAndRowVersion = F)$filepath)

tad.deg <- scores %>% 
  filter(isScored_omics == 'Y', OmicsScore > 0) %>% 
  pull(GeneName) %>% 
  .[!is.na(.)] %>% 
  unique()

###
# Human protein atlas info

tissue.array <- read_csv(synGet('syn51074598')$path)

# Specify brain tissue types
hpa.brain.tissues <- c(
  "caudate", "cerebellum", "cerebral cortex", "hippocampus", 
  "hypothalamus", "pituitary gland", "dorsal raphe", 
  "choroid plexus", "substantia nigra"
)

# ID brain expressed genes
hpa.brain <- tissue.array %>% 
  filter(
    Tissue %in% hpa.brain.tissues, 
    Level %in% c('Low','Medium','High','Ascending','Descending'),
    Reliability %in% c('Enhanced','Supported','Approved') ) %>% 
  pull(Gene.name) %>% unique()

###
# SEA-AD scRNA-seq data

seaad <- read_csv(synGet('syn51658024')$path, guess_max = 7.5e5)
# seaad.broad <- read_csv(synGet('syn51658025')$path, guess_max = 2e5)

seaad.expr <- seaad %>%
  filter( upperQ_exp >= -2.5 * fxn_exp + 2.3 ) %>%
  pull(gene) %>% unique()

# Combined
brain.genes <- union(hpa.brain, tad.deg) %>% union(., seaad.expr)

###
# Generate cell type specific expression table

cellTypeExpr <- seaad %>%
  
  # filter out low expression
  filter( upperQ_exp >= -2.5 * fxn_exp + 2.5) %>%
  
  # # only keep genes that are detected above threshold in the same cell in multiple groups
  # # don't do this for now, only removes 757 genes including some hallmark genes (MTHFR, CTSH, etc)
  # group_by(cellType) %>% filter(duplicated(gene), .preserve = T) %>% ungroup() %>%
  
  # retain only gene and broad cell type specification
  select(gene, broad) %>% distinct() %>%
  
  # pivot the table so there are 1's and 0's to indicate if a gene is expressed in each cell class
  mutate(val = 1) %>%
  pivot_wider(
    id_cols = gene,
    names_from = broad,
    values_from = val,
    values_fill = 0
  ) %>% 
  rename(Micro = `Micro-PVM`)
  

# annotate nodes & edges --------------------------------------------------

# pathway commons graph
net <- igraph::read_graph(synGet('syn51080932')$path, format = 'graphml') 

# edge directionality
directed_edge_types = c("catalysis-precedes",
                        "controls-expression-of",
                        "controls-phosphorylation-of",
                        "controls-state-change-of",
                        "controls-transport-of"
)

x = tibble( edge = E(net)$interaction ) %>% 
  mutate(
    directed = if_else(edge %in% directed_edge_types, 1, 0)
  )

igraph::edge_attr(net, "directed") <- x$directed

# Brain expression evidence
x = tibble( node = V(net) %>% names ) %>% 
  mutate(
    hpa = if_else(node %in% hpa.brain, 1, 0)
    , ampad_deg = if_else(node %in% tad.deg, 1, 0)
    , seaad = if_else(node %in% seaad.expr, 1, 0)
    , brain_exp = if_else(node %in% brain.genes, 1, 0)
  )

for(i in 2:ncol(x)){
  igraph::vertex_attr(net, names(x)[i], index = igraph::V(net)) <- x %>% pull(i)
}

# Overall TREAT-AD Target Risk Score (TRS)
x = tibble( node = V(net) %>% names ) %>% 
  left_join(
    ., scores %>% filter(!duplicated(GeneName)) %>% select(node = GeneName, Overall) %>% distinct()
    , by = 'node'
    , na_matches = "never")
if( all( x$node != names(V(net)) ) ){ cat("Network Node IDs and annotation table IDs dont match") }

igraph::vertex_attr(net, "TargetRiskScore", index = igraph::V(net)) <- x$Overall

# Omics Effect Scores (i.e. directionality)
x = tibble( node = V(net) %>% names ) %>% 
  left_join(
    ., omics %>% filter(!duplicated(GName)) %>% select(node = GName, RNA_TE, Pro_TE) %>% distinct()
    , by = 'node'
    , na_matches = "never")
if( all( x$node != names(V(net)) ) ){ cat("Network Node IDs and annotation table IDs dont match") }

igraph::vertex_attr(net, "RNA_EffectScore", index = igraph::V(net)) <- x$RNA_TE 
igraph::vertex_attr(net, "Pro_EffectScore", index = igraph::V(net)) <- x$Pro_TE 

# Cell Type Expression
x = tibble( node = V(net) %>% names ) %>% 
  left_join(
    ., cellTypeExpr %>% filter(!duplicated(gene)) %>% rename(node = gene) %>% distinct()
    , by = 'node'
    , na_matches = "never")

for(i in 2:ncol(x)){
  igraph::vertex_attr(net, names(x)[i], index = igraph::V(net)) <- x %>% pull(i)
}

# >> save annotated NW graph -------------------------------------------------

write.graph(net, 
            file = paste0(here::here(), '/data/annotated_base_network.graphml'),
            format = 'graphml')

synLogin()
foo <- synStore(File(
  paste0(here::here(), '/data/annotated_base_network.graphml'),
  parent = 'syn51080283'))
