# This Rscript performs two broad functions:
# 1. self-trace leading edge genes from enriched biodomain term


# setup -------------------------------------------------------------------

cat('
##################################
##        PATHWAY TRACING       ##
##################################
    ')

# Package names
packages <- c('synapser','igraph','tidyverse')

# Load packages
suppressPackageStartupMessages(
  invisible(lapply(packages, library, character.only = TRUE))
)

# source path tracing functions
source(paste0(here::here(), '/scripts/igraph_NW_exp_functions.R'))

# # source Mergeomics
# source('https://raw.githubusercontent.com/jessicading/mergeomics/master/Mergeomics_Version_1.99.0.R')

theme_set(theme_bw())

cat('\npackages loaded:\n', packages, '\n')

# read data ---------------------------------------------------------------

synLogin()

# target risk scores
scores <- read_csv(synTableQuery('select * from syn25575156', 
                                 includeRowIdAndRowVersion = F)$filepath,
                   col_types = cols()
                   )

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
enr.bd <- read_csv(synGet('syn45824995')$path, col_types = cols()) %>% 
  mutate(leadingEdge_genes = str_split(leadingEdge_genes, '\\|'))

# base network
net <- igraph::read_graph(synGet('syn51110930')$path, format = 'graphml') 

# read arguments ----------------------------------------------------------

# arguments
# 1. biodomain index: 1..19
# 2. directed
# 3. node filter
# 4. edge filter
# 5. don't trace NW (only wKDA existing NW)

# parse args and establish settings
args <- commandArgs(trailingOnly = TRUE)

cat('\narguments:\n',args, '\n')

bd.idx <- args[ which( grepl('[0-9]{1,2}', args) ) ] %>% as.numeric()
directed <- any( grepl('dir', args) )
filt_nodes <- any( grepl('node', args) )
filt_edges <- any( grepl('edge', args) )
noTrace <- any( grepl('noTrace', args) )

# Specify which biodomain to trace
dom <- domains[bd.idx]
cat('\n\nBiodomain NW to trace: #', bd.idx, ', ', dom, '\n')
cat('Trace directed edges?: ', directed, '\n')
cat('Filter nodes for brain expression?: ', filt_nodes, '\n')
cat('Filter edges for PMID evidence?: ', filt_edges, '\n')
# cat('Skip pathway tracing and only run wKDA?: ', noTrace, '\n\n')

net_filename = dom %>% str_replace_all(.,' ','_')
if( !(dir.exists( paste0(here::here(), '/results/',net_filename) )) ){
  dir.create( paste0(here::here(), '/results/',net_filename) )
}

# filter base network -----------------------------------------------------

cat('\n\n','Filtering Pathway Commons network based on specifications...','\n')

if( filt_edges & filt_nodes ){
  tmp <- igraph::delete_vertices(
    net,
    igraph::V(net)[ igraph::V(net)$brain_exp == 0 ]
  )
  nw <- igraph::subgraph.edges(
    tmp,
    igraph::E(tmp)[ igraph::E(tmp)$n_edge_evidence > 1 ],
    delete.vertices = T 
  )
  filt = '_filt_node_edge'
} else if( filt_nodes ){
  nw <- igraph::delete_vertices(
    net,
    igraph::V(net)[ igraph::V(net)$brain_exp == 0 ]
  )
  filt = '_filt_node'
} else if( filt_edges ){
  nw <- igraph::subgraph.edges(
    net,
    igraph::E(net)[ igraph::E(net)$n_edge_evidence > 1 ],
    delete.vertices = T 
  )
  filt = '_filt_edge'
} else { 
  nw <- net
  filt = '_filt_none'
  }


if( directed ){
  nw <- igraph::subgraph.edges( 
    nw, 
    igraph::E(nw)[ igraph::E(nw)$directed == 1 ],
    delete.vertices = T 
  ) %>% 
    as.directed(mode = 'arbitrary')
  directionality = '_directed'  
} else { 
  directionality = '_undirected'
}

cat( 'Base network filtered.', '\n')

# gene list to trace ------------------------------------------------------

cat('\n','Tracing paths for leading edges from the <<', dom, '>> biological domain \n')

# pull gene list from leading edge genes
input.gene.list <- enr.bd %>% 
  filter(
    Biodomain == dom,
    padj < 0.01,
    NES > 1.7
  ) %>% 
  pull(leadingEdge_genes) %>% 
  unlist() %>% 
  unique()

le.genes <- input.gene.list

if( length(setdiff( input.gene.list, names(V(nw)))) > 0 ){
  cat('\nBiodomain genes missing from filtered NW object: ',
      setdiff( input.gene.list, names(V(nw)) ), sep = '\n')
  cat('\n')
  
  input.gene.list = intersect( 
    le.genes, names(V(nw))
    )
}

cat('Number of leading edge genes to trace: ', 
    length(input.gene.list), ' / ', length(le.genes),
    ' (', length(input.gene.list)/length(le.genes)*100, '%)',
    ' \n')

# path tracing ------------------------------------------------------------

cat('\n','Beginning trace... \n')
Sys.time()

# trace paths in parallel
future::plan(strategy = 'multisession', workers = 10)
trace <- furrr::future_map(
  input.gene.list,
  ~ short_paths(
    tnet = nw,
    target = .x,
    targets = input.gene.list,
    sentinals = input.gene.list,
    cores = 1)
)
future::plan(strategy = 'sequential')

# Filter NW obj for traced nodes
trace_filt <- unlist(trace) %>% unique()

trace.nw <- igraph::induced_subgraph(
  nw,
  v=igraph::V(nw)[ names(igraph::V(nw)) %in% trace_filt ]
)

# Filter NW obj for traced nodes **annotated to biodomain**
bd.genes <- biodom %>% filter(Biodomain == dom) %>% 
  pull(symbol) %>% unlist() %>% unique() %>% .[!is.na(.)]

trace_bd_filt <- unlist(trace) %>% intersect(bd.genes,.)

trace.nw.bdFilt <- igraph::induced_subgraph(
  nw, 
  v=igraph::V(nw)[ names(igraph::V(nw)) %in% trace_bd_filt ]
)

# save NW trace and filtered NW
igraph::write_graph(
  trace.nw,
  paste0( here::here(), '/results/',net_filename,'/',
          net_filename, directionality, filt,'.graphml'),
  format = "graphml"
)

igraph::write_graph(
  trace.nw.bdFilt,
  paste0( here::here(), '/results/',net_filename,'/',
          'bdFiltered_', net_filename, directionality, filt , '.graphml'),
  format = "graphml"
)

##
# Remove redundant edges
nw.simple <- igraph::simplify(
  trace.nw.bdFilt,
  remove.multiple = TRUE,
  remove.loops = FALSE,
  edge.attr.comb = list( 
    interaction = 'concat',
    edge = 'random',
    occurrance = 'concat',
    n_edge = 'max',
    n_edge_types = 'max',
    n_edge_evidence = 'max',
    n_source = 'max',
    sources = 'concat',
    n_evidence = 'sum',
    evidence_pmid = 'concat',
    n_pathways = 'max',
    pathway_names = 'concat',
    directed = 'max'
  )
)

cat('\n','Trace complete and networks saved. \n')
Sys.time()

# generate plots ----------------------------------------------------------

cat('\n\nGenerating plots...\n')

# Network Stats
nw.stats <- tibble(  
  network = c('full','pre-trace','traced','biodom_filtered','simplified'),
  n_nodes = c( 
    net %>% V %>% length,
    nw %>% V %>% length,
    trace.nw %>% V %>% length,
    trace.nw.bdFilt %>% V %>% length,
    nw.simple %>% V %>% length
    ),
  n_edges = c(
    net %>% E %>% length,
    nw %>% E %>% length,
    trace.nw %>% E %>% length,
    trace.nw.bdFilt %>% E %>% length,
    nw.simple %>% E %>% length
    ),
  avg_path_length = c(
    net %>% average.path.length, 
    nw %>% average.path.length, 
    trace.nw %>% average.path.length, 
    trace.nw.bdFilt %>% average.path.length, 
    nw.simple %>% average.path.length
    ),
  assortativity_coef = c(
    net %>% assortativity(., types1 = V(.)),
    nw %>% assortativity(., types1 = V(.)),
    trace.nw %>% assortativity(., types1 = V(.)),
    trace.nw.bdFilt %>% assortativity(., types1 = V(.)),
    nw.simple %>% assortativity(., types1 = V(.))
    ),
  connected_components = c(
    net %>% no.clusters, 
    nw %>% no.clusters,
    trace.nw %>% no.clusters, 
    trace.nw.bdFilt %>% no.clusters, 
    nw.simple %>% no.clusters
    )
)

write_csv(nw.stats,
          paste0( here::here(), '/results/',net_filename,'/',
                  net_filename, directionality, filt, '_netStats.csv' )
          )

nw.stats %>% 
  pivot_longer(cols = -network, names_to = 'properties', values_to = 'val') %>% 
  mutate(properties = factor(properties, 
                             levels = c('n_nodes',
                                        'n_edges', 
                                        'avg_path_length',
                                        'assortativity_coef',
                                        'connected_components')) 
         , network = factor(network, levels = c('full',
                                                'pre-trace',
                                                'traced',
                                                'biodom_filtered',
                                                'simplified'))) %>% 
  ggplot(aes(network, val)) +
  geom_bar(stat = 'identity', position = 'dodge')+
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  labs(y = '',x = '', subtitle = 'base network properties')+
  facet_wrap(~properties, scales = 'free_y', ncol = 2)

ggsave(
  paste0( here::here(), '/results/',net_filename,'/',
          net_filename, directionality, filt, '_netStats.pdf' )
)

# # synapse upload ----------------------------------------------------------
# 
# parent_id <- 'syn51117833'
# 
# d <- list.dirs('results', recursive = F) %>% 
#   str_subset('input_gene_lists|wKDA', negate=T)
# 
# for(j in 4:5){
#   foo <- synStore( Folder(d[j] %>% str_remove_all('results/'), parent = parent_id) )
#   f <- list.files(d[j]) %>% str_subset('kda',negate=T)
#   for(i in 1:length(f)){
#     foo2 <- synStore( File(
#       paste0(d[j],'/',f[i]),
#       parent=foo$properties$id
#     ))
#   }
#   foo3 <- synStore( Folder('kda', parent = foo$properties$id) )
#   f <- list.files( paste0(d[j],'/','kda')) 
#   for(i in 1:length(f)){
#     foo4 <- synStore( File(
#       paste0(d[j],'/kda/',f[i]),
#       parent=foo3$properties$id
#     ))
#   }
# }

# EOF ##