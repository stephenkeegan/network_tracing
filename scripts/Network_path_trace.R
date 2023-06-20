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

# base network
synLogin()
net <- igraph::read_graph(synGet('syn51110930')$path, format = 'graphml') 

# read arguments ----------------------------------------------------------

# arguments
# 1. query filepath
# 2. trace directed 
# 3. node filter
# 4. edge filter
# 5. cell type

# parse args and establish settings
args <- commandArgs(trailingOnly = TRUE)

cat('\narguments:\n',args, '\n')

# parse args
query.file <- args[ which( grepl('\\.txt|\\.tsv', args) ) ] 
directed <- any( grepl('dir', args) )
filt_nodes <- any( grepl('node', args) )
filt_edges <- any( grepl('edge', args) )
noTrace <- any( grepl('noTrace', args) )
filt_cellType <- any( grepl('Exc|Inh|Astro|Micro', args) )
if(filt_cellType){ cellType <- args[ which( grepl('Exc|Inh|Astro|Micro', args) ) ] }

# report args
full_path = normalizePath(query.file) %>% dirname()
working_path = normalizePath(query.file) %>% dirname() %>% basename()
cat('\n\nMove to dir: ',  full_path , '\n')
setwd(full_path)
cat('Query file to trace: ', query.file, '\n')
cat('File found in current path?: ', basename(query.file) %in% list.files(), '\n')
cat('Name of working directory: ', working_path, '\n\n')

cat('PATH TRACE OPTIONS\nTrace directed edges?: ', directed, '\n')
cat('Filter nodes for brain expression?: ', filt_nodes, '\n')
cat('Filter edges for PMID evidence?: ', filt_edges, '\n')
cat('Filter nodes for expression in certain cells?: ', filt_cellType, '\n')
if(filt_cellType){ cat('Trace NW in: ', cellType, '\n') }
# cat('Skip pathway tracing and only run wKDA?: ', noTrace, '\n\n')

# filter base network -----------------------------------------------------

cat('\n\n','Filtering Pathway Commons network based on specifications...','\n')

# nodes & edges
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

# cell type
if( filt_cellType ){
  nw <- igraph::delete_vertices(
    nw,
    V(nw)[ which( vertex_attr(nw, cellType) %in% c(0, NaN) ) ]
  )
  filt = paste0(filt,'_',cellType)
}

# directionality
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

# Read file
file_name = query.file %>% basename()

cat('\n','Tracing paths for genes in the file: ', query.file %>% basename(), 
    '\n At the path ', getwd(), ' \n\n')

query = read_tsv( query.file %>% basename(), col_names = 'gene')

# Intersect with network nodes
input.gene.list <- query$gene
if( length(setdiff( input.gene.list, names(V(nw)))) > 0 ){
  cat('\nQuery genes missing from filtered NW object: ',
      setdiff( input.gene.list, names(V(nw)) ), sep = '\n')
  cat('\n')
  
  input.gene.list = intersect( 
    query$gene, names(V(nw))
    )
}

cat('Number of queried genes to trace: ', 
    length(input.gene.list), ' / ', length(query$gene),
    ' (', length(input.gene.list)/length(query$gene)*100, '%)',
    ' \n')

# path tracing ------------------------------------------------------------

cat('\n','Beginning trace... \n')
Sys.time()

##
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

##
# Filter NW obj for traced nodes
trace_filt <- unlist(trace) %>% unique()

trace.nw <- igraph::induced_subgraph(
  nw,
  v=igraph::V(nw)[ names(igraph::V(nw)) %in% trace_filt ]
)

##
# Annotate nodes that are queried or added by trace
x = tibble( node = V(trace.nw) %>% names ) %>% 
  mutate( node_status = if_else(node %in% input.gene.list, 'query', 'added') )

igraph::vertex_attr(trace.nw, 'node_status', index = igraph::V(trace.nw)) <- x$node_status

##
# save NW trace and filtered NW
igraph::write_graph(
  trace.nw,
  paste0( full_path, '/',
          working_path, directionality, filt,
          '.graphml'),
  format = "graphml"
)

##
# Remove redundant edges
nw.simple <- igraph::simplify(
  trace.nw,
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

# calculate network stats
nw.stats <- tibble(  
  network = c('full','pre-trace','traced','simplified'),
  n_nodes = c( 
    net %>% V %>% length,
    nw %>% V %>% length,
    trace.nw %>% V %>% length,
    nw.simple %>% V %>% length
    ),
  n_edges = c(
    net %>% E %>% length,
    nw %>% E %>% length,
    trace.nw %>% E %>% length,
    nw.simple %>% E %>% length
    ),
  avg_path_length = c(
    net %>% average.path.length, 
    nw %>% average.path.length, 
    trace.nw %>% average.path.length, 
    nw.simple %>% average.path.length
    ),
  assortativity_coef = c(
    net %>% assortativity(., types1 = V(.)),
    nw %>% assortativity(., types1 = V(.)),
    trace.nw %>% assortativity(., types1 = V(.)),
    nw.simple %>% assortativity(., types1 = V(.))
    ),
  connected_components = c(
    net %>% no.clusters, 
    nw %>% no.clusters,
    trace.nw %>% no.clusters, 
    nw.simple %>% no.clusters
    )
)

write_csv(nw.stats,
          paste0( full_path, '/',
                  working_path, directionality, filt,
                  '_netStats.csv' )
          )

# plot network stats
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
  paste0( full_path, '/',
          working_path, directionality, filt,
          '_netStats.pdf' )
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
