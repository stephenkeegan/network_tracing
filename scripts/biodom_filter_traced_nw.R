# filter biodomain path traced networks

# Package names
packages <- c('synapser','igraph','tidyverse')

# Load packages
suppressPackageStartupMessages(
  invisible(lapply(packages, library, character.only = TRUE))
)


# load data ---------------------------------------------------------------

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
enr.bd <- read_csv(synGet('syn45824995')$path, col_types = cols()) %>% 
  mutate(leadingEdge_genes = str_split(leadingEdge_genes, '\\|'))

# read args ---------------------------------------------------------------

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

full_path = normalizePath(query.file) %>% dirname()
working_path = normalizePath(query.file) %>% dirname() %>% basename()

# set filters -------------------------------------------------------------

# nodes & edges
if( filt_edges & filt_nodes ){
  filt = '_filt_node_edge'
} else if( filt_nodes ){
  filt = '_filt_node'
} else if( filt_edges ){
  filt = '_filt_edge'
} else { 
  filt = '_filt_none'
}

# cell type
if( filt_cellType ){
  filt = paste0(filt,'_',cellType)
}

# directionality
if( directed ){
  directionality = '_directed'  
} else { 
  directionality = '_undirected'
}

# filter traced NW ----------------------------------------------------------

cat('Filter traced NW for nodes annotated to biodomain. \n')

# read nw
trace.nw <- igraph::read_graph(
  paste0( full_path, '/',
          working_path, directionality, filt,
          '.graphml')
  , format = 'graphml') 

dom = working_path %>% str_replace_all(.,'_',' ')

# Filter NW obj for traced nodes **annotated to biodomain**
bd_genes <- biodom %>% 
  filter(Biodomain == dom) %>%
  pull(symbol) %>% 
  unlist() %>% 
  unique() %>% 
  .[!is.na(.)]

trace_nw_genes <- names(igraph::V(trace.nw))
trace_nw_filter <- intersect(trace_nw_genes, bd_genes)

trace.nw.bdFilt <- igraph::induced_subgraph(
  trace.nw,
  v=igraph::V(trace.nw)[ names(igraph::V(trace.nw)) %in% trace_nw_filter ]
)

# write filtered nw
igraph::write_graph(
  trace.nw.bdFilt,
  paste0( full_path, '/',
          'bdFiltered_', working_path,directionality, filt,
          '.graphml'),
  format = "graphml"
)

cat('Filtered NW saved to:',
    paste0( full_path, '/',
            'bdFiltered_', working_path,directionality, filt,
            '.graphml')
    )
