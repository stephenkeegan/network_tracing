# This Rscript performs two broad functions:
# 1. self-trace leading edge genes from enriched biodomain terms
# 2. perform wKDA on resulting NWs

# setup -------------------------------------------------------------------

# Package names
packages <- c('synapser','igraph','tidyverse')

# Load packages
suppressPackageStartupMessages(
  invisible(lapply(packages, library, character.only = TRUE))
)

# source path tracing functions
source(paste0(here::here(), '/scripts/igraph_NW_exp_functions.R'))

# source Mergeomics
source('https://raw.githubusercontent.com/jessicading/mergeomics/master/Mergeomics_Version_1.99.0.R')

theme_set(theme_bw())

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

cat('\narguments:\n',args, sep='\n')

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
cat('Skip pathway tracing and only run wKDA?: ', noTrace, '\n\n')

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
    net, 
    igraph::E(net)[ igraph::E(net)$directed == 1 ],
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

net_filename = dom %>% str_replace_all(.,' ','_')
if( !(dir.exists( paste0(here::here(), '/results/',net_filename) )) ){
  dir.create( paste0(here::here(), '/results/',net_filename) )
}

if( noTrace ){
  
  trace.nw <- read.graph(
    paste0( here::here(), '/results/',net_filename,'/',
            net_filename, filt, directionality,'.graphml'),
    format = 'graphml'
  )
  trace <- list(names(V(trace.nw)))
  
} else {
  
  cat('\n','Beginning trace... \n')
  Sys.time()

  # trace paths in parallel
  future::plan(strategy = 'multisession', workers = 10)
  # tictoc::tic()
  # system.time({
    trace <- furrr::future_map(
      input.gene.list,
      ~ short_paths(
        tnet = nw,
        target = .x,
        targets = input.gene.list,
        sentinals = input.gene.list,
        cores = 1)
    )
    # tictoc::toc()
  # })
  future::plan(strategy = 'sequential')

  # Filter NW obj for traced nodes
  trace_filt <- unlist(trace) %>% unique()

  trace.nw <- igraph::induced_subgraph(
    nw,
    v=igraph::V(nw)[ names(igraph::V(nw)) %in% trace_filt ]
  )
  
}

# Filter NW obj for traced nodes **annotated to biodomain**
bd.genes <- biodom %>% filter(Biodomain == dom) %>% 
  pull(hgnc_symbol) %>% unlist() %>% unique() %>% .[!is.na(.)]

trace_bd_filt <- unlist(trace) %>% intersect(bd.genes,.)

trace.nw.bdFilt <- igraph::induced_subgraph(
  nw, 
  v=igraph::V(nw)[ names(igraph::V(nw)) %in% trace_bd_filt ]
)

# save NW trace and filtered NW
igraph::write_graph(
  trace.nw,
  paste0( here::here(), '/results/',net_filename,'/',
          net_filename, filt, directionality,'.graphml'),
  format = "graphml"
)

igraph::write_graph(
  trace.nw.bdFilt,
  paste0( here::here(), '/results/',net_filename,'/',
          'bdFiltered_', net_filename, filt , directionality,'.graphml'),
  format = "graphml"
)

cat('\n','Trace complete and networks saved. \n')
Sys.time()

# wKDA --------------------------------------------------------------------

cat('\n\nStarting wKDA...\n')
Sys.time()

if( !exists('trace.nw.bdFilt') ){
  trace.nw.bdFilt <- read.graph(
    paste0( here::here(), '/results/',net_filename,'/',
            'bdFiltered_', net_filename, filt , directionality,'.graphml' ),
    format = 'graphml'
  )
}

if( !(dir.exists( paste0(here::here(), '/results/',net_filename,'/kda') )) ){
  dir.create( paste0(here::here(), '/results/',net_filename,'/kda') )
}

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

##
# Add TREAT-AD scores to edge attributes
kda.nw <- tibble( ea = edge.attributes(nw.simple) ) %>% 
  t() %>% as_tibble(rownames = NA, .name_repair = 'unique') %>% 
  unnest(everything()) %>% 
  rename_with(., ~names(edge.attributes(nw.simple)), everything()) %>% 
  select(edge, interaction, occurrance, directed, n_edge, 
         n_edge_evidence, n_source, n_edge_types) %>% 
  mutate(HEAD = str_split_fixed(edge, ':',2)[,1], 
         TAIL = str_split_fixed(edge, ':',2)[,2]) %>% 
  relocate(HEAD, TAIL) %>% select(-edge) %>% 
  left_join(., 
            scores %>% 
              select(HEAD = GeneName, h.c = Overall, h.g = GeneticsScore, 
                     h.o = OmicsScore, h.n = NeuropathScore),
            by = 'HEAD') %>% 
  left_join(., 
            scores %>% 
              select(TAIL = GeneName, t.c = Overall, t.g = GeneticsScore, 
                     t.o = OmicsScore, t.n = NeuropathScore),
            by = 'TAIL')

##
# Module file
biodom %>% 
  filter(!is.na(n_hgncSymbol)
         # , Biodomain == dom
  ) %>% 
  select(MODULE = GOterm_Name, NODE = hgnc_symbol) %>% 
  unnest_longer(NODE) %>% 
  write_tsv(paste0(here::here(), '/results/', net_filename, '/kda/',
                   'module_file.tsv'))

##
# generate network file
kda.nw %>% 
  mutate(
    across(.cols = c(starts_with('h.'), starts_with('t.')), 
           ~ if_else(is.na(.x), 0, .x)),
    WEIGHT = h.c+t.c
  ) %>%
  select(HEAD, TAIL, WEIGHT) %>% 
  distinct() %>% 
  write_tsv(paste0(here::here(), '/results/', net_filename, '/kda/',
                   net_filename, filt , directionality,'network_file.tsv'))

edgybois <- c(0,1)

for(ef in edgybois){
  
  ### Setup KDA job
  job.kda <- list()
  # job.kda$label<-paste0(weight_type,'_',ef)  #filename
  job.kda$label<- paste0(filt, directionality,'_addWeights_edgeFactor_',ef)  #filename
  job.kda$folder<- paste0(here::here(), '/results/', net_filename,'/')  #path  , '_depth2'
  job.kda$netfile <- paste0(here::here(), '/results/', net_filename, '/kda/', 
                            net_filename, filt , directionality,'network_file.tsv')
  job.kda$modfile <- paste0(here::here(), '/results/', net_filename, '/kda/',
                            'module_file.tsv')
  job.kda$edgefactor<- ef  #edgybois
  job.kda$depth<- 1
  if(directed) {job.kda$direction <- 1} else {job.kda$direction <- 0}
  job.kda$nperm <- 100000  #100000  
  moddata <- tool.read(job.kda$modfile)
  ## save this to a temporary file and set its path as new job.kda$modfile:
  tool.save(moddata, "subsetof.supersets.txt")
  job.kda$modfile <- "subsetof.supersets.txt"
  
  ## Running KDA
  job.kda <- kda.configure(job.kda)
  job.kda <- kda.start(job.kda)
  job.kda <- kda.prepare(job.kda)
  job.kda <- kda.analyze(job.kda)
  job.kda <- kda.finish(job.kda)
  
}

cat('\n\nwKDA complete. \n')
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
                  net_filename, filt, directionality, '_netStats.csv' )
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
          net_filename, filt, directionality, '_netStats.pdf' )
)

# wKDA results
kda.res  <- bind_rows(
  read_tsv(paste0(here::here(),'/results/', net_filename,
                  '/kda/',filt, directionality,'_addWeights_edgeFactor_0.results.txt')) %>% 
    mutate(ef = 'ef0_fdr'),
  read_tsv(paste0(here::here(),'/results/', net_filename,
                  '/kda/',filt, directionality,'_addWeights_edgeFactor_1.results.txt')) %>% 
    mutate(ef = 'ef1_fdr') 
  ) %>%
  mutate(FDR = -log10(FDR)) %>%
  pivot_wider(id_cols = c(MODULE, NODE), names_from = ef, values_from = FDR) %>%
  mutate(delta_fdr = ef1_fdr - ef0_fdr) %>% arrange(desc(delta_fdr))

write_csv(kda.res,
          paste0(here::here(),'/results/', net_filename,'/kda/',
                 filt, directionality,'_delta_edgeFactor.txt')
          )

kda.res %>% 
  inner_join(., biodom %>% select(GOterm_Name, Biodomain, abbr, label, color, n_hgncSymbol),
             by = c('MODULE'='GOterm_Name')) %>%
  filter( ef1_fdr > -log10(0.1) ) %>%
  ggplot(aes( ef1_fdr, delta_fdr ))+
  labs(x = 'wKDA FDR ef1, -log10', y = 'delta FDR, ef1 - ef0')+
  # geom_abline(intercept = 0, slope =1, lty= 2, lwd = .5)+
  geom_smooth(method = 'glm', lty = 2, lwd = .5, color = 'grey20')+
  geom_vline(xintercept = 0, lwd = .5)+ geom_hline(yintercept = 0, lwd = .5) +
  geom_point(aes(color = color, text = paste0(NODE, '\n', MODULE, '\n', Biodomain)))+
  ggrepel::geom_label_repel(aes(label = paste0(NODE, '\n', MODULE, '\n', Biodomain)),
                            size = 2, alpha = .7, min.segment.length = 0)+
  scale_color_identity(guide = 'none')

ggsave(
  paste0( here::here(), '/results/',net_filename,'/kda/',filt, directionality,'_deltaEF.pdf' )
)

# EOF ##