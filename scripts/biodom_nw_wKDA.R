# This Rscript performs two broad functions:
# 2. perform wKDA on resulting NWs

# setup -------------------------------------------------------------------

cat('
##################################
## WEIGHTED KEY DRIVER ANALYSIS ##
##################################
    ')

# Package names
packages <- c('synapser','igraph','tidyverse')

# Load packages
suppressPackageStartupMessages(
  invisible(lapply(packages, library, character.only = TRUE))
)

# # source path tracing functions
# source(paste0(here::here(), '/scripts/igraph_NW_exp_functions.R'))

# source Mergeomics
source('https://raw.githubusercontent.com/jessicading/mergeomics/master/Mergeomics_Version_1.99.0.R')

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
# 6. wKDA edge factor

# parse args and establish settings
args <- commandArgs(trailingOnly = TRUE)

cat('\narguments:\n',args, sep='\n')

bd.idx <- args[ which( grepl('[0-9]{1,2}', args) ) ] %>% str_subset('ef',negate = T) %>% as.numeric()
directed <- any( grepl('dir', args) )
filt_nodes <- any( grepl('node', args) )
filt_edges <- any( grepl('edge', args) )
noTrace <- any( grepl('noTrace', args) )
edgyboi <- args[ which( grepl('ef1|ef0', args) ) ] %>% str_remove_all('ef') %>% as.numeric()

# Specify which biodomain to trace
dom <- domains[bd.idx]
cat('\n\nBiodomain NW to trace: #', bd.idx, ', ', dom, '\n')
cat('Trace directed edges?: ', directed, '\n')
cat('Filter nodes for brain expression?: ', filt_nodes, '\n')
cat('Filter edges for PMID evidence?: ', filt_edges, '\n')
cat('Edge Factor for KDA weights?: ', edgyboi, '\n')
# cat('Skip pathway tracing and only run wKDA?: ', noTrace, '\n\n')

net_filename = dom %>% str_replace_all(.,' ','_')

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

# read traced NW objects
trace.nw <- igraph::read_graph(
  paste0( here::here(), '/results/',net_filename,'/',
          net_filename, directionality, filt, '.graphml'),
  format = "graphml"
)
trace.nw.bdFilt <- igraph::read_graph(
  paste0( here::here(), '/results/',net_filename,'/',
          'bdFiltered_', net_filename, directionality, filt , '.graphml'),
  format = "graphml"
)

# wKDA --------------------------------------------------------------------

cat('\n\nStarting wKDA...\n')
Sys.time()

if( !(dir.exists( paste0(here::here(), '/results/',net_filename,'/kda') )) ){
  dir.create( paste0(here::here(), '/results/',net_filename,'/kda') )
}

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
  filter(!is.na(n_symbol),
         n_symbol > 0
         # , Biodomain == dom
  ) %>% 
  select(MODULE = GOterm_Name, NODE = symbol) %>% 
  unnest_longer(NODE) %>% 
  filter(NODE != '') %>% 
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
                   net_filename, directionality, filt , 'network_file.tsv'))

edgybois <- c(0,1)

for(ef in edgyboi){
  
  ### Setup KDA job
  job.kda <- list()
  # job.kda$label<-paste0(weight_type,'_',ef)  #filename
  job.kda$label<- paste0(net_filename, directionality, filt, '_addWeights_edgeFactor_',ef)  #filename
  job.kda$folder<- paste0(here::here(), '/results/', net_filename,'/')  #path  , '_depth2'
  job.kda$netfile <- paste0(here::here(), '/results/', net_filename, '/kda/', 
                            net_filename, directionality, filt,'network_file.tsv')
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

# # generate plots ----------------------------------------------------------
# 
# cat('\n\nGenerating plots...\n')
# 
# # wKDA results
# kda.res  <- bind_rows(
#   read_tsv(paste0(here::here(),'/results/', net_filename, '/kda/',
#                   net_filename, directionality, filt,'_addWeights_edgeFactor_0.results.txt')) %>% 
#     mutate(ef = 'ef0_fdr'),
#   read_tsv(paste0(here::here(),'/results/', net_filename,'/kda/',
#                   net_filename, directionality, filt, '_addWeights_edgeFactor_1.results.txt')) %>% 
#     mutate(ef = 'ef1_fdr') 
#   ) %>%
#   mutate(FDR = -log10(FDR)) %>%
#   pivot_wider(id_cols = c(MODULE, NODE, MEMBER), names_from = ef, values_from = FDR) %>%
#   mutate(
#     across(.cols = ends_with('fdr'), ~if_else(is.na(.x), 0, .x)),
#     delta_fdr = ef1_fdr - ef0_fdr) %>% arrange(desc(delta_fdr))
# 
# write_csv(kda.res,
#           paste0(here::here(),'/results/', net_filename,'/kda/',
#                  net_filename, directionality, filt, '_delta_edgeFactor.txt')
#           )
# 
# n_kd = kda.res %>% filter(ef1_fdr > -log10(0.1)) %>% nrow()
# 
# kda.res %>% 
#   inner_join(., biodom %>% select(GOterm_Name, Biodomain, abbr, label, color, n_symbol),
#              by = c('MODULE'='GOterm_Name')) %>%
#   filter( if(n_kd > 0) ef1_fdr > -log10(0.1) else !is.na(ef1_fdr)) %>%
#   ggplot(aes( ef1_fdr, delta_fdr ))+
#   labs(x = 'wKDA FDR ef1, -log10', y = 'delta FDR, ef1 - ef0')+
#   # geom_abline(intercept = 0, slope =1, lty= 2, lwd = .5)+
#   geom_smooth(method = 'glm', lty = 2, lwd = .5, color = 'grey20')+
#   geom_vline(xintercept = 0, lwd = .5)+ geom_hline(yintercept = 0, lwd = .5) +
#   geom_point(aes(color = color, text = paste0(NODE, '\n', MODULE, '\n', Biodomain)))+
#   ggrepel::geom_label_repel(aes(label = paste0(NODE, '\n', MODULE, '\n', Biodomain)),
#                             size = 2, alpha = .7, min.segment.length = 0)+
#   scale_color_identity(guide = 'none')
# 
# ggsave(
#   paste0( here::here(), '/results/',net_filename,'/kda/',net_filename, directionality, filt, '_deltaEF.pdf' )
# )

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