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

# risk scores
scores <- read_csv(synTableQuery('select * from syn25575156', includeRowIdAndRowVersion = F)$filepath)

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

# base network
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
edgyboi <- args[ which( grepl('ef1|ef0', args) ) ] %>% str_remove_all('ef') %>% as.numeric()

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
if(filt_cellType){ cat('Analyze NW in: ', cellType, '\n') }
cat('Edge Factor for KDA weights?: ', edgyboi, '\n')
# cat('Skip pathway tracing and only run wKDA?: ', noTrace, '\n\n')

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

# wKDA --------------------------------------------------------------------

cat('\n\nStarting wKDA...\n')
Sys.time()

# start a directory
if( !(dir.exists( paste0(full_path,'/kda') )) ){
  dir.create( paste0(full_path,'/kda') )
}

# read traced, filtered NW objects
trace.nw.bdFilt <- igraph::read_graph(
  paste0( full_path, '/',
          'bdFiltered_', working_path,directionality, filt,
          '.graphml'),
  format = "graphml"
)

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
  distinct() %>% 
  write_tsv(paste0(full_path, '/kda/', 'module_file.tsv'))

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
  write_tsv(paste0(full_path, '/kda/', 
                   working_path, directionality, filt , 
                   '_network_file.tsv'))

edgybois <- c(0,1)
for(ef in edgyboi){
  
  ### Setup KDA job
  job.kda <- list()
  # job.kda$label<-paste0(weight_type,'_',ef)  #filename
  job.kda$label<- paste0(working_path, directionality, filt, '_addWeights_edgeFactor_',ef)  #filename
  job.kda$folder<- paste0(full_path)  #path  , '_depth2'
  job.kda$netfile <- paste0(full_path, '/kda/', 
                            working_path, directionality, filt,'_network_file.tsv')
  job.kda$modfile <- paste0(full_path, '/kda/', 'module_file.tsv')
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
