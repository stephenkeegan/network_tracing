# setup -------------------------------------------------------------------
library(synapser)
library(igraph)
library(tidyverse)

synLogin()
directed_edges <- T

# Base networks -----------------------------------------------------------

if(directed_edges == T){
  net <- read_graph(paste0(here::here(),'/data/path_commons_directed.graphml'), format = 'graphml')
} else {
  net <- read_graph(paste0(here::here(),'/data/path_commons_undirected.graphml'), format = 'graphml')
  }

# Filter edges
test_net <-  igraph::subgraph.edges( 
  graph = net,
  eids = igraph::E(net)[ ( 
    igraph::E(net)$n_edge == 1  &
    igraph::E(net)$n_edge_evidence < 3 
    # TODO: other filters? co-expression? 
    ) == F ],
  delete.vertices = TRUE
)
net <- test_net

# Annotate vertecies
omics <- read_csv( synTableQuery('select * from syn22758536', includeRowIdAndRowVersion = F)$filepath )
scores <- read_csv( synTableQuery('select * from syn25575156', includeRowIdAndRowVersion = F)$filepath )

idx = match( names(V(net)), omics$GName ) 
igraph::vertex_attr(net, "RNA_EffectScore", index = igraph::V(net)) <- omics$RNA_TE[idx]
igraph::vertex_attr(net, "Pro_EffectScore", index = igraph::V(net)) <- omics$Pro_TE[idx]

idx = match( names(V(net)), scores$GeneName ) 
# igraph::vertex_attr(net, "target_risk_score", index = igraph::V(net)) <- scores$Overall[idx]
igraph::vertex_attr(net, "weight", index = igraph::V(net)) <- scores$Overall[idx]

# Biodomain & enrichment results ------------------------------------------

biodom <- full_join(
  # biodomains
  readRDS(synGet('syn25428992')$path),
  # domain labels
  read_csv(synGet('syn26856828')$path),
  by = c('Biodomain'='domain')
) %>%
  mutate(Biodomain = case_when(Biodomain == 'none' ~ NA_character_, T ~ Biodomain))

enr_biodom <- read_csv(synGet('syn45824995')$path) %>% 
  # TODO: other filteres here? more stringent? filter on NES > median(NES) of sig terms
  filter(padj < 0.05, Biodomain != 'none') %>%
  mutate(leadingEdge_genes = str_split(leadingEdge_genes, '\\|'))

leading_edge <- enr_biodom %>% 
  group_by(Biodomain) %>% 
  summarise(le = leadingEdge_genes %>% unlist() %>% unique() %>% list() ) %>% 
  pull(le, name = Biodomain)

bds = purrr::map_dfr(1:length(leading_edge), 
                     ~ tidyr::tibble( 
                       biodom = names(leading_edge)[.x], 
                       n_leadingEdge = length(leading_edge[[.x]]))) %>% 
  dplyr::arrange(n_leadingEdge)
bds

# trace -------------------------------------------------------------------

bd.targets = leading_edge[['APP Metabolism']]

future::plan(strategy = 'multisession', workers = 10)
tictoc::tic()
# trace <- furrr::future_map(
#   bd.targets, 
#   ~ igraphNetworkExpansion::short_paths(
#       tnet = net,
#       target = .x,
#       targets = bd.targets,
#       sentinals = bd.targets,
#       cores = 1)
# )

trace <- map(
  bd.targets, 
  ~ igraphNetworkExpansion::short_paths(
    tnet = net,
    target = .x,
    targets = bd.targets,
    sentinals = bd.targets,
    cores = 1)
)
tictoc::toc()
future::plan(strategy = 'sequential')

targ <- NULL
sent <- NULL
for( i in 1:length(trace) ){
  targ <- c(targ,trace[[i]]$Inter)
  sent <- c(sent,trace[[i]]$Sentinal)
}

u_targ <- targ[!duplicated(targ)]
u_sent <- sent[!duplicated(sent)]

# Only genes found in both target and sentinal traces
opt_1 <- u_targ[u_targ%in%u_sent]

# All genes present more than twice
opt_2 <- c(
  names(table(sent)[table(sent)>1]),
  names(table(targ)[table(targ)>1])
)
opt_2<-opt_2[!duplicated(opt_2)]

# All Genes 
opt_3 <- c(
  names(table(sent)),
  names(table(targ))
)
opt_3<-opt_3[!duplicated(opt_3)]

# only biodom gene

#trace_filt <- igraphNetworkExpansion::trace_filter(trace)
trace_filt <- opt_3

nw <- igraph::induced_subgraph(
  net, v=igraph::V(net)[ names(igraph::V(net)) %in% trace_filt ])


# subnet_simple <- igraph::simplify(
#   nw,
#   remove.multiple = TRUE,
#   remove.loops = FALSE,
#   edge.attr.comb = list( interaction = "concat", 
#                          Occurance = "concat",
#                          UniqCol = "concat",
#                          pathway = "concat", 
#                          EdgeRep = "mean",
#                          Edge = "random",
#                          SumOccurence = "mean",
#                          DLPFC_CE = "mean",
#                          CBE_CE = "mean",
#                          FP_CE = "mean",
#                          IFG_CE = "mean",
#                          PHG_CE = "mean",
#                          STG_CE = "mean",
#                          TCX_CE = "mean",
#                          Avg_Cortex_CE = "mean",
#                          Avg_All_CE = "mean"
#   )
# )

# store and upload --------------------------------------------------------
if(directed_edges == F){
  directionality = '_undirected'
} else { 
  directionality = '_directed'  }

net_filename = bds$biodom[idx] %>% str_replace_all(.,' ','_')
igraph::write_graph(
  nw,
  paste0( net_filename, directionality,'.graphml'),
  format = "graphml"
)

# Store output on Synapse
# foo = Folder(name = 'Greg_sandbox', parent = 'syn25190666')
# f = synStore(foo)
# dir = f$properties$id
# dir = 'syn29884372'

if(Directed == 'NO'){
  dir = 'syn45831456' # undirected
} else { 
  dir = 'syn45831445' # directed  
}

net_synname = paste0(bds$biodom[idx],' Leading Edge Self Trace')
obj <- synapser::synStore(
  synapser::File(paste0( net_filename,directionality,'.graphml'), name = net_synname, parent = dir),
  used = syns_used,
  #executed = prov,
  activityName = paste0(bds$biodom[idx], ' Subnetwork'),
  activityDescription = paste0('Tracces Leading edge genes of significcant GO-Terms of the', 
                               bds$biodom[idx],' Biodomain to themselves')
  
) 



