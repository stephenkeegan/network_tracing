# setup -------------------------------------------------------------------
library(synapser)
library(igraph)
library(tidyverse)

synLogin()

 
cat(args)

directed_edges <- F
if( args[1] == 'directed' ){
  directed_edges <- T
}

input_path = args[2]
file_name = input_path %>% 
  str_remove_all(., 'results/input_gene_lists/')

# Base networks -----------------------------------------------------------
 
# if(directed_edges == T){
#   net <- read_graph(paste0(here::here(),'/data/path_commons_directed.graphml'), format = 'graphml')
# } else {
#   net <- read_graph(paste0(here::here(),'/data/path_commons_undirected.graphml'), format = 'graphml')
#   }
# 
# # Filter edges
# test_net <-  igraph::subgraph.edges( 
#   graph = net,
#   eids = igraph::E(net)[ ( 
#     igraph::E(net)$n_edge == 1  &
#     igraph::E(net)$n_edge_evidence < 3 
#     # TODO: other filters? co-expression? 
#     ) == F ],
#   delete.vertices = TRUE
# )
# net <- test_net
# 
# # Annotate vertecies
# omics <- read_csv( synTableQuery('select * from syn22758536', includeRowIdAndRowVersion = F)$filepath )
# scores <- read_csv( synTableQuery('select * from syn25575156', includeRowIdAndRowVersion = F)$filepath )
# 
# idx = match( names(V(net)), omics$GName ) 
# igraph::vertex_attr(net, "RNA_EffectScore", index = igraph::V(net)) <- omics$RNA_TE[idx]
# igraph::vertex_attr(net, "Pro_EffectScore", index = igraph::V(net)) <- omics$Pro_TE[idx]
# 
# idx = match( names(V(net)), scores$GeneName ) 
# # igraph::vertex_attr(net, "target_risk_score", index = igraph::V(net)) <- scores$Overall[idx]
# igraph::vertex_attr(net, "weight", index = igraph::V(net)) <- scores$Overall[idx]


#Pull Networks:
load( synapser::synGet('syn22992753')$path )
  # load( synapser::synGet('syn22992709')$path )
  #Dummy a network (Could choose JS or Non-JS here)
net <- net_undirected

# Filter Edges
  # lose ( 3691 vertices ) -- A == 10176 Vertacies and 151912 Edges
test_net <-  igraph::subgraph.edges( 
  net,
  igraph::E(net)[ ( igraph::E(net)$Occurance == 1  &
                      igraph::E(net)$Avg_Cortex_CE == 0 ) == F ],
  delete.vertices = TRUE
)
net <- test_net

# prune NW for only directed edges
directed_edge_types = c("catalysis-precedes",
                        "controls-expression-of",
                        "controls-phosphorylation-of",
                        "controls-state-change-of", 
                        "controls-transport-of")

if( directed_edges == T ){
  test_net <- igraph::subgraph.edges( net, 
                                      igraph::E(net)[ igraph::E(net)$interaction %in% 
                                                        directed_edge_types ],
                                      delete.vertices = T
  )
  net <- test_net
}


# Annotate base network ---------------------------------------------------

#Annotate vertices on Omics Weights:
OMICS_dep <- read.csv(synapser::synTableQuery('select * from syn22758536')$filepath)
OMICS <- read.csv( synapser::synTableQuery( 
  paste0( 'SELECT * FROM syn25575156 WHERE GeneName in (\'',
          paste( names(igraph::V(net)), collapse = '\',\'' ),
          '\')'),
  resultsAs = 'csv' )$filepath)
OMICS <-  OMICS[ , c('ENSG', 'GeneName', 'OmicsScore', 'GeneticsScore', 'Overall')]
colnames(OMICS)[ colnames(OMICS) == 'GeneName' ] <- 'GName'
# colnames(OMICS)[ colnames(OMICS) == 'Logsdon' ] <- 'Overall'

OMICS$GName <- as.character( OMICS$GName )
OMICS$ENSG <- as.character(OMICS$ENSG)
OMICS_dep$GName <- as.character( OMICS_dep$GName )
OMICS_dep$ENSG <- as.character(OMICS_dep$ENSG)

OMICS_alt <- OMICS[ (OMICS$GName %in% "") == F, ]
OMICS_dep <- OMICS_dep[ (OMICS_dep$GName %in% "") == F, ]

#Pull out pseduo genes and NAs,  also ENSG00000281123 is a diff ENSG for RNA and Protein...:
OMICS_alt <- OMICS_alt[ (OMICS_alt$ENSG %in% c(
  'ENSG00000272655',
  'ENSG00000284770',
  'ENSG00000168255',
  'ENSG00000281123'
)) == F,]
OMICS_dep <- OMICS_dep[ (OMICS_dep$ENSG %in% c(
  'ENSG00000272655',
  'ENSG00000284770',
  'ENSG00000168255',
  'ENSG00000281123'
)) == F,]

OMICS_alt <- OMICS_alt[ is.na(OMICS_alt$GName)==F,]
OMICS_alt <- OMICS_alt[ OMICS_alt$ENSG %in% as.character(OMICS_dep$ENSG), ] 

OMICS_alt  <- OMICS_alt[ !duplicated(OMICS_alt$ENSG), ]
OMICS_alt <- OMICS_alt[!duplicated(OMICS_alt$GName),]
OMICS_alt <- OMICS_alt[ is.na(OMICS_alt$GName) == F, ]
row.names( OMICS_alt ) <- OMICS_alt$GName

OMICS_dep <- OMICS_dep[ is.na(OMICS_dep$GName) == F, ]
OMICS_dep <- OMICS_dep[ !duplicated(OMICS_dep$GName), ]
row.names( OMICS_dep ) <- OMICS_dep$GName

OMICS_alt$RNA_TE <- OMICS_dep[ row.names(OMICS_alt), ]$RNA_TE
OMICS_alt$Pro_TE <- OMICS_dep[ row.names(OMICS_alt), ]$Pro_TE

#vertex_attr(net, "weight", index = V(net)) <- OMICS_alt[ names( V(net)), ]$Final_Weight 
OMICS_alt[ names( igraph::V(net)), ]$Overall[ is.na(OMICS_alt[ names( igraph::V(net)), ]$Overall) ] <- 0
igraph::vertex_attr(net, "weight", index = igraph::V(net)) <- OMICS_alt[ names( igraph::V(net)), ]$Overall
igraph::vertex_attr(net, "RNA_EffectScore", index = igraph::V(net)) <- OMICS_alt[ names( igraph::V(net)), ]$RNA_TE 
igraph::vertex_attr(net, "Pro_EffectScore", index = igraph::V(net)) <- OMICS_alt[ names( igraph::V(net)), ]$Pro_TE 

# Zero out the TEs
OMICS_alt$PRO_TE_Cor <- OMICS_alt$Pro_TE
OMICS_alt$RNA_TE_Cor <- OMICS_alt$RNA_TE

igraph::vertex_attr(net, "RNA_Cor_EffectScore", index = igraph::V(net)) <- OMICS_alt[ names( igraph::V(net)), ]$RNA_TE_Cor
igraph::vertex_attr(net, "Pro_Cor_EffectScore", index = igraph::V(net)) <- OMICS_alt[ names( igraph::V(net)), ]$PRO_TE_Cor 

# trace -------------------------------------------------------------------

query <- read_tsv( paste0(here::here(), '/results/input_gene_lists/', file_name), col_names = F) %>% 
  filter( X1 %in% names(igraph::V(net)) ) %>% 
  pull( X1 )

future::plan(strategy = 'multisession', workers = 10)
tictoc::tic()

trace <- map(
  query, 
  ~ igraphNetworkExpansion::short_paths(
    tnet = net,
    target = .x,
    targets = query,
    sentinals = query,
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


# All Genes 
opt_3 <- c(
  names(table(sent)),
  names(table(targ))
)
opt_3<-opt_3[!duplicated(opt_3)]

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

# net_filename = bds$biodom[idx] %>% str_replace_all(.,' ','_')
net_filename = file_name %>% str_remove_all(., '_leadingEdge.tsv')
igraph::write_graph(
  nw,
  paste0(here::here(), '/results/NW/', net_filename, directionality,'.graphml'),
  format = "graphml"
)

# Store output on Synapse
# foo = Folder(name = 'Greg_sandbox', parent = 'syn25190666')
# f = synStore(foo)
# dir = f$properties$id
# dir = 'syn29884372'

# if(Directed == 'NO'){
#   dir = 'syn45831456' # undirected
# } else { 
#   dir = 'syn45831445' # directed  
# }
# 
# net_synname = paste0(bds$biodom[idx],' Leading Edge Self Trace')
# obj <- synapser::synStore(
#   synapser::File(paste0( net_filename,directionality,'.graphml'), name = net_synname, parent = dir),
#   used = syns_used,
#   #executed = prov,
#   activityName = paste0(bds$biodom[idx], ' Subnetwork'),
#   activityDescription = paste0('Tracces Leading edge genes of significcant GO-Terms of the', 
#                                bds$biodom[idx],' Biodomain to themselves')
#   
# ) 

# EOF