library(synapser)
library(igraph)
library(ggraph)
library(tidyverse)


# read data ---------------------------------------------------------------

synLogin()

nw = read_graph(synGet('syn51121574')$path, format = 'graphml')
kda_res = read_csv( synGet('syn51121650')$path, col_types = cols() ) %>% 
  rownames_to_column(var = 'idx')

kda_res %>% arrange(desc(delta_fdr)) %>% filter(!duplicated(MODULE)) %>% View()

# specify KD/term to plot -------------------------------------------------

idx = 542
kd = kda_res$NODE[idx] #e.g. 'HDAC1'
term = kda_res$MODULE[idx] #e.g. 'positive regulation of transcription by RNA polymerase II'

vert_num <- match(kd, V(nw)$name)
vert_neigh <- neighbors(nw, vert_num)
term_num <- biodom %>% 
  filter(GOterm_Name == term) %>% 
  pull(hgnc_symbol) %>% unlist() %>% 
  match(., V(nw)$name) %>% .[!is.na(.)]

g2 <- induced_subgraph(nw, c(vert_num,intersect(term_num, vert_neigh)) )
g3 <- induced_subgraph(nw, c(vert_num,term_num) )

# plot -- circular layout -------------------------------------------------

cowplot::plot_grid(
  
  ggraph(g3, layout = 'linear', circular = TRUE) + 
    geom_edge_arc(alpha = .3, color = 'grey75') + 
    geom_node_point(
      shape = 21, color = "black", alpha = .3, 
      aes(size = TargetRiskScore, fill = Pro_EffectScore, filter = .N()$name != kd))+
    geom_node_point(
      shape = 23, color = "black", alpha = .3, 
      aes(size = TargetRiskScore, fill = Pro_EffectScore, filter = .N()$name == kd))+
    geom_node_text(aes(label = name), size = 2) + 
    coord_fixed()+
    scale_size_continuous(range = c(.5,10), limits = c(0,5),trans = 'exp')+
    scale_fill_gradient2(low = 'purple',mid = 'black',high= 'green', na.value = 'grey85',
                          limits = range(omics$Pro_TE, na.rm =T))+
    labs(title = 'All term genes in NW', 
         subtitle = paste0('KD: ',kd,'\nTerm: ',term)       )+
    theme_graph(background = 'white' )
  ,
  ggraph(g2, layout = 'linear', circular = TRUE) + 
    geom_edge_arc(alpha = .3, color = 'grey75') + 
    geom_node_point(
      shape = 21, color = "black", alpha = .3, 
      aes(size = TargetRiskScore, fill = Pro_EffectScore, filter = .N()$name != kd))+
    geom_node_point(
      shape = 23, color = "black", alpha = .3, 
      aes(size = TargetRiskScore, fill = Pro_EffectScore, filter = .N()$name == kd))+
    geom_node_text(aes(label = name), size = 2) + 
    coord_fixed()+
    scale_size_continuous(range = c(.5,10), limits = c(0,5), trans = 'exp')+
    scale_fill_gradient2(low = 'purple',mid = 'black',high= 'green', na.value = 'grey85',
                         limits = range(omics$Pro_TE, na.rm =T))+
    labs(title = 'Overlap of KD neighbors & term genes ', 
         subtitle = paste0('KD: ',kd,'\nTerm: ',term)       )+
    theme_graph(background = 'white' )
  
  , ncol = 2
  
)

# plot -- kk layout -------------------------------------------------------

ggraph(g3, layout = 'kk') + 
  geom_edge_link(alpha = .3, color = 'grey85') + 
  geom_node_point(
    shape = 21, color = "black", alpha = .3, 
    aes(size = TargetRiskScore, fill = Pro_EffectScore, filter = .N()$name != kd))+
  geom_node_point(
    shape = 23, color = "black", alpha = .3, 
    aes(size = TargetRiskScore, fill = Pro_EffectScore, filter = .N()$name == kd))+
  geom_node_text(aes(label = name), size = 2) + 
  scale_fill_gradient2(low = 'purple',mid = 'black',high= 'green', na.value = 'grey85',
                       limits = range(omics$Pro_TE, na.rm =T))+
  scale_size_continuous(range = c(.5,10), limits = c(0,5), trans = 'exp')+
  labs(title = 'All term genes in NW', 
       subtitle = paste0('KD: ',kd,'\nTerm: ',term)        )+
  theme_graph(background = 'white' )

ggraph(g2, layout = 'kk') + 
  geom_edge_link(alpha = .3, color = 'grey85') + 
  geom_node_point(
    shape = 21, color = "black", alpha = .3, 
    aes(size = TargetRiskScore, fill = Pro_EffectScore, filter = .N()$name != kd))+
  geom_node_point(
    shape = 23, color = "black", alpha = .3, 
    aes(size = TargetRiskScore, fill = Pro_EffectScore, filter = .N()$name == kd))+
  geom_node_text(aes(label = name), size = 2) + 
  scale_fill_gradient2(low = 'purple',mid = 'black',high= 'green', na.value = 'grey85',
                       limits = range(omics$Pro_TE, na.rm =T))+
  scale_size_continuous(range = c(.5,10), limits = c(0,5), trans = 'exp')+
  labs(title = 'Overlap of KD neighbors & term genes ',
       subtitle = paste0('KD: ',kd,'\nTerm: ',term)        )+
  theme_graph(background = 'white' )


# plot -- other -----------------------------------------------------------


# ggraph(g3, 'focus', focus = node_is_center()) + 
#   ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = r), data.frame(r = 1:5), colour = 'grey') + 
#   geom_edge_link(alpha = .3, color = 'grey75') + 
#   geom_node_point(aes(size = TargetRiskScore^2, color = Pro_EffectScore), alpha = .5)+
#   geom_node_text(aes(label = name), size = 2,  position = position_dodge(width = .5) ) + 
#   coord_fixed()

# ggraph(g3, 'fabric', sort.by = node_rank_fabric()) + 
#   geom_node_range(colour = 'grey') + 
#   geom_edge_span(end_shape = 'square') + 
#   coord_fixed()

# ggraph(g3, 'fabric', sort.by = node_rank_fabric(), shadow.edges =TRUE) + 
#   geom_node_range(colour = 'grey') + 
#   geom_edge_span(aes(filter = shadow_edge), colour ='lightblue' , end_shape = 'square') + 
#   geom_edge_span(aes(filter = !shadow_edge), end_shape = 'square') + 
#   coord_fixed()
