library(igraph)
library(tidyverse)

# base net filtering ------------------------------------------------------

tmp = nw.dat %>% 
  select(filt) %>% distinct() %>% 
  rowwise() %>% 
  mutate(
    n_nodes = case_when( filt == 'none' ~ length(names(igraph::V(undir))),
                         filt == 'node' ~ length(names(igraph::V(undir_node))),
                         filt == 'edge' ~ length(names(igraph::V(undir_edge))),
                         filt == 'node_edge' ~ length(names(igraph::V(undir_node_edge))) ),
    n_edges = case_when( filt == 'none' ~ length((igraph::E(undir))),
                         filt == 'node' ~ length((igraph::E(undir_node))),
                         filt == 'edge' ~ length((igraph::E(undir_edge))),
                         filt == 'node_edge' ~ length((igraph::E(undir_node_edge))) )
  ) %>% 
  pivot_longer(cols = c(n_nodes, n_edges), 
               names_to = 'properties', 
               values_to = 'val') %>% 
  mutate(filt = factor(filt, 
                       levels = c('node_edge', 'edge','node','none')),
         properties = factor(properties, 
                             levels = c('n_nodes',
                                        'n_edges')) 
  ) 

p1 = tmp %>% mutate(val = as.double(val)) %>% 
  ggplot(aes(filt, val)) +
  geom_bar(stat = 'identity', position = 'dodge')+
  labs(x = 'network filter',y = '')+
  # scale_y_log10()+
  scale_y_continuous( labels = scales::scientific)+
  coord_flip()+ theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  facet_wrap(~properties, scales = 'free_x', ncol = 2)


# LE gene overlap ---------------------------------------------------------

olap1 = map_dfr(
  unique('none'),
  ~{
    tmp <- nw.dat %>% 
      filter(dir == 'undirected', 
             filt == .x)
    
    olap <- crossing(nw1 = tmp$bd, nw2 = tmp$bd) %>% 
      group_by(nw1,nw2) %>% 
      mutate(.keep = 'all', 
             intersect = length(intersect(
               tmp %>% filter(bd == nw1) %>% pull(le_genes) %>% unlist(),
               tmp %>% filter(bd == nw2) %>% pull(le_genes) %>% unlist()
             ))
      ) %>% 
      ungroup() %>% group_by(nw1) %>% 
      mutate(.keep = 'all', pct_overlap = intersect / intersect[which(nw2 == nw1)]) %>% 
      ungroup() %>% 
      mutate(filt = .x)
  }
) %>%
  mutate(
    # filt = fct_relevel(filt, c('none','node','edge','node_edge')),
    nw1 = fct_relevel(
      nw1, 
      nw.dat %>% select(bd, le_genes) %>% distinct() %>% 
        mutate(len = length(unlist(le_genes))) %>% arrange(len) %>% pull(bd)
    ),
    nw2 = fct_relevel(
      nw2, 
      nw.dat %>% select(bd, le_genes) %>% distinct() %>% 
        mutate(len = length(unlist(le_genes))) %>% arrange(desc(len)) %>% pull(bd)
    ),
    set = 'Leading Edge Genes'
  ) 


p2 = ggplot(olap1, aes(x = nw2, y = nw1)) +
  geom_tile(aes(fill=pct_overlap), alpha = .8) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  viridis::scale_fill_viridis("% Total", option = "C", direction = -1)+
  geom_text(aes(label = paste0(intersect,'\n', signif(100*pct_overlap,3), '%') ), size = 2.5) +
  facet_wrap(~set)+
  labs(x='',y='') 


# input vs added scores ---------------------------------------------------

tmp <- nw.dat %>% 
  select(bd, filt, dir, input_genes, added_genes) %>% 
  pivot_longer(cols = c(input_genes, added_genes)) %>% 
  unnest(value) %>% 
  left_join(
    .,
    scores %>% select(value = GeneName, Overall, 
                      GeneticsScore, OmicsScore, LiteratureScore),
    by = 'value', na_matches = 'never'
  )

p3 = tmp %>% 
  filter(filt =='node_edge', dir == 'undirected',
         bd %in% c('Endolysosome','Immune Response',
                   'Lipid Metabolism','Mitochondrial Metabolism')) %>% 
  rename(set = name, gene = value) %>% 
  select(bd, set, gene, value = Overall) %>% 
  filter(value > 0) %>% 
  distinct() %>% 
  ggplot(aes(value))+
  theme(legend.position = 'top') + 
  # geom_density(aes(fill = set, color = set), alpha = .3)+
  geom_histogram(aes(fill = set), alpha = .8, position='dodge')+
  scale_fill_discrete('Gene Set')+ #,breaks = c('Input Genes','Genes Added By Trace')
  labs(x = 'TREAT-AD Target Risk Score')+
  facet_wrap(~bd, scales = 'free')


# lap of nodes added by trace ---------------------------------------------

olap2 = map_dfr(
  # unique(nw.dat$filt),
  c('none'),
  ~{
    tmp <- nw.dat %>% 
      filter(dir == 'undirected', 
             filt == .x)
    
    olap <- crossing(nw1 = tmp$bd, nw2 = tmp$bd) %>% 
      group_by(nw1,nw2) %>% 
      mutate(.keep = 'all', 
             intersect = length(intersect(
               setdiff( tmp %>% filter(bd == nw1) %>% pull(nodes) %>% unlist(), 
                        tmp %>% filter(bd == nw1) %>% pull(le_genes) %>% unlist() ),
               setdiff( tmp %>% filter(bd == nw2) %>% pull(nodes) %>% unlist(), 
                        tmp %>% filter(bd == nw2) %>% pull(le_genes) %>% unlist() )
             ))
      ) %>% 
      ungroup() %>% group_by(nw1) %>% 
      mutate(.keep = 'all', pct_overlap = intersect / intersect[which(nw2 == nw1)]) %>% 
      ungroup() %>% 
      mutate(filt = .x)
  }
)

ord = olap1 %>% filter(nw1==nw2) %>% arrange(desc(intersect)) %>% pull(nw1)

olap <- olap2 %>% 
  mutate(set = 'Nodes Added by Trace') %>% 
  bind_rows(olap1,.) %>% 
  mutate(
    nw1 = factor(nw1, levels = rev(ord) ),
    nw2 = factor(nw2, levels = (ord))
    )

p4 = olap %>% 
  # filter(filt == 'node_edge') %>% 
  ggplot(., aes(x = nw2, y = nw1)) +
  geom_tile(aes(fill=pct_overlap), alpha = .8) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  # theme(axis.text.x = element_blank(), legend.position = 'bottom') +
  # scale_fill_gradient2('% Total', low = 'white', high = '4a66ac', na.value = 'grey90') +
  viridis::scale_fill_viridis("% Total", option = "C", direction = -1)+
  # scale_x_discrete(position = "top") +
  # geom_text(aes(label = paste0(intersect, '\n', signif(pct_overlap*100, 3), '%') ), size = 1) +
  geom_text(aes(label = intersect ), size = 2.5) +
  labs(x='',y='') + 
  facet_wrap(~set)


# enrichment of NW nodes --------------------------------------------------

tmp <- gl %>% 
  filter(name == 'nodes',
         filt == 'node_edge',
         dir == 'directed', 
         Biodomain != 'none')

bdt <- bd.tally(tmp %>% filter(p_value <= 0.05) %>% pull(term_name), biodom)

tmp = tmp %>% 
  filter(bd %in% c('APP Metabolism','Endolysosome','Immune Response','Mitochondrial Metabolism'))%>% 
  rowwise() %>% 
  mutate(bd = paste0(bd,'\n', length(unique(unlist(value))), ' total nodes') ) %>% 
  mutate(Biodomain = factor(Biodomain, 
                            levels = bdt$domain %>% 
                              str_subset('none', negate = T) %>%
                              c('none',.))
  ) 

p5 = tmp  %>% 
  ggplot(aes(-log10(p_value), Biodomain))+
  geom_violin(data = subset(tmp, p_value <= 0.05), aes(color = color), scale = 'width')+
  geom_jitter(data = subset(tmp, p_value <= 0.05), aes(color = color), alpha = .7)+
  scale_color_identity()+
  # xlim(c(0,30))+
  labs(y = '')+
  facet_wrap(~bd, scales = 'free_x', ncol = 4)

# figure ------------------------------------------------------------------

cowplot::plot_grid(
  cowplot::plot_grid(
    cowplot::plot_grid(p1,p3,nrow = 2, rel_heights = c(.6,1) ),
    cowplot::plot_grid(NULL,p4,NULL,nrow = 3, rel_heights = c(.1,1,.1))
    , nrow =1, rel_widths = c(.5,1)),
  p5,nrow = 2, rel_heights = c(1, .7))


# kd comp -----------------------------------------------------------------

tmp <- bind_rows(
  am.kd.res, en.kd.res, ir.kd.res
) %>% 
  filter(ef1_fdr > -log10(0.05)
         # , MODULE %in% enr.bd$pathway[enr.bd$padj < 0.05]
         , delta_fdr > 0
         , top_kd == 1
         , filt == 'none'
         , dir == 'undirected'
         , !is.na(Biodomain)
  ) %>%
  mutate(pair = paste0(NODE,':',MODULE)) %>% 
  group_by(
    # filt, dir, 
          bd, Biodomain) %>% 
  summarise(n_term = length(unique(MODULE)),
            n_node = length(unique(NODE)),
            n_kd = length(unique(pair)),
            # terms = list(unique(MODULE)),
            min_fdr = signif( 10^-max(ef1_fdr), digits = 3),
            max_delta = signif( max(delta_fdr), digits = 3) ) %>% 
  full_join(., 
            biodom %>% group_by(Biodomain) %>% mutate(n_go = length(unique(GO_ID))) %>% 
              select(Biodomain, n_go, abbr:color) %>% distinct(), 
            by = 'Biodomain') %>% 
  ungroup() %>%
  mutate(fxn_bd_term = n_node / n_go) %>% 
  filter(!is.na(bd), !is.na(Biodomain)) %>%
  mutate(#dir = fct_relevel(dir, c('undirected','directed')),
    # filt = fct_relevel(filt, c('none','node','edge','node_edge')),
    Biodomain = fct_reorder(Biodomain, n_node)) %>% 
  arrange(Biodomain, desc(n_node), desc(max_delta)) 

tmp %>% 
  ggplot(aes(n_node, Biodomain))+
  geom_segment(aes(yend=Biodomain, xend=0), color = 'grey50')+
  geom_point(aes(size = max_delta, color = color))+
  scale_color_identity()+
  labs(y='', x = '# key drivers')+ 
  theme(legend.position = 'bottom')+
  scale_size_continuous('max \u0394 FDR', range = c(.05,5))+
  facet_wrap(~bd, ncol = 3, scales = 'free_x')

