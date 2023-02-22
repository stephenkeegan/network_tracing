library(synapser)
library(Mergeomics)
library(igraph)
library(tidyverse)

synLogin()

scores <- read_csv(synTableQuery('select * from syn25575156', includeRowIdAndRowVersion = F)$filepath)

biodom <- full_join(
  # biodomains
  readRDS(synGet('syn25428992')$path),
  # domain labels
  read_csv(synGet('syn26856828')$path),
  by = c('Biodomain'='domain')
) %>%
  mutate(Biodomain = case_when(Biodomain == 'none' ~ NA_character_, T ~ Biodomain))

bd_enr_res = read_csv(synapser::synGet('syn45824995')$path) %>% 
  # filter(Biodomain != 'none') %>% 
  mutate(leadingEdge_genes = str_split(leadingEdge_genes, '\\|'))


kda_res_dir <- map_dfr(
  synGetChildren('syn48199344')$asList() %>% tibble(f = .) %>% unnest_wider(f) %>% pull(id),
  ~{ synGetChildren(  parent = .x  )$asList() %>% 
      tibble(f = .) %>% unnest_wider(f) %>% 
      mutate(bd = str_extract(name, '^.+(?=_[un]*directed)') %>% str_replace_all('_',' '),
             dir = str_extract(name, '[un]*directed') ,
             ef = str_extract(name, '(?<=_)e[01](?=_)')) %>% 
      relocate(bd, dir, ef)}
)

kd = map_dfr(
  kda_res_dir %>% filter(bd == 'Synapse') %>% pull(id),
  ~{f = synGetChildren(  parent = .x  )$asList() %>% 
      tibble(f = .) %>% unnest_wider(f) %>%
      filter(grepl('results.txt', name)) 
    
    read_tsv(synGet(f$id)$path, col_types = cols()) %>% 
      mutate(bd = str_extract(f$name, '^.+(?=_[un]*directed)') %>% str_replace_all('_',' '),
             dir = str_extract(f$name, '[un]*directed') ,
             ef = str_extract(f$name, '(?<=_)e[01](?=_)')) %>% 
      relocate(bd, dir, ef)}) %>% 
  mutate(FDR = -log10(FDR)) %>% 
  pivot_wider(id_cols = c(bd,dir,MODULE,NODE), names_from = ef, values_from = FDR) %>% 
  mutate(d = e1 - e0) %>% 
  arrange(desc(d)) %>% 
  left_join(., biodom %>% select(GOterm_Name, Biodomain, n_hgncSymbol, abbr, label, color), 
            by = c('MODULE'='GOterm_Name'))

kd_summary = kd %>% 
  filter(e1 > -log10(0.05)
              , MODULE %in% bd_enr_res$pathway[bd_enr_res$padj < 0.05]
              # , d > 0
         ) %>% 
  group_by(NODE) %>% 
  summarise(n_bd = length(unique(Biodomain)),
            n_term = length(unique(MODULE))) %>% 
  arrange(desc(n_term))

p2= kd %>% 
  filter( e1 > -log10(0.1)
          , MODULE %in% bd_enr_res$pathway[bd_enr_res$padj < 0.05]
          ) %>%
  ggplot(aes( d, e1 ))+
  labs(y = 'wKDA FDR ef1, -log10', x = '\u0394 FDR, ef1 - ef0', title = unique(kd$bd))+
  # geom_abline(intercept = 0, slope =1, lty= 2, lwd = .5)+
  # geom_smooth(method = 'glm', lty = 2, lwd = .5, color = 'grey20')+
  geom_vline(xintercept = 0, lwd = .5)+ geom_hline(yintercept = 0, lwd = .5) +
  geom_point(aes(color = color, text = paste0(NODE, '\n', MODULE, '\n', Biodomain)))+
  ggrepel::geom_label_repel(aes(label = paste0(NODE, '\n', MODULE, '\n', Biodomain)),
                            size = 2, alpha = .7, min.segment.length = 0)+
  facet_wrap(~dir)+ #, scales = 'free'
  scale_color_identity(guide = 'none') #; plotly::ggplotly(p)

plotly::ggplotly(p2)

