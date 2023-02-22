# # https://github.com/jessicading/mergeomics
# source('https://raw.githubusercontent.com/jessicading/mergeomics/master/Mergeomics_Version_1.99.0.R')

# setup -------------------------------------------------------------------

# .libPaths('/projects/carter-lab/caryg/rlib/')

types = c('flat','add','mult','comp')

library(synapser)
# library(Mergeomics)
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

args = commandArgs(trailingOnly=TRUE)
args = 7 # APP Metab; 1 # Immune; 4 # SS; 14 # Mito; 15 # Synapse
cat(args)
cat('\n')
dom = biodom %>% pull(Biodomain) %>% unique() %>% 
  # str_subset(., 'Cell Cycle|Metal Binding|DNA Repair', negate = T) %>% 
  .[!is.na(.)] %>% .[as.numeric(args[1])]
cat(dom)
cat('\n')

directed = T
if(directed) {directionality = '_directed'} else {directionality = ''}

# dom = 'APP Metabolism'
# dom = 'Tau Homeostasis'
# dom = 'Endolysosome'
# dom = 'Lipid Metabolism'
# dom = 'Structural Stabilization'
# dom = 'Immune Response'
# dom = 'Synapse'
# dom = 'Mitochondrial Metabolism'

##
# point to directory on synapse containing NW files
# syn45831428: full traced biodom networks
# syn48170345: filtered for biodom annotation
top.dir <- synGetChildren('syn48170345')$asList() %>% 
  tibble(f = .) %>% unnest_wider(f) 
nw.dir <- {if(directed){ 
  synGetChildren(  parent = top.dir$id[top.dir$name == 'directed']  )$asList() 
  } else { 
    synGetChildren(  parent = top.dir$id[top.dir$name == 'undirected']  )$asList() 
    }} %>% 
  tibble(f = .) %>% unnest_wider(f) %>% 
  mutate(bd = str_extract(name, '^.+(?= Leading)')) %>% 
  relocate(bd)

nw.id <- nw.dir %>% filter(bd == dom) %>% pull(id)
bd.nw <- read_graph(synGet(nw.id)$path, format = 'graphml')

node_attr <- tibble( na = vertex.attributes(bd.nw) ) %>%
  t() %>% as_tibble(rownames = NA) %>% unnest(everything()) %>%
  rename_with(., ~names(vertex.attributes(bd.nw)), everything())

edge_attr <- tibble( ea = edge.attributes(bd.nw) ) %>% 
  t() %>% as_tibble(rownames = NA) %>% unnest(everything()) %>% 
  rename_with(., ~names(edge.attributes(bd.nw)), everything())

kda.nw <- edge_attr %>% select(Edge, interaction, Occurance, EdgeRep, SumOccurence, starts_with('Avg_')) %>% 
  mutate(HEAD = str_split_fixed(Edge, ':',2)[,1], 
         TAIL = str_split_fixed(Edge, ':',2)[,2]) %>% 
  relocate(HEAD, TAIL) %>% select(-Edge) %>% 
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



# Explore Weighting Options -----------------------------------------------

# x = kda.nw %>%
#   mutate(
#     across(.cols = c(starts_with('h.'), starts_with('t.')), ~ if_else(is.na(.x), 0, .x)),
#     wA = h.c + t.c,
#     w1 = h.g * t.g ,
#     w2 = h.o * t.o,
#     w3 = h.c * t.c,
#     w3.1 = w1+w2,
#     w4 = w1/max(w1) + w2/max(w2) + Avg_All_CE/max(Avg_All_CE, na.rm =T) + SumOccurence/max(SumOccurence) ) %>%
#   relocate(wA, w1,w2,w3,w3.1,w4) %>% arrange(desc(wA))

# p <- kda.nw %>% mutate(w1 = h.c + t.c , w2 = (h.g * t.g)+(h.o * t.o)) %>%
#   ggplot(aes(w1,w2, text = paste0(HEAD, ' > ', TAIL))) + geom_point() +
#   geom_abline(intercept = 0, slope = 1, lty =2); 
# plotly::ggplotly(p)

# prep wKDA ---------------------------------------------------------------

##
# generate module file
biodom %>% 
  filter(!is.na(n_hgncSymbol)
         # , Biodomain == dom
  ) %>% 
  select(MODULE = GOterm_Name, NODE = hgnc_symbol) %>% 
  unnest_longer(NODE) %>% 
  write_tsv('module_file.tsv')

for(weight_type in types){
  # for(weight_type in c('comp')){
  
  ##
  # generate network file
  kda.nw %>% 
    mutate(
      across(.cols = c(starts_with('h.'), starts_with('t.')), ~ if_else(is.na(.x), 0, .x)),
      across(.cols = c(starts_with('Avg')), ~ if_else(is.na(.x), 0, .x)),
      WEIGHT = case_when(weight_type == 'flat' ~ 1,
                         weight_type == 'add' ~ h.c+t.c,
                         weight_type == 'mult' ~ (h.c+1)*(t.c+1),
                         weight_type == 'comp' ~ (h.g*t.g)/max(h.g*t.g, na.rm = T) +
                           (h.o*t.o)/max(h.o*t.o, na.rm = T) +
                           SumOccurence / max(SumOccurence, na.rm = T) 
                         # + Avg_All_CE/max(Avg_All_CE, na.rm =T)
                         )
    ) %>%
    select(HEAD, TAIL, WEIGHT) %>% 
    distinct() %>% 
    write_tsv('network_file.tsv')
  
  edgybois <- c(0,1)
  
  for(ef in edgybois){
    
    ### Setup KDA job
    job.kda <- list()
    job.kda$label<-paste0(weight_type,'_',ef)  #filename
    job.kda$folder<- paste0(here::here(), '/results/wKDA/', dom, directionality, '_depth2')  #path
    job.kda$netfile <- 'network_file.tsv'
    job.kda$modfile <- 'module_file.tsv'
    job.kda$edgefactor<- ef  #edgybois
    job.kda$depth<- 2  # 1
    if(directed) {job.kda$direction <- 1} else {job.kda$direction <- 0}
    job.kda$nperm <- 100000  #100000  
    moddata <- tool.read(job.kda$modfile)
    ## save this to a temporary file and set its path as new job.kda$modfile:
    tool.save(moddata, "subsetof.supersets.txt")
    job.kda$modfile <- "subsetof.supersets.txt"
    
    # run wKDA ----------------------------------------------------------------
    
    ## Running KDA
    job.kda <- kda.configure(job.kda)
    job.kda <- kda.start(job.kda)
    job.kda <- kda.prepare(job.kda)
    job.kda <- kda.analyze(job.kda)
    job.kda <- kda.finish(job.kda)
    
  }}

# read and analyze results ------------------------------------------------

dom = biodom %>% pull(Biodomain) %>% unique() %>% .[!is.na(.)] %>%
  # str_subset(., 'Cell Cycle|Metal Binding|DNA Repair', negate = T) %>%
  str_subset(., 'APP')

directed = T
if(directed) {directionality = '_directed'} else {directionality = ''}

flat <- bind_rows(
  read_tsv(paste0(here::here(),'/results/wKDA/',dom,directionality,'/kda/flat_0.results.txt')) %>% mutate(ef = 'e0'),
  read_tsv(paste0(here::here(),'/results/wKDA/',dom,directionality,'/kda/flat_1.results.txt')) %>% mutate(ef = 'e1') ) %>%
  mutate(FDR = -log10(FDR)) %>%
  pivot_wider(id_cols = c(MODULE, NODE), names_from = ef, values_from = FDR) %>%
  mutate(d = e1-e0) %>% arrange(desc(d))

add  <- bind_rows(
  read_tsv(paste0(here::here(),'/results/wKDA/',dom,directionality,'/kda/add_0.results.txt')) %>% mutate(ef = 'e0'),
  read_tsv(paste0(here::here(),'/results/wKDA/',dom,directionality,'/kda/add_1.results.txt')) %>% mutate(ef = 'e1') ) %>%
  mutate(FDR = -log10(FDR)) %>%
  pivot_wider(id_cols = c(MODULE, NODE), names_from = ef, values_from = FDR) %>%
  mutate(d = e1-e0) %>% arrange(desc(d))

mult  <- bind_rows(
  read_tsv(paste0(here::here(),'/results/wKDA/',dom,directionality,'/kda/mult_0.results.txt')) %>% mutate(ef = 'e0'),
  read_tsv(paste0(here::here(),'/results/wKDA/',dom,directionality,'/kda/mult_1.results.txt')) %>% mutate(ef = 'e1') ) %>%
  mutate(FDR = -log10(FDR)) %>%
  pivot_wider(id_cols = c(MODULE, NODE), names_from = ef, values_from = FDR) %>%
  mutate(d = e1-e0) %>% arrange(desc(d))

comp <- bind_rows(
  read_tsv(paste0(here::here(),'/results/wKDA/',dom,directionality,'/kda/comp_0.results.txt')) %>% mutate(ef = 'e0'),
  read_tsv(paste0(here::here(),'/results/wKDA/',dom,directionality,'/kda/comp_1.results.txt')) %>% mutate(ef = 'e1') ) %>%
  mutate(FDR = -log10(FDR)) %>%
  pivot_wider(id_cols = c(MODULE, NODE), names_from = ef, values_from = FDR) %>%
  mutate(d = e1-e0) %>% arrange(desc(d))


# plots -------------------------------------------------------------------

inner_join(flat, add, by =c('MODULE','NODE')) %>%
  inner_join(., biodom %>% select(GOterm_Name, Biodomain, abbr, label, color, n_hgncSymbol),
             by = c('MODULE'='GOterm_Name')) %>%
  mutate(d = e1.y-e1.x,
         d1 = e1.y-e0.y) %>%
  arrange(desc(e1.y)) %>% View()

theme_set(theme_bw())

p2=
  inner_join(flat, mult, by =c('MODULE','NODE')) %>%
  inner_join(., biodom %>% select(GOterm_Name, Biodomain, abbr, label, color, n_hgncSymbol),
             by = c('MODULE'='GOterm_Name')) %>%
  mutate(d = e1.y-e1.x,
         d1 = e1.y-e0.y) %>%
  filter(e1.y > -log10(0.1) | e1.x > -log10(0.1)) %>%
  ggplot(aes( d.y, e1.y ))+
  labs(y = 'wKDA FDR ef1, -log10', x = '\u0394 FDR, ef1 - ef0', title = dom)+
  # geom_abline(intercept = 0, slope =1, lty= 2, lwd = .5)+
  geom_smooth(method = 'glm', lty = 2, lwd = .5, color = 'grey20')+
  geom_vline(xintercept = 0, lwd = .5)+ geom_hline(yintercept = 0, lwd = .5) +
  geom_point(aes(color = color, text = paste0(NODE, '\n', MODULE, '\n', Biodomain)))+
  ggrepel::geom_label_repel(aes(label = paste0(NODE, '\n', MODULE, '\n', Biodomain)),
                            size = 2, alpha = .7, min.segment.length = 0)+
  scale_color_identity(guide = 'none') #; plotly::ggplotly(p)

plotly::ggplotly(p2)

# 
# off <- bind_rows( read_tsv('kda/off_0.results.txt') %>% mutate(ef = 'e0'), 
#                      read_tsv('kda/off_1.results.txt') %>% mutate(ef = 'e1') ) %>% 
#   mutate(FDR = -log10(FDR)) %>% 
#   pivot_wider(id_cols = c(MODULE, NODE), names_from = ef, values_from = FDR) %>% 
#   mutate(d = e1-e0) %>% arrange(desc(d))
# 
# add <- bind_rows( read_tsv('kda/add_0.results.txt') %>% mutate(ef = 'e0'), 
#                   read_tsv('kda/add_1.results.txt') %>% mutate(ef = 'e1') ) %>% 
#   mutate(FDR = -log10(FDR)) %>% 
#   pivot_wider(id_cols = c(MODULE, NODE), names_from = ef, values_from = FDR) %>% 
#   mutate(d = e1-e0) %>% arrange(desc(d))
# 
# mult <- bind_rows( read_tsv('kda/mult_0.results.txt') %>% mutate(ef = 'e0'), 
#                    read_tsv('kda/mult_1.results.txt') %>% mutate(ef = 'e1') ) %>% 
#   mutate(FDR = -log10(FDR)) %>% 
#   pivot_wider(id_cols = c(MODULE, NODE), names_from = ef, values_from = FDR) %>% 
#   mutate(d = e1-e0) %>% arrange(desc(d))
# 
# edgeRep <- read_tsv('kda/edgeRep.results.txt')
# sumOcc <- read_tsv('kda/sumOcc.results.txt')
# 
# add %>% 
#   mutate(FDR = -log10(FDR)) %>% 
#   pivot_wider(id_cols = c(MODULE, NODE), names_from = ef, values_from = FDR) %>% 
#   mutate(d = e1-e0) %>% arrange(desc(d)) %>% head()