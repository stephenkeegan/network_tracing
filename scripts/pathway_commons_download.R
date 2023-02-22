library(tidyverse)

## FTP base url
base = 'https://www.pathwaycommons.org/archives/PC2/v13/'

## list files available via FTP
files = httr::GET(base) %>% httr::content(., as = 'text') %>% 
  str_replace_all(.,"<.*?>","") %>% 
  read_delim(., delim = ' ',skip = 7) %>% select(1:4)

## select file format to download
# format = 'sif'
format = 'txt'

## pull file names to use
sources = files %>% 
  filter(grepl(paste0(format,'.gz'), Name)) %>% 
  pull(Name)

## check for destination directory existence
if( !(format %in% list.dirs(paste0(here::here(),'/data/path_commons_data'), full.names = F)) ){
  dir.create( path = paste0('mkdir ', here::here(), '/data/path_commons_data/',format) )
}

## download specified files
for(i in 1:length(sources)){
  download.file(
    url = paste0(base,sources[i]),
    destfile = paste0(here::here(), '/data/path_commons_data/',format,'/',sources[i]),
    method = 'wget', quiet = T
      )
}

## list downloaded files
f = list.files(paste0(here::here(),'/data/path_commons_data/', format), full.names = T)

## upload to synapse
synapser::synLogin()
for(i in 1:length(f)){
  if( format == 'sif' ){ foo = synapser::synStore( synapser::File(f[i], parent = 'syn46947105') ) }
  if( format == 'txt' ){ foo = synapser::synStore( synapser::File(f[i], parent = 'syn47045122') ) }
}

