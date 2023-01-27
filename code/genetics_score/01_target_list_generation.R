## HEAD: Generate TEP set to score ----

## Setup & Libraries ----

# Clear environment
rm(list = ls())

# # Working directory
# setwd('/Users/caryg/projects/OpenAD')

# Load libraries
suppressMessages({ 
  library(biomaRt) 
  # library(SNPlocs.Hsapiens.dbSNP151.GRCh38) 
  library(tidyverse) 
})
cat(' - Packages Loaded: tidyverse, biomaRt \n')

# # Read arguments
# args = commandArgs(trailingOnly=TRUE)
# set = args[2]
# cat(paste0(" - Processing set: ",set," \n"))

# # Handle warning messages
# warning_file = file(args[3], open = "at") 
# sink(warning_file, append=T, type = "message")

## TEP - Target Set Import ----

system('wget -q -O Homo_sapiens.GRCh38.104.gff3.gz http://ftp.ensembl.org/pub/release-104/gff3/homo_sapiens/Homo_sapiens.GRCh38.104.gff3.gz ')
system('gunzip Homo_sapiens.GRCh38.104.gff3.gz')

chr.sizes <- read_tsv('Homo_sapiens.GRCh38.104.gff3', skip = 1, n_max = 194, col_names = F) %>% 
  mutate(X1 = str_split(X1, ' ')) %>% 
  unnest_wider(X1, names_sep = '_') %>% 
  select(chr = X1_4, start = X1_5, end = X1_6)

gff3 <- read_tsv('Homo_sapiens.GRCh38.104.gff3',
                            skip = 200,
                            col_names = F) %>% 
  filter(X3 %in% c('gene','pseudogene','ncRNA_gene')) %>% 
  mutate(gene_id = str_extract(X9, 'gene_id=.+(?=;logic)') %>% str_remove(.,'gene_id=|;'),
         gene_name = str_extract(X9, 'Name=.+(?=;biotype)') %>% str_remove(.,'Name=|;'),
         biotype = str_extract(X9, 'biotype=.+(?=;description)') %>% str_remove(.,'biotype=|;')) 

tg.list <- gff3 %>% 
  select(ENSG = gene_id, 
         GeneName = gene_name,
         biotype, 
         chromosome_name = X1,
         start_position = X4,
         end_position = X5, 
         strand = X7)

# Add target gene locus range and convert into GRanges object
tg.list$window_start <-
  as.numeric(as.character(tg.list$start_position)) - 100000
tg.list$window_end <-
  as.numeric(as.character(tg.list$end_position)) + 100000

tg.list <- tg.list %>% 
  left_join(., chr.sizes, by = c('chromosome_name'='chr')) %>% 
  mutate(window_start = case_when(window_start < 1 ~ 1, T ~ window_start)
         , window_end = case_when(window_end > as.numeric(end) ~ as.numeric(end), T ~ window_end)
         ) %>% 
  select(-start, -end)

write.csv(tg.list, paste0(here::here(),'/tep_sets/allENSG_GRCh38-104_tg_list.csv'), row.names = F)
cat(paste0(" - Processed target list output to:  tep_sets/allENSG_GRCh38-104_tg_list.csv \n"))

synapser::synLogin()
foo <- synapser::synStore( synapser::File(paste0(here::here(), '/tep_sets/allENSG_GRCh38-104_tg_list.csv'), , parent = 'syn26409209' ))

system('rm Homo_sapiens.GRCh38.104.gff3')

#EOF