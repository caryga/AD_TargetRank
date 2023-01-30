
# # point to rlibs
# host <- system('hostname', intern = T)
# if(substr(host,1,6)=='sumner'){ .libPaths('/projects/carter-lab/caryg/rlib') }

# set local directory
proj.dir <- here::here()

# setup
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(tidyverse)

# current defs:
synID <- 'syn38840892'

# Read biodom - GO ID mapping XLSX from Jesse; sheet 1 = 'complete set'
synapser::synLogin()
bd.sheets <- synapser::synGet(synID)$path %>%  readxl::excel_sheets(.)

# ## not this one!
# bd.complete <-  synapser::synGet('syn26560254')$path %>%
#   readxl::read_xlsx( ., sheet = which(bd.sheets=='Complete Set'),
#                      col_names = c('Biodomain', 'GO_Root','GO_ID', 'GOterm_Name') )

# > biodom$Biodomain %>% unique()
# [1] "Immune Response"          "Endolysosome"             "Autophagy"                "Structural Stabilization"
# [5] "Epigenetic"               "Oxidative Stress"         "APP Metabolism"           "Tau Homeostasis"         
# [9] "Apoptosis"                "Vasculature"              "Myelination"              "Lipid Metabolism"        
# [13] "RNA Spliceosome"          "Mitochondrial Metabolism" "Synapse"                  "Proteostasis"    

# > bd.sheets
# [1] "Immune response"               "Endolysosomal Process"         "Autophagy"                     "Structural Stabilization"     
# [5] "Epigenetic"                    "Oxidative stress"              "APP Metabolism"                "Tau Homeostasis"              
# [9] "Apoptosis"                     "Vascular Function"             "Myelination"                   "Lipid Metabolism"             
# [13] "RNA Spliceosome"               "Mitochondria Metabolism"       "Synaptic Function"             "Cell Cycle"                   
# [17] "Metal Binding and Homeostasis" "DNA Repair"                    "Proteastasis"                 

bd.individ <- map_dfr(
    1:length(bd.sheets),
    ~ synapser::synGet(synID)$path %>%
      readxl::read_xlsx( ., sheet = .x) %>%
      select( GO_ID = contains('ID'), GOterm_Name = contains('Term') ) %>%
      mutate(Biodomain = bd.sheets[.x])
    ) %>%
  filter( !is.na(GO_ID) ) %>%
  distinct() %>%
  mutate(
    Biodomain = case_when(
      Biodomain == 'Immune response' ~ 'Immune Response',
      Biodomain == 'Endolysosomal Process' ~ 'Endolysosome',
      Biodomain == 'Oxidative stress' ~ 'Oxidative Stress',
      Biodomain == 'Vascular Function' ~ 'Vasculature',
      Biodomain == 'Mitochondria Metabolism' ~ 'Mitochondrial Metabolism',
      Biodomain == 'Synaptic Function' ~ 'Synapse',
      Biodomain == 'Proteastasis' ~ 'Proteostasis',
      T ~ Biodomain
  ))

# Read Biodom GO ID from annotated RDS
synapser::synLogin()
biodom.current <- readRDS( synapser::synGet('syn25428992')$path )
# biodom <- biodom %>% select(1:3)

# Differences between current biodom defs and updated set
bd.individ %>% filter(!(GO_ID %in% biodom.current$GO_ID)) %>% 
  group_by(Biodomain) %>% summarise(new_term = length(unique(GO_ID))) %>% 
  arrange(desc(new_term))

# biomart query
ensembl_gene <- useMart("ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl")
# ensembl_gene <- useEnsembl("ENSEMBL_MART_ENSEMBL", 
#                            dataset = "mmusculus_gene_ensembl",
#                            mirror = 'useast')

# goid <- biodom %>% pull(GO_ID) %>% unique()
goid <- bd.individ %>% pull(GO_ID) %>% unique()
bm_data <- c()

tictoc::tic()

idx = 1102
for(i in seq(idx, length(goid), by = 367)) {
  bm_query <-
    getBM(
      attributes <- c(
        "hgnc_symbol",
        # "mgi_symbol",
        "ensembl_gene_id",
        "entrezgene_id",
        "go_id"
      ),
      filters = "go",
      values = goid[i:(i+366)],
      mart = ensembl_gene
      , useCache = T
    ) 
  bm_data <- rbind(bm_data, bm_query)
}

tictoc::toc()

# Join and summarize biomart results and biodomain table
biodom <- bm_data %>% distinct() %>% group_by(go_id) %>%
  summarise( n_ensGene = n_distinct(ensembl_gene_id),
             ensembl_id = unique(ensembl_gene_id) %>% list(), 
             n_entrezGene = n_distinct(entrezgene_id),
             entrez_id = unique(entrezgene_id) %>% list(),
             n_hgncSymbol = n_distinct(hgnc_symbol),
             hgnc_symbol = unique(hgnc_symbol) %>% list()
             # n_mgiSymbol = n_distinct(mgi_symbol),
             # mgi_symbol = unique(mgi_symbol) %>% list()
  ) %>% 
  left_join(bd.individ, ., by = c('GO_ID' = 'go_id'))
# biodom$ensembl_id <- str_split(biodom$ensembl_id, ' ')
# biodom$entrez_id <- str_split(biodom$entrez_id, ' ')
# # biodom$hgnc_symbol <- str_split(biodom$hgnc_symbol, ' ')
# biodom$mgi_symbol <- str_split(biodom$mgi_symbol, ' ')
# biodom$n_mgiSymbol <- map_dbl(1:nrow(biodom), ~ biodom$mgi_symbol[[.x]] %>% unique() %>% length() )

# Save annotated biodomains
saveRDS(biodom, file = 'annotated_biodomains_hsap.RDS')
# saveRDS(biodom, file = 'mmus_annotated_biodomains.RDS')

synapser::synLogin()
foo <- synapser::synStore( synapser::File(
  'annotated_biodomains_hsap.RDS',
  parent = 'syn23002938'
))

# # Write gene-biodomain relationships
# biodom %>% 
#   select(Biodomain, ensembl_id) %>% 
#   unnest_longer(ensembl_id) %>% 
#   distinct() %>% 
#   filter(ensembl_id != '') %>% 
#   group_by(ensembl_id) %>% 
#   summarise(n_biodoms = length(unique(Biodomain)),
#             biodoms = paste0(sort(unique(Biodomain)), collapse = ' | ')) %>% 
#   arrange(desc(n_biodoms)) %>%
#   write_csv('mmus_ensg_biodom.csv')
# 
# foo <- synapser::synStore( synapser::File(
#   'mmus_ensg_biodom.csv',
#   parent = 'syn23002938'
# ))
# 
# biodom %>% 
#   select(Biodomain, entrez_id) %>% 
#   unnest_longer(entrez_id) %>% 
#   distinct() %>% 
#   filter(entrez_id != '') %>% 
#   group_by(entrez_id) %>% 
#   summarise(n_biodoms = length(unique(Biodomain)),
#             biodoms = paste0(sort(unique(Biodomain)), collapse = ' | ')) %>% 
#   arrange(desc(n_biodoms)) %>%
#   write_csv('mmus_entrez_biodom.csv')
# 
# foo <- synapser::synStore( synapser::File(
#   'mmus_entrez_biodom.csv',
#   parent = 'syn23002938'
# ))
# 
# # biodom %>% 
# #   select(Biodomain, hgnc_symbol) %>% 
# #   unnest_longer(hgnc_symbol) %>% 
# #   distinct() %>% 
# #   filter(hgnc_symbol != '') %>% 
# #   group_by(hgnc_symbol) %>% 
# #   summarise(n_biodoms = length(unique(Biodomain)),
# #             biodoms = paste0(sort(unique(Biodomain)), collapse = ' | ')) %>% 
# #   arrange(desc(n_biodoms)) %>%
# #   write_csv('hgnc_biodom.csv')
# # 
# # foo <- synapser::synStore( synapser::File(
# #   'hgnc_biodom.csv',
# #   parent = 'syn23002938'
# # ))
# 
# biodom %>% 
#   select(Biodomain, mgi_symbol) %>% 
#   unnest_longer(mgi_symbol) %>% 
#   distinct() %>% 
#   filter(mgi_symbol != '') %>% 
#   group_by(mgi_symbol) %>% 
#   summarise(n_biodoms = length(unique(Biodomain)),
#             biodoms = paste0(sort(unique(Biodomain)), collapse = ' | ')) %>% 
#   arrange(desc(n_biodoms)) %>%
#   write_csv('mgi_biodom.csv')
# 
# foo <- synapser::synStore( synapser::File(
#   'mgi_biodom.csv',
#   parent = 'syn23002938'
# ))


