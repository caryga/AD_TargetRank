## HEAD: Score targets based on observed/known human phenotypes ----

## SETUP: Load libraries and data ----

# Load libraries
suppressMessages({
  library('org.Hs.eg.db')
  library(httr)
  library(jsonlite)
  library(tidyverse)
})
cat(' - Packages Loaded: org.Hs.eg.db, httr, jsonlite, tidyverse \n')


# Get Dz phenotypes from Monarch/MONDO ----

api_path <- 'https://api.monarchinitiative.org/api/'

ad_phenos <- fromJSON( 
  content(
    GET(
      paste0(api_path,
             "bioentity/disease/",
          # ad MONDO:0004975
             "MONDO%3A0004975",
             "/phenotypes?rows=-1" )), 
    as='text', encoding = 'UTF-8'), flatten=T)$associations %>% 
  select(object.iri, object.id, object.label) %>% distinct()

dem_phenos <- fromJSON( 
  content(
    GET(
      paste0(api_path,
             "bioentity/disease/",
      # dementia MONDO:0001627
             "MONDO%3A0001627",
             "/phenotypes?rows=-1" )), 
    as='text', encoding = 'UTF-8'), flatten=T)$associations %>% 
  select(object.iri, object.id, object.label) %>% distinct()


## See also:
# cognitive disorder MONDO:0002039 ## super-class of dementia
# tauopathy MONDO:0005574 ## super-class of AD
# neurodegen dz MONDO:0005559 ## super-class of tauopathy

# Get gene phenotypes from Monarch/uPheno ----

# source(paste0(here::here(), '/scripts/monarchAPI_genePheno.R')) ## ~12 hours to run 
orthoPheno <- readRDS(paste0(here::here(),'/results/monarch_orthoPheno.rds'))

# Calculate Orholog pheno score ----

# load target list
synapser::synLogin()
op <- read_csv( synapser::synGet('syn26474902')$path ) %>%
  dplyr::select(ENSG, GeneName, biotype) 

if(nrow(op) != length(orthoPheno)) { message('Warning: gene phenotype list and target list are different lengths and may not be in sync!')}

best_matches <- bind_rows(
  read_tsv(paste0(here::here(),'/data/hp-to-mp-bestmatches.tsv'), 
           col_names = c('hp_id','hp_label','orth_id','orth_label','fuzzy_equivalence', 'fuzzy_subclass')),
  read_tsv(paste0(here::here(),'/data/hp-to-zp-bestmatches.tsv'), 
           col_names = c('hp_id','hp_label','orth_id','orth_label','fuzzy_equivalence', 'fuzzy_subclass')),
  read_tsv(paste0(here::here(),'/data/hp-to-wbphenotype-bestmatches.tsv'), 
           col_names = c('hp_id','hp_label','orth_id','orth_label','fuzzy_equivalence', 'fuzzy_subclass'))
  )

op$ortho_phenos <- map_if(1:length(orthoPheno), ~(length(orthoPheno[[.x]]) != 0),
                          ~ orthoPheno[[.x]] %>% select(object.id, object.label) %>% distinct() %>% 
                            left_join(., best_matches, by = c('object.id'='orth_id','object.label'='orth_label') ),
                          .else = ~NA_character_)

op <- op %>% 
  mutate(
    n_ortho_pheno = map_if(1:nrow(op), ~ (length(op$ortho_phenos[[.x]]) > 1),
                           ~ op$ortho_phenos[[.x]]$object.id %>% unique() %>% length(),
                           .else = ~NA_real_) %>% unlist(),
    n_ad_pheno = map_if(1:nrow(op), ~ (length(op$ortho_phenos[[.x]]) > 1),
                        ~ op$ortho_phenos[[.x]] %>% filter(object.id %in% ad_phenos$object.id | hp_id %in% ad_phenos$object.id) %>% 
                          pull(object.id) %>% unique() %>% length(),
                        .else = ~NA_real_) %>% unlist(),
    n_dem_pheno = map_if(1:nrow(op), ~ (length(op$ortho_phenos[[.x]]) > 1),
                         ~ op$ortho_phenos[[.x]] %>% filter(object.id %in% dem_phenos$object.id | hp_id %in% dem_phenos$object.id) %>% 
                           pull(object.id) %>% unique() %>% length(),
                         .else = ~NA_real_) %>% unlist() ,
    ad_phenotypes = map_if( 1:nrow(op), ~ (length(op$ortho_phenos[[.x]]) > 1),
                            ~ op$ortho_phenos[[.x]] %>% filter(object.id %in% ad_phenos$object.id | hp_id %in% ad_phenos$object.id) %>%  
                              pull(object.label) %>% unique() %>% sort(),
                            .else = ~NA_character_) ,
    dementia_phenotypes = map_if( 1:nrow(op), ~ (length(op$ortho_phenos[[.x]]) > 1),
                                  ~ op$ortho_phenos[[.x]] %>% filter(object.id %in% dem_phenos$object.id | hp_id %in% dem_phenos$object.id) %>%  
                                    pull(object.label) %>% unique() %>% sort(),
                                  .else = ~NA_character_)
    ) %>% 
  arrange(desc(n_ad_pheno)) %>% 
  mutate(fx_ad_pheno = n_ad_pheno / n_ortho_pheno,
         n_ad_pheno = case_when(n_ad_pheno == 0 ~ NA_real_, T ~ n_ad_pheno),
         ad_rank = rank(n_ad_pheno, ties.method = 'min', na.last = 'keep'),
         fx_dem_pheno = n_dem_pheno / n_ortho_pheno,
         n_dem_pheno = case_when(n_dem_pheno == 0 ~ NA_real_, T ~ n_dem_pheno),
         dem_rank = rank(n_dem_pheno, ties.method = 'min', na.last = 'keep')
         ) 

op$ad_rank <- op$ad_rank / max(op$ad_rank, na.rm = T)
op$dem_rank <- op$dem_rank / max(op$dem_rank, na.rm = T)

op$score <- op %>% 
  select(ad_rank, fx_ad_pheno, dem_rank, fx_dem_pheno) %>% 
  rowSums(na.rm = T)
op$score <- op$score / max(op$score, na.rm = T)
op$score[which(is.na(op$n_ortho_pheno))] <- NA_real_

# # Get Dz model info from Monarch ----
# 
# ad_model <- fromJSON( 
#   content(
#     GET(
#       paste0(api_path,
#              "bioentity/disease/",
#              "MONDO%3A0004975",
#              "/models?rows=-1" )), 
#     as='text', encoding = 'UTF-8'), flatten=T)$associations %>% 
#   select(object.iri, object.id, object.label, contains('evidence')) %>% distinct()
# 
# dem_model <- fromJSON( 
#   content(
#     GET(
#       paste0(api_path,
#              "bioentity/disease/",
#              "MONDO%3A0001627",
#              "/models?rows=-1" )), 
#     as='text', encoding = 'UTF-8'), flatten=T)$associations %>% 
#   select(object.iri, object.id, object.label, contains('evidence')) %>% distinct()
# 

## Reorder table and output summary ----

# reorder table and rename columns
op <- op %>% 
  dplyr::select( GeneName, ENSG, 
                 Ortho_pheno_score = score,
                 n_ortho_pheno, 
                 n_ortho_adRel= n_ad_pheno, 
                 n_ortho_demRel = n_dem_pheno,
                 ortho_ad_phenotypes = ad_phenotypes,
                 ortho_dementia_phenotypes = dementia_phenotypes) %>% 
  distinct() %>% arrange( desc(Ortho_pheno_score) ) %>% 
  rowwise() %>% 
  mutate(
    ortho_ad_phenotypes = unlist(ortho_ad_phenotypes) %>% paste0(., collapse = ' | '),
    ortho_dementia_phenotypes = unlist(ortho_dementia_phenotypes) %>% paste0(., collapse = '|')
  )

# write output
write_csv( op, paste0(here::here(),'/results/allENSG_tg_orthologPhenotype_summary.csv') )

# upload to synapse
synapser::synLogin()
foo <- synapser::synStore(
  synapser::File( paste0(here::here(),'/results/allENSG_tg_orthologPhenotype_summary.csv'),
                  parent = 'syn26409209')
)


# # plot --------------------------------------------------------------------
# 
# ad_gwas <- c('APP', 'PSEN1', 'PSEN2', 'APOE', 'SORL1', 
#              # 'CR1', 'IL34','IQCK','CASS4', 'HLA-DRB1','PILRA',  'PTK2B',
#              'CD33', 'CLU', 'PICALM', 'BIN1', 
#              # 'MS4A4A', 
#              'ABCA7', 'SPI1', 'CD2AP', 'SLC2A4', 'INPP5D', 'EPHA1',  'FERMT2', 'TREM2', 
#              # 'TREML2', 'ECHDC3', 'HS3ST1', 
#              'MAPT', 'ABI3', 'PLCG2', 'ACE', 'ADAM10',  'APH1B', 'SCIMP', 
#              'MINK1', 'PRKD3',  'SHARPIN',
#              
#              'TMEM106B','GRN',
#              'PRNP', 'BCKDK','WDR81','UNC5C','WWOX','ABCA1','AKAP9','PICALM')
# 
# x <- inner_join( 
#   hp %>% select(ENSG, GeneName, Hsap_pheno_score, n_hs_phenotypes),
#   op %>% select(ENSG, GeneName, Ortho_pheno_score, n_ortho_pheno),
#   by = c('ENSG','GeneName'))
# 
# ggplot(x, aes(Hsap_pheno_score, Ortho_pheno_score)) + #, text = GeneName
#   geom_point()+ geom_abline(intercept = 0, slope = 1, lty = 2)+ theme_bw() +
#   ggrepel::geom_label_repel(data = subset(x, GeneName %in% ad_gwas), 
#   aes(label = GeneName), size = 3, max.overlaps = 20)

# EOF ####