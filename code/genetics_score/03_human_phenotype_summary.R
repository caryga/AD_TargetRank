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
genePheno <- readRDS(paste0(here::here(),'/results/monarch_genePheno.rds'))

# Calculate Hsap pheno score ----

# load target list
synapser::synLogin()
hp <- read_csv( synapser::synGet('syn26474902')$path ) %>%
  dplyr::select(ENSG, GeneName, biotype) 

if(nrow(hp) != length(genePheno)) { message('Warning: gene phenotype list and target list are different lengths and may not be in sync!')}

hp <- hp %>% 
  mutate(
    n_phenotype = map_if( 1:nrow(hp), ~ (length(genePheno[[.x]]) != 0),
                          ~ genePheno[[.x]]$object.id %>% unique() %>% length(),
                          .else = ~NA_real_) %>% unlist(),
    n_ad_phenotype = map_if( 1:nrow(hp), ~ (length(genePheno[[.x]]) != 0),
                             ~ intersect(genePheno[[.x]]$object.label, ad_phenos$object.label) %>% length(),
                             .else = ~NA_real_) %>% unlist(),
    ad_phenotypes = map_if( 1:nrow(hp), ~ (length(genePheno[[.x]]) != 0),
                            ~ intersect(genePheno[[.x]]$object.label, ad_phenos$object.label) %>% sort(), 
                            .else = ~NA_character_) ,
    n_dementia_phenotype = map_if( 1:nrow(hp), ~ (length(genePheno[[.x]]) != 0),
                             ~ intersect(genePheno[[.x]]$object.label, dem_phenos$object.label) %>% length(),
                             .else = ~NA_real_) %>% unlist(),
    dementia_phenotypes = map_if( 1:nrow(hp), ~ (length(genePheno[[.x]]) != 0),
                                  ~ intersect(genePheno[[.x]]$object.label, dem_phenos$object.label) %>% sort(),
                                  .else = ~NA_character_)
    ) %>% 
  arrange(desc(n_ad_phenotype)) %>% 
  mutate(fx_ad_pheno = n_ad_phenotype / n_phenotype,
         n_ad_phenotype = case_when(n_ad_phenotype == 0 ~ NA_real_, T ~ n_ad_phenotype),
         ad_rank = rank(n_ad_phenotype, ties.method = 'min', na.last = 'keep'),
         fx_dem_pheno = n_dementia_phenotype / n_phenotype,
         n_dementia_phenotype = case_when(n_dementia_phenotype == 0 ~ NA_real_, T ~ n_dementia_phenotype),
         dem_rank = rank(n_dementia_phenotype, ties.method = 'min', na.last = 'keep')
         ) 

hp$ad_rank <- hp$ad_rank / max(hp$ad_rank, na.rm = T)
hp$dem_rank <- hp$dem_rank / max(hp$dem_rank, na.rm = T)

hp$score <- hp %>% 
  select(ad_rank, fx_ad_pheno, dem_rank, fx_dem_pheno) %>% 
  rowSums(na.rm = T)
hp$score <- hp$score / max(hp$score, na.rm = T)
hp$score[which(is.na(hp$n_phenotype))] <- NA_real_

# Retrieve OMIM titles related to TREAT-AD targets ----

# Read in OMIM titles; Jun2022
url='https://data.omim.org/downloads/ZcP5CzooRQSPLDxtmYysUQ/mimTitles.txt'
download.file(url, paste0(here::here(), '/data/mimTitles.txt'), method='libcurl', quiet = TRUE, mode = "w", cacheOK = TRUE)
# system('wget -q -O data/mimTitles.txt https://data.omim.org/downloads/ZcP5CzooRQSPLDxtmYysUQ/mimTitles.txt')
mimTitle <- read_tsv(paste0(here::here(),'/data/mimTitles.txt'), skip=2)

# Associate MIM titles to OpenAD tg
hp$omim_titles <- NA
for(i in 1:nrow(hp)){
  entrez <- org.Hs.egENSEMBL2EG[[ as.character(hp$ENSG[i]) ]]
  x <- c()
  for(e in entrez){
    x <- c(x, org.Hs.egOMIM[[e]])
  }
  hp$omim_titles[i] <- mimTitle %>% filter(`MIM Number` %in% x & `# Prefix` == 'Number Sign') %>%
    dplyr::select(`Preferred Title; symbol`) %>% as.matrix() %>% paste0(collapse=' | ')
}

## Reorder table and output summary ----

# reorder table and rename columns
hp <- hp %>% 
  dplyr::select( GeneName, ENSG, 
                 Hsap_pheno_score = score,
                 n_hs_phenotypes=n_phenotype, 
                 n_hs_adRel= n_ad_phenotype, 
                 n_hs_demRel = n_dementia_phenotype,
                 hs_ad_phenotypes= ad_phenotypes,
                 hs_dementia_phenotypes = dementia_phenotypes,
                 omim_titles) %>% 
  distinct() %>% arrange( desc(Hsap_pheno_score) ) %>% 
  rowwise() %>% 
  mutate(
    hs_ad_phenotypes = unlist(hs_ad_phenotypes) %>% paste0(., collapse = ' | '),
    hs_dementia_phenotypes = unlist(hs_dementia_phenotypes) %>% paste0(., collapse = '|')
  )

# write output
write_csv( hp, paste0(here::here(),'/results/allENSG_tg_humanPhenotype_summary.csv') )

# upload to synapse
synapser::synLogin()
foo <- synapser::synStore(
  synapser::File( paste0(here::here(),'/results/allENSG_tg_humanPhenotype_summary.csv'),
                  parent = 'syn26409209')
)

# EOF ####