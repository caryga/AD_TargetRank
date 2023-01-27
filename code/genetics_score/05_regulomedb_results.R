## HEAD: Download RegulomeDB tables of scores for eQTL variants ----

## SETUP: Load libraries and data ----

# # Working directory
proj.dir <- here::here()

# Load libraries
suppressMessages({ 
  library(tidyverse)
  # library(httr)
  # library(jsonlite)
})
  cat(' - Packages Loaded: tidyverse \n')

# Read data files
qtl.snps <- read_tsv(paste0(proj.dir,'/results/allENSG_qtl_rsid.txt'), col_names = 'snp')

# Declare path to the API interface ----
api_path <- "https://www.regulomedb.org/regulome-summary/"

##  Query RegulomeDB and download results ----

# start clock
ptm <- proc.time()[3]

last = which(qtl.snps$snp == 'rs34414297') # 'rs77552996')


for(i in seq(1, to=nrow(qtl.snps), by=200) ){
# for(i in seq(last+1, to=nrow(qtl.snps), by=200) ){
url=paste0(api_path,'?regions=',
            paste0(qtl.snps$snp[i:(i+199)] %>% na.omit() ,collapse='%0D%0A'),
            '&genome=GRCh37&format=tsv', collapse='')
download.file(url, paste0(proj.dir, '/results/allENSG_tmp_regulomedb.tsv'), method='libcurl', quiet = TRUE, mode = "a", cacheOK = TRUE)
}
  # proc.time()[3] - ptm
  # # ~8 minutes to download 10k snps; ~1 MB total file size
  # # 3.5 hr to download 273k; 27.2 MB

## Open results file and clean up the header rows ----

reg <- read_tsv( paste0(proj.dir, '/results/allENSG_tmp_regulomedb.tsv') ) %>% 
  mutate(QTL = str_remove_all(QTL, 'chrom.*')) %>% 
  filter(chrom!='chrom') %>% distinct()

  cat( paste0(' - Downloaded info for ', nrow(reg), ' variants from RegulomeDB \n'))
  cat( paste0(' - Processed in ', round((proc.time()[3]-ptm)/60, digits = 2), ' minutes \n'))

  ## Write the results file for further processing ----

write_csv(reg, paste0(proj.dir,'/results/allENSG_regulomeDB_results.csv'))
  
synapser::synLogin()
foo <- synapser::synStore(
    synapser::File( paste0(proj.dir,'/results/allENSG_regulomeDB_results.csv'),
                    parent = 'syn26409209')
  )

#### EOF ####