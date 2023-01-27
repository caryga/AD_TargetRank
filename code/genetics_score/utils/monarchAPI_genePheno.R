# project directory
proj.dir <- here::here()

# packages
# suppressMessages({ 
library(httr)
library(jsonlite)
# library(furrr)
library(tidyverse)
# })

# data -- TREAT-AD model system data
synapser::synLogin()
tg <- read_csv( synapser::synGet('syn26474902')$path )

# monarch api
api_path <- 'https://api.monarchinitiative.org/api/'

# # set parallel processing
# plan(multisession, workers = 20)

# query gene phenotypes
genePheno <- list()
# genePheno <- readRDS(paste0(here::here(),'/results/monarch_genePheno.rds'))

tictoc::tic(); 
for(i in 1:nrow(tg)){
  genePheno[[i]] <- fromJSON( 
    content(
      GET(
        paste0(api_path, 
               "bioentity/gene/ENSEMBL%3A",
               tg$ENSG[i],
               "/phenotypes?rows=-1" )), 
      as='text', encoding = 'UTF-8'), flatten=T)$associations
  
  if( i %in% seq(1000, nrow(tg), by = 1000) ) {
    saveRDS(genePheno, 
            paste0(proj.dir,'/results/monarch_genePheno.rds'))
    cat(paste0(i,'...\n'))
  } 
}
tictoc::toc()

# save phenotype results
saveRDS(genePheno, 
        paste0(proj.dir,'/results/monarch_genePheno.rds'))