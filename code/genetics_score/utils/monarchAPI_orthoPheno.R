# project directory
proj.dir <- here::here()

# packages
# suppressMessages({ 
library(httr)
library(jsonlite)
library(furrr)
library(tidyverse)
# })

# data -- TREAT-AD model system data
synapser::synLogin()
tg <- read_csv( synapser::synGet('syn26474902')$path )

# monarch api
api_path <- 'https://api.monarchinitiative.org/api/'

# # set parallel processing
# plan(multisession, workers = 20)

# query gene ortholog phenotypes
tictoc::tic(); 
orthoPheno <- list()
for(i in 1:nrow(tg)){
  orthoPheno[[i]] <- fromJSON( 
    content(
      GET(
        paste0(api_path, 
               "bioentity/gene/ENSEMBL%3A",
               tg$ENSG[i],
               "/ortholog/phenotypes?rows=-1" )), 
      as='text', encoding = 'UTF-8'), flatten=T)$associations
  
  if( i %in% seq(1000, nrow(tg), by = 1000) ) {
    saveRDS(orthoPheno, 
            paste0(proj.dir,'/results/monarch_orthoPheno.rds'))
    cat(paste0(i,'...\n'))
  } 
}
tictoc::toc()

# save phenotype results
saveRDS(orthoPheno, 
        paste0(proj.dir,'/results/monarch_orthoPheno.rds'))