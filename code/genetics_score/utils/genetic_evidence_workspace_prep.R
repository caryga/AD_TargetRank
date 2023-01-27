### HEAD: Prepare genetic evidence datasets ----

## Setup & Libraries ----

# Clear environment
rm(list = ls())

# Working directory
proj.dir <- here::here()
data.dir <- paste0(proj.dir,'/data/source_data/filtered_GRCh38/')
# setwd('/Users/caryg/projects/OpenAD')

# Load libraries
suppressMessages({ 
  # library(SNPlocs.Hsapiens.dbSNP151.GRCh38) 
  # library(biomaRt) 
  library(tidyverse) 
})

## Load pre-filtered GWAS datasets and snp coords ----

cat(paste0(" - Start reading GWAS datasets ---"))

gwas <- list.files(data.dir) %>% str_subset('qtl', negate=T)

for(f in gwas){
  x <- assign(paste0(f %>% str_remove(., '.tsv'), '.filter'), read_tsv(paste0(data.dir,f)))
  assign(paste0(f %>% str_remove(., '.tsv'), '.coords'),
         x %>% 
           filter(!is.na(seqnames)) %>%
           select(chrom = seqnames, start = pos, end = pos) %>% 
           GenomicRanges::makeGRangesFromDataFrame()
         )
}

## Import AD e/pQTL datasets ====

# load pre-filtered eQTL data (FDR < 0.01) and filter for genes in Hit list

cat(paste0(" - Start reading e/pQTL datasets --- "))

qtl <- list.files(data.dir) %>% str_subset('qtl', negate=F)

for(f in qtl){
  assign(f %>% str_remove(., '.tsv'), read_tsv(paste0(data.dir,f)))
  }

## save image for reuse in scoring ----

rm(x, y, f, data.dir, proj.dir, gwas,qtl)

save(list = ls(),
     file='data/source_data/genetic_evidence_workspace_GRCh38.Rdata')

synapser::synLogin()
foo <- synapser::synStore(
  synapser::File('data/source_data/genetic_evidence_workspace_GRCh38.Rdata',
                 parent = 'syn26470873'))

# EOF ####