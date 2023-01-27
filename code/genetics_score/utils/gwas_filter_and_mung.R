## HEAD: Pre-filter and prep genetic evidence for TREAT-AD scoring -----

# SETUP: load packages, synapse login -------------------------------------------------------------------------

# Working directory
proj.dir <- here::here()

# Load libraries
suppressMessages({ 
  # library(synapser)
  # library(Homo.sapiens)
  # library("org.Hs.eg.db")
  # require(SNPlocs.Hsapiens.dbSNP144.GRCh37) # old version
  require(SNPlocs.Hsapiens.dbSNP151.GRCh38)
  # library(biomaRt) 
  library(tidyverse)
})

synapser::synLogin()

# # SETUP: function definitions ---------------------------------------------------------------------------------
# 
# # snp.coords: 
# #  - updates the genomic position of reported variants to correspond to GRCh38
# #  - returns a GRanges object to be used for overlap with target gene set
# snp.coords <- function(dbSNP,set) {
#   query.snps <- snpsById(dbSNP,id=as.character(set), 
#                          ifnotfound="drop") %>% as.data.frame()
#   query.snps$end <- query.snps$pos
#   coords <- query.snps %>%
#     dplyr::select(chrom=seqnames, start=pos, end=end) %>%
#     makeGRangesFromDataFrame()
#   result.list <- list(query.snps, coords)
#   return(result.list)
# }

# # SETUP: ensure source data (GRCh38) on synapse ------------------------------------
# 
# foo = synapser::synStore( 
#   synapser::File('data/source_data/GCST90027158_buildGRCh38.tsv.gz', 
#                  'syn26470873'))

# SETUP: load snps, set threshold, pull list of source data from synapse --------------------------------------

## Set significance threshold
sig_threshold = 0.05

## Load SNPs from dbSNP
# snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
snps <- SNPlocs.Hsapiens.dbSNP151.GRCh38
snpcount(snps) %>% sum()

# syanpse source data folder
tmp <- tibble(children = as.list(synapser::synGetChildren('syn26470873'))) %>% 
  hoist(., children, 'name','id', 'modifiedOn')

gwas.sources <- tmp %>% filter(name != 'qtl')

qtl.sources <- tibble(children = as.list(synapser::synGetChildren(tmp$id[tmp$name == 'qtl']))) %>% 
  hoist(., children, 'name','id', 'modifiedOn')


# Process GWAS Data -------------------------------------------------------------------------------------------

## IGAP data -- hg19
i <- which(substr(gwas.sources$name, 1,4)=='IGAP')
gwas <- read_tsv( gzfile( synapser::synGet(gwas.sources$id[i])$path ) )
gwas.filter <- filter(gwas , substr(MarkerName,1,2)=='rs' & Pvalue < sig_threshold) 
# 7,055,881 total variants
#    19,831 positionally mapped variants (no rsid) -- 0.28%
# 6,617,243 nonsignificant variants -- 93.8%
#   437,302 filtered variants -- 6.20%
nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
tictoc::tic()
snpcoords <-  snpsById(snps,id=as.character(gwas.filter$MarkerName), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by = c('RefSNP_id'='MarkerName'))
write_tsv(gwas.filter, 'data/source_data/filtered_GRCh38/igap.tsv')
rm(gwas, gwas.filter, snpcoords)

## ROSMAP data -- hg19
i <- which(substr(gwas.sources$name, 1,6)=='rosmap')
gwas <-  read_tsv( gzfile( synapser::synGet(gwas.sources$id[i])$path ) )
gwas.filter <- filter(gwas, substr(snpId,1,2)=='rs' & pvalue < sig_threshold)
gwas.filter <- dplyr::distinct( dplyr::select( gwas.filter, snpId, effectEstimate,pvalue,alleles,minorAlleleFrequency, phenotype, riskAllele ) )
# 22,146,787 variants reported
#      8,985 positionally mapped variants (no rsid) -- 0.04%
# 20,996,296 nonsignificant variants -- 94.8%
#  1,097,979 filtered variatns -- 4.96%
nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$snpId), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by = c('RefSNP_id'='snpId'))
write_tsv(gwas.filter, 'data/source_data/filtered_GRCh38/rosmap.tsv')
rm(gwas, gwas.filter, snpcoords)

## ADSP data -- hg19
i <- which(substr(gwas.sources$name, 1,4)=='adsp')
gwas <-  read_tsv( gzfile( synapser::synGet(gwas.sources$id[i])$path ) )
gwas.filter <- filter(gwas, substr(snpId,1,2)=='rs' & pvalue < sig_threshold) 
gwas.filter <- dplyr::distinct( dplyr::select( gwas.filter, snpId, effectEstimate, pvalue, alleles, minorAllele, minorAlleleFrequency ) )
# 92,768 variants reported
#    862 positionally mapped variants (no rsid) -- 0.93%
# 86,169 nonsignificant variants -- 92.9%
#    784 filtered variants -- 0.845%
nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$snpId), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='snpId'))
write_tsv(gwas.filter, 'data/source_data/filtered_GRCh38/adsp.tsv')
rm(gwas, gwas.filter, snpcoords)

## Jansen et al data -- hg19
i <- which(substr(gwas.sources$name, 1,22)=='AD_sumstats_Jansenetal')
gwas <-  read_tsv( gzfile( synapser::synGet(gwas.sources$id[i])$path ) )
gwas.filter <- filter(gwas, substr(SNP,1,2)=='rs' & P < sig_threshold) 
# 13,367,299 variants reported
#    177,460 positionally mapped variants (no rsid) -- 1.33%
# 12,604,124 nonsignificant variants -- 94.3%
#    751,361 filtered variants -- 5.62%
nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNP), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNP'))
write_tsv(gwas.filter, 'data/source_data/filtered_GRCh38/jansen.tsv')
rm(gwas, gwas.filter, snpcoords)

## Kunkle et al data -- hg19
i <- which(substr(gwas.sources$name, 1,7)=='Kunkle_')
gwas <-  read_delim( gzfile( synapser::synGet(gwas.sources$id[i])$path ), delim=' ' )
gwas.filter <- filter(gwas, substr(MarkerName,1,2)=='rs' & Pvalue < sig_threshold) 
# 11,480,632 variants reported
#    943,919 positionally mapped variants (no rsid) -- 8.2%
# 10,813,879  nonsignificant variants -- 94.2%
#.   611,306 filtered variants -- 5.32%
nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
gwas.filter <- gwas.filter %>% 
  mutate(SNP = str_split(MarkerName,';')) %>% 
  unnest_longer(., SNP)
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNP), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNP'))
write_tsv(gwas.filter, 'data/source_data/filtered_GRCh38/kunkle.tsv')
rm(gwas, gwas.filter, snpcoords)

## ADNI data -- hg19
i <- which(substr(gwas.sources$name, 1,7)=='ADNI_Li')
gwas <- read_tsv( gzfile( synapser::synGet(gwas.sources$id[i])$path ))
gwas.filter <- filter(gwas, substr(SNP,1,2)=='rs' & p < sig_threshold) 
# 721,224 variants reported -- all have rsid and p<0.05 (Yi Li re-threshold = 5e-2)
nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNP), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNP'))
gwas.filter <- rename(gwas.filter, pos = pos.x)
write_tsv(gwas.filter, 'data/source_data/filtered_GRCh38/adni.tsv')

## GWAS Catalog data -- Nov 2021, GRCh38
i <- which(substr(gwas.sources$name, 1,12)=='gwas_catalog')
gwas <- read_tsv( gzfile( synapser::synGet(gwas.sources$id[i])$path ))
traits <- unique(gwas$MAPPED_TRAIT)
ad.traits <- traits[grep('alzheimer',traits, ignore.case = T)]
names(gwas) <- names(gwas) %>% str_replace_all(., '-| ','_')
gwas.filter <- filter(gwas, substr(SNPS,1,2) == 'rs' & MAPPED_TRAIT %in% ad.traits & P_VALUE < sig_threshold)
gwas.filter <-  dplyr::distinct( dplyr::select( gwas.filter, PUBMEDID, FIRST_AUTHOR, SNPS, CONTEXT, INTERGENIC, P_VALUE, OR_or_BETA, MAPPED_TRAIT, `REPORTED_GENE(S)`, MAPPED_GENE ) )
# 318,587 variants total
#   1,605 filtered variants: AD-associated (incl "alzheimer" in $MAPPED_TRAIT), p < 0.05 & RefSnp
gwas.filter <- gwas.filter %>% 
  mutate(SNPS = str_split(SNPS,' x |; ')) %>% 
  unnest_longer(., SNPS)
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNPS), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNPS'))
write_tsv(gwas.filter, 'data/source_data/filtered_GRCh38/gwascat.tsv')

## Marioni et al data 2018 - GCST005922
i <- which(substr(gwas.sources$name, 1,10)=='GCST005922')
gwas <-  read_delim( gzfile( synapser::synGet(gwas.sources$id[i])$path ), delim=' ' )
gwas.filter <- filter(gwas, substr(SNP,1,2)=='rs' & P < sig_threshold) 
# 7,795,605 variants reported
#   466,727 filtered variants -- 5.99%
nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNP), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNP'))
write_tsv(gwas.filter, 'data/source_data/filtered_GRCh38/marioni.tsv')
rm(gwas, gwas.filter, snpcoords)

# Schwartzentruber et al 2021 - GCST90012877
i <- which(substr(gwas.sources$name, 1,21)=='33589840-GCST90012877')
gwas <-  read_delim( gzfile( synapser::synGet(gwas.sources$id[i])$path ), delim=' ' )
gwas.filter <- filter(gwas, substr(hm_rsid,1,2)=='rs' & p_value < sig_threshold) 
# 7,795,605 variants reported
#   669,852 filtered variants -- 6.29%
nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$hm_rsid), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='hm_rsid'))
write_tsv(gwas.filter, 'data/source_data/filtered_GRCh38/schwartzentruber.tsv')
rm(gwas, gwas.filter, snpcoords)

# Moreno-Grau et al 2019 - GCST009019
i <- which(substr(gwas.sources$name, 1,10)=='GCST009019')
gwas <-  read_tsv( gzfile( synapser::synGet(gwas.sources$id[i])$path ))
gwas.filter <- filter(gwas, substr(rsID,1,2)=='rs' & P < sig_threshold) 
# 9,075,785 variants reported
#   502,044 filtered variants -- 5.53%
nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$rsID), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='rsID'))
write_tsv(gwas.filter, 'data/source_data/filtered_GRCh38/moreno_grau.tsv')
rm(gwas, gwas.filter, snpcoords)

## restarted Nov 24, 2021

# Wightman et al 2021 - **ALL POSITIONALLY MAPPED**
i <- which(substr(gwas.sources$name, 1,15)=='wightman_GRCh38')
gwas <-  read_tsv( gzfile( synapser::synGet(gwas.sources$id[i])$path ))
gwas.filter <- filter(gwas, substr(RefSNP_id,1,2)=='rs' & pval < sig_threshold)
#   719,290 filtered variants
gwas.filter <- rename(gwas.filter, pos_GRCh37 = pos, pos = pos_GRCh38)
write_tsv(gwas.filter, 'data/source_data/wightman.tsv')
rm(gwas, gwas.filter)

# NG00075 Kunkle rare
i <- which(substr(gwas.sources$name, 1,22)=='NG00075_RareVar_Kunkle')
gwas <-  read_delim( ( synapser::synGet(gwas.sources$id[i])$path ), delim =' ')
gwas.filter <- filter(gwas, substr(MarkerName,1,2)=='rs' & Pvalue < sig_threshold)
# 11,480,632 variants reported
#    611,306 filtered variants -- 5.32%
gwas.filter <- gwas.filter %>% 
  mutate(SNP = str_split(MarkerName,';')) %>% 
  unnest_longer(., SNP)
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNP), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNP'))
write_tsv(gwas.filter, 'data/source_data/kunkle_rare.tsv')
rm(gwas, gwas.filter)

# NG00100 Kunkle african american
i <- which(substr(gwas.sources$name, 1,21)=='NG00100_GRCh38_kunkle')
gwas <-  read_delim( ( synapser::synGet(gwas.sources$id[i])$path ), delim ='\t')
gwas.filter <- filter(gwas, substr(RefSNP_id,1,2)=='rs' & Pvalue < sig_threshold)
#  1,048,548 filtered variants 
gwas.filter <- rename(gwas.filter, pos = start_GRCh38)
write_tsv(gwas.filter, 'data/source_data/kunkle_aa.tsv')
rm(gwas, gwas.filter)

# NG00039 Reitz african american
i <- which(substr(gwas.sources$name, 1,20)=='NG00039_GRCh38_reitz')
gwas <-  read_delim( ( synapser::synGet(gwas.sources$id[i])$path ), delim ='\t')
gwas.filter <- filter(gwas, substr(RefSNP_id,1,2)=='rs' & pval < sig_threshold)
#  1,048,548 filtered variants 
gwas.filter <- rename(gwas.filter, pos_GRCh37 = pos, pos = pos_GRCh38)
write_tsv(gwas.filter, 'data/source_data/reitz_aa.tsv')
rm(gwas, gwas.filter)

# NG00056 Jun transethnic
i <- which(substr(gwas.sources$name, 1,18)=='NG00056_GRCh38_jun')
gwas <-  read_delim( ( synapser::synGet(gwas.sources$id[i])$path ), delim ='\t')
gwas.filter <- filter(gwas, substr(RefSNP_id,1,2)=='rs' & pval < sig_threshold)
#  384,584 filtered variants 
gwas.filter <- rename(gwas.filter, pos_GRCh37 = pos, pos = pos_GRCh38)
write_tsv(gwas.filter, 'data/source_data/jun_transeth.tsv')
rm(gwas, gwas.filter)

# NG00073 LOAD subgroups - 
i <- which(substr(gwas.sources$name, 1,22)=='NG00073_LOAD_subgroups')
f <- synapser::synGet(gwas.sources$id[i])$path
uf <- unzip(f)

gwas <-  read_delim( uf[1] , delim ='\t')
gwas.filter <- filter(gwas, substr(SNP,1,2)=='rs' & Pvalue < sig_threshold)
#  6,386,521 variants
#    316,343 filtered variants 4.95%
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNP), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNP'))
write_tsv(gwas.filter, 'data/source_data/load_subtype_lan.tsv')
rm(gwas, gwas.filter)

gwas <-  read_delim( uf[1] , delim ='\t')
gwas.filter <- filter(gwas, substr(SNP,1,2)=='rs' & Pvalue < sig_threshold)
#  6,386,521 variants
#    316,343 filtered variants 4.95%
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNP), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNP'))
write_tsv(gwas.filter, 'data/source_data/load_subtype_lan.tsv')
rm(gwas, gwas.filter)

gwas <-  read_delim( uf[2] , delim ='\t')
gwas.filter <- filter(gwas, substr(SNP,1,2)=='rs' & Pvalue < sig_threshold)
nrow(gwas); nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
#  6,396,339
#    321,770 - 5.03%
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNP), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNP'))
write_tsv(gwas.filter, 'data/source_data/load_subtype_mem.tsv')
rm(gwas, gwas.filter)

gwas <-  read_delim( uf[3] , delim ='\t')
gwas.filter <- filter(gwas, substr(SNP,1,2)=='rs' & Pvalue < sig_threshold)
nrow(gwas); nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
#  6,266,218
#    298,781 4.77%
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNP), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNP'))
write_tsv(gwas.filter, 'data/source_data/load_subtype_mix.tsv')
rm(gwas, gwas.filter)

gwas <-  read_delim( uf[4] , delim ='\t')
gwas.filter <- filter(gwas, substr(SNP,1,2)=='rs' & Pvalue < sig_threshold)
nrow(gwas); nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
#  6,398,526
#    324,058  5.06%
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNP), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNP'))
write_tsv(gwas.filter, 'data/source_data/load_subtype_none.tsv')
rm(gwas, gwas.filter)

gwas <-  read_delim( uf[5] , delim ='\t')
gwas.filter <- filter(gwas, substr(SNP,1,2)=='rs' & Pvalue < sig_threshold)
nrow(gwas); nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
#  6,389,871
#    310,243  4.86%
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNP), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNP'))
write_tsv(gwas.filter, 'data/source_data/load_subtype_vsp.tsv')
rm(gwas, gwas.filter)

# NG00055 CSF biomarkers - abeta42, tau, ptau181
i <- which(substr(gwas.sources$name, 1,7)=='NG00055')
f <- synapser::synGet(gwas.sources$id[i])$path
uf <- unzip(f)

gwas <-  read_delim( uf[1] , delim ='\t')
gwas.filter <- filter(gwas, substr(SNP,1,2)=='rs' & P < sig_threshold)
nrow(gwas); nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
#  5,991,781
#.   327,188 5.46%
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNP), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNP'))
write_tsv(gwas.filter, 'data/source_data/csf_abeta42.tsv')
rm(gwas, gwas.filter)

gwas <-  read_delim( uf[2] , delim ='\t')
gwas.filter <- filter(gwas, substr(SNP,1,2)=='rs' & P < sig_threshold)
nrow(gwas); nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
#  6,409,839
#.   331,687 5.17%
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNP), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNP'))
write_tsv(gwas.filter, 'data/source_data/csf_ptau181.tsv')
rm(gwas, gwas.filter)

gwas <-  read_delim( uf[3] , delim ='\t')
gwas.filter <- filter(gwas, substr(SNP,1,2)=='rs' & P < sig_threshold)
nrow(gwas); nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
#  5,991,780
#.   318,429 5.31%
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$SNP), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc()
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='SNP'))
write_tsv(gwas.filter, 'data/source_data/csf_tau.tsv')
rm(gwas, gwas.filter)

## ADD Bellenguez Jun 2022

# Bellenguez et al 2022 - GCST90027158
i <- gwas.sources$name %>% str_which(., 'GCST90027158')
gwas <-  read_tsv( gzfile( synapser::synGet(gwas.sources$id[i])$path ) )
gwas.filter <- filter(gwas, substr(variant_id,1,2)=='rs' & p_value < sig_threshold) 
nrow(gwas); nrow(gwas.filter); nrow(gwas.filter)/nrow(gwas)
# 21,101,114 variants reported
#  1,276,233 filtered variants -- 6.05%
tictoc::tic()
snpcoords <- snpsById(snps, id = as.character(gwas.filter$variant_id), ifnotfound="drop") %>% 
  as.data.frame() %>% as_tibble()
tictoc::toc() # 5.46 min / 327.947 sec
gwas.filter <- full_join(snpcoords, gwas.filter, by =  c('RefSNP_id'='variant_id'))
write_tsv(gwas.filter, 'data/source_data/filtered_GRCh38/bellenguez.tsv')
rm(gwas, gwas.filter, snpcoords)

# TODO: # NG00041 path


# Process QTL data --------------------------------------------------------------------------------------------

# load pre-filtered eQTL data (FDR < 0.01) and filter for genes in Hit list

cat(paste0(" - Start reading e/pQTL datasets --- "))

## ROSMAP DLPFC eQTL data -- hg19
i <- which(substr(qtl.sources$name, 1,12)=='DLPFC_ROSMAP')
qtl <-  read_csv( gzfile( synapser::synGet(qtl.sources$id[i])$path ),
                  col_names = c("chromosome", "snpLocation", "snpid", "gene", 
                                "geneSymbol", "statistic", "pvalue", "FDR", 
                                "beta", "A1", "A2", "A2freq", 
                                "expressionIncreasingAllele",
                                "strand", "geneBiotype", 
                                "geneStartPosition", "geneEndPosition" ))
write_tsv(qtl, 'data/source_data/filtered_GRCh38/rosmap_eqtl.tsv')

## MAYO TCX eQTL data -- hg19
i <- which(substr(qtl.sources$name, 1,8)=='TCX_Mayo')
qtl <- read_csv( gzfile( synapser::synGet(qtl.sources$id[i])$path ),
                 col_names = c("chromosome", "snpLocation", "snpid", "gene",
                               "geneSymbol", "statistic", "pvalue", "FDR", 
                               "beta", "A1", "A2", "A2freq",
                               "expressionIncreasingAllele",
                               "strand", "geneBiotype",
                               "geneStartPosition", "geneEndPosition")) 
write_tsv(qtl, 'data/source_data/filtered_GRCh38/mayo_eqtl.tsv')


## Emory pQTL data -- hg19
i <- which(substr(qtl.sources$name, 1,17)=='pQTLresults_fdr05')
qtl <- read_csv( gzfile( synapser::synGet(qtl.sources$id[i])$path ),
                 col_names = c("Uniprot.ID", "Chromosome", 
                               "SNP.Genomic.Position", "SNP.Genomic.Coordinate",
                               "SNP.Absolute.Genomic.Position", "Estimate", 
                               "t.value", "Standard.Error", "p.value", 
                               "log10.p.value", "Bonferroni.p.value", "FDR",
                               "Sample.Size", "Samples.with.Minor.Alleles"))

# pull ENSEMBL IDs from Uniprot IDs provided
ensembl_gene <-
  useMart("ENSEMBL_MART_ENSEMBL", 
          host = "http://www.ensembl.org",
          # host="grch37.ensembl.org",
          # path="/biomart/martservice", 
          dataset = "hsapiens_gene_ensembl")

pqtl.ens <-
  getBM(
    attributes <- c("ensembl_gene_id", "hgnc_symbol", "uniprot_gn_id"),
    filters = 'uniprot_gn_id',
    values = unique(qtl$Uniprot.ID),
    mart = ensembl_gene,
    useCache = FALSE
  )

names(pqtl.ens)[which(names(pqtl.ens) == 'uniprot_gn_id')] <- 'Uniprot.ID'
qtl <- suppressWarnings( left_join(qtl, pqtl.ens, by = 'Uniprot.ID') ) 
qtl$Chromosome <- str_remove(qtl$Chromosome, 'chr')
qtl$SNP.END <- qtl$SNP.Genomic.Position
qtl$loc <- paste0(qtl$Chromosome, ':', qtl$SNP.Genomic.Position)

pqtl.coords <- qtl %>%
  dplyr::select(chrom = Chromosome, start = SNP.Genomic.Position, end = SNP.END) %>%
  makeGRangesFromDataFrame()

suppressWarnings( library(SNPlocs.Hsapiens.dbSNP144.GRCh37) ) # old version
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
pqtl.rsid <- snpsByOverlaps(snps, pqtl.coords) %>% as.data.frame()
### hg19!!!! ###
pqtl.rsid$loc <- paste0(pqtl.rsid$seqnames, ":", pqtl.rsid$pos)

qtl <- left_join(qtl, pqtl.rsid, by = 'loc')

write_tsv(qtl, 'data/source_data/filtered_GRCh38/emory_pqtl.tsv')



# prep workspace ----------------------------------------------------------

rm(list = ls())

data.dir <- paste0(here::here(),'/data/source_data/filtered_GRCh38/')

gwas <- list.files(data.dir) %>% str_subset('qtl', negate=T)

for(f in gwas[1]){
  x <- assign(paste0(f %>% str_remove(., '.tsv'), '.filter'), read_tsv(paste0(data.dir,f)))
  assign(paste0(f %>% str_remove(., '.tsv'), '.coords'),
         x %>% 
           filter(!is.na(seqnames)) %>%
           select(chrom = seqnames, start = pos, end = pos) %>% 
           GenomicRanges::makeGRangesFromDataFrame()
  )
}

rm(x, f, data.dir,gwas )

# load current workspace file
synapser::synLogin()
load( synapser::synGet('syn26474865')$path )

save(list = ls(),
     file='data/source_data/genetic_evidence_workspace_GRCh38.Rdata')

foo <- synapser::synStore(
  synapser::File('data/source_data/genetic_evidence_workspace_GRCh38.Rdata',
                 parent = 'syn26470873'))


## EOF ####
