# Setup: pkg, synapse login -----------------------------------------------

setwd('/projects/carter-lab/caryg/hg19_liftover/')

.libPaths('/projects/carter-lab/caryg/rlib')
library(SNPlocs.Hsapiens.dbSNP151.GRCh38); library(tidyverse)

synapser::synLogin()


# GWAS to re-process ------------------------------------------------------


## DONE 
# syn = 'syn26470891' # Wightman 

## DONE 
# syn = 'syn26470907' # NG00100 AA (zipped, space delim)
#  # **1 = Model1 2 = Model2

## has rsid
# syn = 'syn26470905' # NG00075 rare (txt, space delim) -- 

## has rsid
# syn = 'syn26470900' # NG00073 LOAD_subgroups
#   #1 = lan 2 = mem 3 = mix 4 = none 5 = vsp

## DONE 
# syn = 'syn26470897' # NG00056 transethnic (zipped tbx, 4x, space delim) --
#   #1 = INT_wt_APOE **2 = in ALL sample 3 = in APOE e4 carriers 4 = in APOE e4 non-carriers

## has rsid
# syn = 'syn26470893' # NG00055 CSFbiomarkers (zipped tbx, 4x, space delim) --
#  #1 = Aβ42, 2 = pTau, 3 = Tau

## DONE
# syn = 'syn26470895' # NG00039 AfricanAmerican Reitz
#  #1 = age, sex & PCAs 2 = age, sex, PCAs, & APOE


syn = 'syn27043618' # NG00041 Beecham Neuropath
 #1 = age, sex & PCAs 2 = age, sex, PCAs, & APOE


# Pull sumstats from synapse, filter, and write bedfile -------------------

# f = synapser::synGet( syn )$path

# uf = unzip(f)
# gwas = read_delim( uf[1], delim = ' ')

# gwas = read_delim( f, delim = ' ')

# gwas = map_dfr(f, ~ read_delim(paste0('./p-value only/Results_Metaanalysis(age,sex,PCAs)_pvalueOnly/',.x), delim='\t'))

# f.gwas = gwas %>% 
#   mutate(chr = str_split_fixed(MarkerName,'-',2)[,1],
#          pos = str_split_fixed(MarkerName,'-',2)[,2]) %>%
#   rename( pval = `P-value` ) %>% 
#   # rename(chr = ` chr`, pval = p, pos = PosGRCh37) %>% 
#   filter(pval < 0.05)
# 
# f.gwas %>%
#  select(chr, st = pos, end = pos) %>%
#  mutate(chr = paste0('chr',chr)) %>%
#  write_tsv(. , paste0(syn,'_sig.bed'), col_names = F)


# # for NG00041 files
# f = synapser::synGet( syn )$path
# system(paste0('tar -xzf ',f))
f = list.files(".") %>% str_subset(., 'pvalue.')
for(i in f){
  id = i %>% str_split_fixed(., '\\.', 3) %>% as_tibble(.name_repair='universal') %>% pull(...2)
  gwas = read_delim(i, delim = ' ')
  f.gwas = gwas %>% rename(pval = Pval) %>% filter(pval < 0.05)
  f.gwas %>% 
    select(chr, st = pos, end = pos) %>% 
    mutate(chr = paste0('chr',chr)) %>% 
    write_tsv(., paste0(syn,'.',id,'.sig.bed'), col_names = F)
  }


# UCSC liftOver hg19 > GRCh38 ---------------------------------------------

# cmd = paste0('./liftOver ',syn,'_sig.bed hg19ToHg38.over.chain ', syn,'_sig_h38.bed ', syn, '_sig_unmapped.bed' )
# system(cmd)

ids <- map_chr(
  f,
  ~ .x %>% str_split_fixed(., '\\.', 3) %>% as_tibble(.name_repair='universal') %>% pull(...2)
)

for(i in ids){
  cmd = paste0('./liftOver ',syn,'.',i,'.sig.bed hg19ToHg38.over.chain ', syn,'.',i,'.sig_h38.bed ', syn, '.',i,'.sig_unmapped.bed' )
  system(cmd)
}

# system(paste0('wc ',syn,'_sig_h38.bed'))
# system(paste0('wc ',syn,'_sig_unmapped.bed'))

    ## Using rtracklayer::liftOver implementation
    #chain = import.chain('hg19ToHg38.over.chain')
    #
    #hg19_coords = f.gwas %>% 
    #  select(Chr, start = Pos, end = Pos) %>% 
    #  mutate(Chr = paste0('chr',Chr)) %>%
    #  makeGRangesFromDataFrame()
    #
    #new_coords = liftOver(hg19_coords, chain)


# Pull snp coords and map RSid --------------------------------------------


snps <- SNPlocs.Hsapiens.dbSNP151.GRCh38

for( i in ids ){
  
  newcoords = read_tsv(paste0(syn,'.',i,'.sig_h38.bed'), 
                       col_names = c('chr','start','end') )
  
  d = newcoords %>% mutate(chr = str_remove(chr, 'chr')) %>% makeGRangesFromDataFrame()
  
  genome(d) <- genome(snps)
  
  tictoc::tic()
  rs = snpsByOverlaps(snps, d) %>% 
    as.data.frame() %>% 
    mutate(chr = paste0('chr',seqnames))
  tictoc::toc()
  
  write_tsv(rs, paste0(syn,'.',i,'.mapped_rs.txt'))
  
}

# filter unmapped and join newcoords & rsid -------------------------------

for( i in ids[13:14] ){
  
  fi <- f %>% str_subset(.,i)
  f.gwas = read_delim(fi, delim = ' ') %>% rename(pval = Pval) %>% filter(pval < 0.05)
  
  newcoords = read_tsv(paste0(syn,'.',i,'.sig_h38.bed'), 
                       col_names = c('chr','start','end') )
  
  unmapped = read_tsv(paste0(syn,'.',i,'.sig_unmapped.bed'), 
                      col_names = F, comment = '#') %>% 
    mutate(Chr = str_remove(X1, 'chr'),
           pos = paste0(Chr,':',X2))
  
  rs <- read_tsv( paste0(syn,'.',i,'.mapped_rs.txt') )
  
  f.gwas.new = f.gwas %>% 
    mutate( pos2 = paste0(chr,':',pos)) %>% 
    filter( !(pos2 %in% unmapped$pos) ) %>% 
    select(-pos2) %>% 
    bind_cols(
      ., newcoords %>% select(chr_GRCh38 = chr, 
                              pos_GRCh38 = start) ) 
  
  f.gwas.new <- 
    left_join(f.gwas.new, rs, 
              by = c('chr_GRCh38' = 'chr','pos_GRCh38' = 'pos')) 
  
  # write out & upload to synapse -------------------------------------------
  
  write_tsv(f.gwas.new, paste0(syn,'.',i,'.gwas_filtered_remapped_GRCh38.tsv'))
  
  foo <- synapser::synStore( synapser::File(
    paste0(syn,'.',i,'.gwas_filtered_remapped_GRCh38.tsv'), 
    name = paste0('NG00041_GRCh38_beecham_neuropath_',i,'.tsv'),
    parent = 'syn26470873'))
  
  }

# EOF ---------------------------------------------------------------------