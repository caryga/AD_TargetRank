## HEAD: Overlap Gene Coordinates with GWAS/QTL SNP coordinates ----

# Setup & Libraries -------------------------------------------------------

# # Working directory
proj.dir <- here::here()

# Load libraries
suppressMessages({ 
  library(synapser)
  library(GenomicRanges)
  library(tidyverse) 
})
cat(' - Packages Loaded: tidyverse, GenomicRanges, synapser \n')

# # Read arguments
# args = commandArgs(trailingOnly=TRUE)
# set = args[1]
# cat(paste0(" - Processing set: ",set," \n"))

# # Handle warning messages
# # warning_file = file(args[2], open = "at") 
# warning_file = file(
#   paste0('logs/',Sys.time() %>% str_remove(., ' .*'), '_log.txt')
#   , open = "at") 
# sink(warning_file, append=T, type = "message")

synLogin()

# > load( synapser::synGet('syn26474865')$path )
# > save(list = ls(),
#        +      file='data/source_data/genetic_evidence_workspace_GRCh38.Rdata')
# > foo <- synapser::synStore(
#   +   synapser::File('data/source_data/genetic_evidence_workspace_GRCh38.Rdata',
#                      +                  parent = 'syn26470873'))

# Synapse download ----------------------------------------------------------

load( file= synGet('syn26474865')$path )
tg.list <- read_csv( synGet('syn26474902')$path )
results <- read_csv( synGet( 'syn25556226' )$path ) # allENSG_tg_genetic_evidence.csv 
rsid <- read_csv( synGet( 'syn26475185' )$path, col_names = 'id' ) 
rsid_qtl <- read_csv( synGet( 'syn26475190' )$path, col_names = 'id' ) 

# Load target list and genetic data ---------------------------------------

tg.coords <- tg.list %>%
  select(chrom = chromosome_name, start = window_start, end = window_end) %>%
  makeGRangesFromDataFrame()

gwas <- ls() %>% str_subset('.filter') %>% str_remove('.filter')
qtl <- ls() %>% str_subset('qtl')

# IDs of new sets to process
gwas.toProcess <- setdiff( gwas, names(results) %>% str_remove_all('_minPval'))

# Overlap target regions with GWAS SNP coordinates ------------------------

for(i in gwas.toProcess){
  assign(paste0('tg.',i,'.overlaps'), 
         suppressWarnings( 
           findOverlaps(tg.coords, 
                        eval(parse(text = paste0(i,'.coords'))))) %>% 
           as.data.frame()
         )
}

cat(paste0("--- GWAS summary stats read and overlapped \n"))

# Summarize variant-level data --------------------------------------------

cat(paste0(" - Start variant-level summarization --- "))

# TODO: generalize the variant-level summarization
# new.data <- map_dfr(
#   gwas.toProcess,
#   ~ {
#     x <- tibble(
#       study = .x,
#       data = tibble(
#         ENSG = tg.list %>% 
#           slice(eval(parse( text = paste0('tg.',.x,'.overlaps$queryHits')))) %>%
#           pull(ENSG) ,
#         tg = tg.list %>%
#           slice(eval(parse( text = paste0('tg.',.x,'.overlaps$queryHits')))) %>%
#           pull(GeneName) ,
#         snp = eval(parse( text = paste0(.x,'.filter'))) %>%
#           filter(!is.na(seqnames)) %>%
#           slice(eval(parse(text = paste0('tg.',.x,'.overlaps$subjectHits')))) %>%
#           pull(RefSNP_id)))
#     
#     eval(parse( text = paste0(.x,'.filter'))) %>%
#       rename_with(., ~ str_replace_all( .x, '.*val.*', 'pval'   ))
#       select(snp, pval, beta)
#     }
# )

bellenguez <-
  tibble(
    ENSG = tg.list %>% slice(tg.bellenguez.overlaps$queryHits) %>% pull(ENSG),
    GeneName = tg.list %>% slice(tg.bellenguez.overlaps$queryHits) %>% pull(GeneName),
    study = 'bellenguez',
    snp = bellenguez.filter %>% filter(!is.na(seqnames)) %>% 
      slice(tg.bellenguez.overlaps$subjectHits) %>% pull(RefSNP_id) 
      )
bellenguez <- bellenguez.filter %>% 
  select(snp = RefSNP_id, 
         pval = p_value,
         beta = beta) %>%
  bind_cols(.,   phenotype = rep(NA, dim(bellenguez.filter)[1])) %>% 
  left_join(bellenguez, . , by = 'snp')


cat(paste0(" --- variant-level summarization complete \n"))


# join to target-snp list -------------------------------------------------

# read exisiting list
tg.snp <- read_csv( synGet( 'syn25556224' )$path )

# merge new data tg.snp and output
tg.snp <-
  rbind(bellenguez, 
        tg.snp
        )

tg.snp <- distinct(tg.snp)

write_csv(tg.snp,'results/allENSG_tg_snps.csv')

# make rsid and output
tg.snp %>% 
  filter( !(snp %in% rsid$id) ) %>% 
  select(snp) %>% 
  distinct() %>% 
  write_tsv(.,paste0('results/NEW_rsid.txt'),col_names = F)

cat(paste0(" - Target-variants processed and summary table written to: results/allENSG_tg_snp.csv \n"))
cat(paste0(" - RSID list written to: results/allENSG_rsid.txt \n"))

# Make target-snp summary table -------------------------------------------

# tg.snp.summary <- read_csv( synGet( 'syn25556197' )$path ) 

cat(" - Start target-snp summary --- ")

tictoc::tic()
tg.snp.summary <- tg.snp %>%
  mutate( beta = case_when(study %in% c('rosmap_eqtl','mayo_eqtl','emory_pqtl') ~ beta, T ~ NaN)) %>%
  group_by(ENSG, GeneName, snp) %>%
  summarise(
    n_study = length( unique(study)),
    study = paste0( unique(study), collapse = '|'),
    min_pval = min(pval),
    qtl_beta = paste0(unique(beta %>% signif(digits = 3)), collapse=' | '),
    phenotype = paste0(sort(unique(phenotype[which(!is.na(phenotype))])), collapse=' | ')
  )
tictoc::toc()
tg.snp.summary$n_qtl <- str_count(tg.snp.summary$study, 'qtl')
tg.snp.summary <- distinct(tg.snp.summary)

tmp <- tg.snp.summary$qtl_beta %>% str_remove('NA') %>% str_trim('both') %>% str_split(pattern = '\\|')
tg.snp.summary$qtl_beta <- purrr::map_chr(tmp, . %>% unlist() %>% str_trim() %>% str_c(collapse=' | ')); rm(tmp)

cat(" --- done \n")

# write output
write_csv(
  tg.snp.summary,
  paste0('results/allENSG_tg_snp_summary.csv'),
)


cat(paste0(" - Target-snp summary file written to: results/allENSG_tg_snp_summary.csv \n"))

tg.snp.summary <- read_csv('results/allENSG_tg_snp_summary.csv')

# Summarize target-level data ---------------------------------------------

cat(paste0(" - Start target-level summarization --- "))

# results <- tg.list
results <- read_csv( synGet( 'syn25556226' )$path ) # allENSG_tg_genetic_evidence.csv 
tg.snp <- read_csv( synGet( 'syn25556224' )$path )
d = tg.snp %>% filter(study == 'bellenguez')

# loop through data, collecting results -- use parallel processing with <furrr>; set max size to 1 GB
library(furrr)

n_cores <- 10
options(future.globals.maxSize = 7e8)
plan(multisession, workers = n_cores)

x <- future_map_dfr(
  1:nrow(results),
  ~ d %>% 
    filter(ENSG == results$ENSG[.x]) %>% 
    summarise(bellenguez_nSigVar = length(unique(snp)),
              bellenguez_minPval = min(pval)),
  .progress = T )

results <- bind_cols(results,x) %>% 
  relocate(starts_with('bellenguez'), .after = 'igap_minPval')

results$n_gwas <- map_dbl(1:nrow(results),
  ~ length(which(results[.x, c(
    "bellenguez_nSigVar",
    "igap_nSigVar",
    "rosmap_nSigVar",
    "adsp_nSigVar",
    "jansen_nSigVar",
    "kunkle_nSigVar",
    "gwascat_nSigVar",
    "adni_nSigVar",
    "marioni_nSigVar",
    "moreno_grau_nSigVar",
    "schwartzentruber_nSigVar",
    "wightman_nSigVar",
    "kunkle_rare_nSigVar",
    "kunkle_aa_nSigVar",
    "reitz_aa_nSigVar",
    "jun_transeth_nSigVar",
    "load_subtype_lan_nSigVar",
    "load_subtype_mem_nSigVar",
    "load_subtype_mix_nSigVar",
    "load_subtype_none_nSigVar",
    "load_subtype_vsp_nSigVar",
    "csf_abeta42_nSigVar",
    "csf_ptau181_nSigVar",
    "csf_tau_nSigVar"
  )] != 0)))

# results$avg_nSig[i] <-
#   mean(results[i, c(
#     "igap_nSigVar",
#     "rosmap_nSigVar",
#     "adsp_nSigVar",
#     "jansen_nSigVar",
#     "kunkle_nSigVar",
#     "gwascat_nSigVar",
#     "adni_nSigVar",
#     "marioni_nSigVar",
#     "moreno_grau_nSigVar",
#     "schwartzentruber_nSigVar"
#   )] %>% as.numeric(.))

results$min_gwasP <- map_dbl(
  1:nrow(results),
  ~ results %>% slice(.x) %>% select(
    "bellenguez_minPval",
    "igap_minPval",
    "rosmap_minPval",
    "adsp_minPval",
    "jansen_minPval",
    "kunkle_minPval",
    "gwascat_minPval",
    "adni_minPval",
    "marioni_minPval",
    "moreno_grau_minPval",
    "schwartzentruber_minPval",
    "wightman_minPval",
    "kunkle_rare_minPval",
    "kunkle_aa_minPval",
    "reitz_aa_minPval",
    "jun_transeth_minPval",
    "load_subtype_lan_minPval",
    "load_subtype_mem_minPval",
    "load_subtype_mix_minPval",
    "load_subtype_none_minPval",
    "load_subtype_vsp_minPval",
    "csf_abeta42_minPval",
    "csf_ptau181_minPval",
    "csf_tau_minPval") %>% 
    min())

results$n_study <- map_dbl(
  1:nrow(results),
  ~ results %>% slice(.x) %>% select(
    'n_emory_pqtlVar',
    'n_rosmap_eqtlVar',
    'n_mayo_eqtlVar',
    "bellenguez_nSigVar",
    'igap_nSigVar',
    'rosmap_nSigVar',
    'adsp_nSigVar',
    'jansen_nSigVar',
    'kunkle_nSigVar',
    'gwascat_nSigVar',
    'adni_nSigVar',
    'marioni_nSigVar',
    'moreno_grau_nSigVar',
    'schwartzentruber_nSigVar',
    "wightman_nSigVar",
    "kunkle_rare_nSigVar",
    "kunkle_aa_nSigVar",
    "reitz_aa_nSigVar",
    "jun_transeth_nSigVar",
    "load_subtype_lan_nSigVar",
    "load_subtype_mem_nSigVar",
    "load_subtype_mix_nSigVar",
    "load_subtype_none_nSigVar",
    "load_subtype_vsp_nSigVar",
    "csf_abeta42_nSigVar",
    "csf_ptau181_nSigVar",
    "csf_tau_nSigVar") %>% t() %>% 
    as_tibble(.name_repair =  
                ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE) ) %>%     
    filter(...1 != 0) %>% nrow()    )

results$min_signif <- map_dbl(
  1:nrow(results),
  ~ results %>% slice(.x) %>% select(
    "min_emory_pqtl_FDR",
    "min_rosmap_eqtl_FDR",
    "min_mayo_eqtl_FDR",
    "bellenguez_minPval",
    "igap_minPval",
    "rosmap_minPval",
    "adsp_minPval",
    "jansen_minPval",
    "kunkle_minPval",
    "gwascat_minPval",
    "adni_minPval",
    "marioni_minPval",
    "moreno_grau_minPval",
    "schwartzentruber_minPval",
    "wightman_minPval",
    "kunkle_rare_minPval",
    "kunkle_aa_minPval",
    "reitz_aa_minPval",
    "jun_transeth_minPval",
    "load_subtype_lan_minPval",
    "load_subtype_mem_minPval",
    "load_subtype_mix_minPval",
    "load_subtype_none_minPval",
    "load_subtype_vsp_minPval",
    "csf_abeta42_minPval",
    "csf_ptau181_minPval",
    "csf_tau_minPval") %>% 
    min()    )

cat(paste0(" --- target-level summarization complete \n\n"))

write_csv( results, paste0('results/allENSG_tg_genetic_evidence.csv'))

cat(paste0(" - Target genetic evidence summary table written to: results/allENSG_tg_genetic_evidence.csv \n"))


# Synapse upload ----------------------------------------------------------

foo <- synapser::synStore(
  synapser::File('results/allENSG_tg_snps.csv', parent = 'syn26409209')
)

foo <- synapser::synStore(
  synapser::File('results/allENSG_rsid.txt', parent = 'syn26409209')
)

foo <- synapser::synStore(
  synapser::File('results/allENSG_qtl_rsid.txt', parent = 'syn26409209')
)

foo <- synapser::synStore(
  synapser::File('results/allENSG_tg_snp_summary.csv', parent = 'syn26409209')
)

foo <- synapser::synStore(
  synapser::File('results/allENSG_tg_genetic_evidence.csv', parent = 'syn26409209')
)

## EOF ####