## HEAD: Integrate VEP data into variant-level information ----

# Working directory
proj.dir <- here::here()

# Load libraries
suppressMessages({
  library(furrr)
  library(tidyverse) 
  })
cat(' - Packages Loaded: tidyverse \n')

plan('multisession')


# synapser::synLogin()
# tg.snp.summary <- read_csv( synapser::synGet('syn25556197')$path )
tg.snp.summary <- read_csv( paste0(here::here(),'/results/allENSG_tg_snp_summary.csv'))

##  ANNOVAR results ----

x = list.files(path=paste0(here::here(),'/results/annovar'), pattern='hg38_multianno.txt', full.name=T)
annovar.results <- read_tsv(x[1])
annovar.results <- bind_rows( annovar.results,
                              read_tsv(x[2], col_types = spec(annovar.results)))

# annovar.file <- list.files(path=paste0(proj.dir,'/results/annovar'), pattern='allENSG.hg38_multianno.txt', full.name=T)
# annovar.results <- read_tsv(annovar.file)

names(annovar.results) <- names(annovar.results) %>% str_replace_all('-','_')

# Add specific variant information to OpenAD list
tg.snp.summary <- annovar.results %>%
  dplyr::select(. , snp = Otherinfo1,
                Gene_ensGene, Func_ensGene, ExonicFunc_ensGene, AAChange_ensGene, Interpro_domain,
              GTEx_V8_gene, GTEx_V8_tissue,
              SIFT_converted_rankscore, SIFT_pred,
                SIFT4G_converted_rankscore, SIFT4G_pred,
              Polyphen2_HDIV_rankscore, Polyphen2_HDIV_pred,
                Polyphen2_HVAR_rankscore, Polyphen2_HVAR_pred,
              LRT_converted_rankscore, LRT_pred,
              MutationTaster_converted_rankscore, MutationTaster_pred,
              MutationAssessor_rankscore, MutationAssessor_pred,
              FATHMM_converted_rankscore, FATHMM_pred,
              PROVEAN_converted_rankscore, PROVEAN_pred,
                VEST4_rankscore,
              MetaSVM_rankscore, MetaSVM_pred,
              MetaLR_rankscore, MetaLR_pred,
                MetaRNN_rankscore, MetaRNN_pred,
              M_CAP_rankscore, M_CAP_pred,
                REVEL_rankscore,
                MutPred_rankscore,
                MVP_rankscore,
                MPC_rankscore,
                PrimateAI_rankscore, PrimateAI_pred,
                DEOGEN2_rankscore, DEOGEN2_pred,
                BayesDel_addAF_rankscore, BayesDel_addAF_pred,
                BayesDel_noAF_rankscore, BayesDel_noAF_pred,
                ClinPred_rankscore, ClinPred_pred,
                LIST_S2_rankscore, LIST_S2_pred,
                Aloft_pred, Aloft_Confidence,
              CADD_raw_rankscore, CADD_phred,
              DANN_rankscore, DANN_score,
                fathmm_MKL_coding_rankscore, fathmm_MKL_coding_pred,
                fathmm_XF_coding_rankscore, fathmm_XF_coding_pred,
                Eigen_raw_coding, Eigen_raw_coding_rankscore,
                Eigen_PC_raw_coding, Eigen_PC_raw_coding_rankscore,
                GenoCanyon_score, GenoCanyon_rankscore,
                integrated_fitCons_score, integrated_fitCons_rankscore, integrated_confidence_value,
                LINSIGHT, LINSIGHT_rankscore,
              dbscSNV_ADA_SCORE, dbscSNV_RF_SCORE) %>%
  left_join(tg.snp.summary,.,by='snp') %>% distinct()

saveRDS(tg.snp.summary, 'results/tg_snp_summary.rds')
# write_csv(  tg.snp.summary, paste0(proj.dir,'/results/allENSG_tg_snp_summary.csv'))
# cat(paste0(" - IN PROGRESS Target-snp summary file written to: results/allENSG_tg_snp_summary.csv"))

##
# Coding & Splcing variants ----
##

tg.snp.summary <- readRDS(paste0(here::here(), '/results/tg_snp_summary.rds'))

# tg.snp.summary <- read_csv( 
#   paste0(proj.dir,'/results/allENSG_tg_snp_summary.csv')
#   # ,col_types = 'cccdcddcdcccccccdcdcdcdcccdcdcdcdcddcdcdcccdccddcdcdcdcdcdcccdddddcdcdddddddddcdddd' 
#   )

tictoc::tic()
tg.snp.summary$tg.coding <- 0
tg.snp.summary$tg.coding[which(
  sapply( 1:nrow(tg.snp.summary),
          function(i) grepl(tg.snp.summary$GeneName[i], as.character(tg.snp.summary$Gene_ensGene[i]))) &
    tg.snp.summary$Func_ensGene %in% c('exonic', 'exonic;splicing', 'splicing') )] <- 1
tictoc::toc() # 7.198 min

# Count up predictions of deleterious coding variants
x <- dplyr::select(tg.snp.summary, SIFT_pred, Polyphen2_HDIV_pred, LRT_pred,
                   MutationTaster_pred, MutationAssessor_pred, FATHMM_pred, PROVEAN_pred,
                   MetaSVM_pred, MetaLR_pred, M_CAP_pred
                   , SIFT4G_pred, Polyphen2_HVAR_pred,  MetaRNN_pred,
                   PrimateAI_pred,DEOGEN2_pred,
                   BayesDel_addAF_pred,BayesDel_noAF_pred,
                   ClinPred_pred, LIST_S2_pred, Aloft_pred,
                   fathmm_MKL_coding_pred, fathmm_XF_coding_pred
                   )
x <- x %>% mutate( across(everything(), ~ as.character(.x) ))
x[is.na(x)] <- '.'
y <- str_glue_data(x, "{SIFT_pred}{Polyphen2_HDIV_pred}{LRT_pred}{MutationTaster_pred}{MutationAssessor_pred}",
              "{FATHMM_pred}{PROVEAN_pred}{MetaSVM_pred}{MetaLR_pred}{M_CAP_pred}{PROVEAN_pred}"
              ,"{SIFT4G_pred}{Polyphen2_HVAR_pred}{MetaRNN_pred}{PrimateAI_pred}{DEOGEN2_pred}"
              , "{BayesDel_addAF_pred}{BayesDel_noAF_pred}{ClinPred_pred}{LIST_S2_pred}" # {Aloft_pred}
              , "{fathmm_MKL_coding_pred}{fathmm_XF_coding_pred}"
              ) %>% str_count("D|A|H|M|TRUE")

tg.snp.summary$deleterious <- y 
rm(x)

# convert numeric scores to numeric vectors; were being handled as factor levels
tg.snp.summary$CADD_phred <- tg.snp.summary$CADD_phred %>% as.matrix() %>% as.numeric()
tg.snp.summary$DANN_score <- tg.snp.summary$DANN_score %>% as.matrix() %>% as.numeric()
tg.snp.summary$Eigen.raw <- tg.snp.summary$Eigen_raw_coding %>% as.matrix() %>% as.numeric()

# add CADD, DANN and Eigen.raw score cutoffs to deleterious
tg.snp.summary$deleterious[which(tg.snp.summary$CADD_phred > 10)] = 1+ tg.snp.summary$deleterious[which(tg.snp.summary$CADD_phred > 10)]
tg.snp.summary$deleterious[which(tg.snp.summary$DANN_score > 0.95)] = 1+ tg.snp.summary$deleterious[which(tg.snp.summary$DANN_score > 0.95)]
tg.snp.summary$deleterious[which(tg.snp.summary$Eigen_raw_coding > 0.5)] = 1+ tg.snp.summary$deleterious[which(tg.snp.summary$Eigen_raw_coding > 0.5)]

# Score variants based on the average rankscores for deleteriousness
x <- tg.snp.summary %>%
  dplyr::select(
  SIFT_converted_rankscore ,
  SIFT4G_converted_rankscore,
  Polyphen2_HDIV_rankscore ,
  Polyphen2_HVAR_rankscore,
    # LRT_converted_rankscore ,
    # MutationTaster_converted_rankscore ,
  MutationAssessor_rankscore ,
  FATHMM_converted_rankscore ,
  PROVEAN_converted_rankscore ,
  VEST4_rankscore,
  MetaSVM_rankscore ,
  MetaLR_rankscore ,
  MetaRNN_rankscore,
  M_CAP_rankscore ,
  CADD_raw_rankscore ,
  DANN_rankscore,
  REVEL_rankscore,
  MutPred_rankscore,
  MVP_rankscore,
  MPC_rankscore,
  PrimateAI_rankscore,
  DEOGEN2_rankscore,
  BayesDel_addAF_rankscore,
  BayesDel_noAF_rankscore,
  ClinPred_rankscore,
  LIST_S2_rankscore,
  fathmm_MKL_coding_rankscore,
  fathmm_XF_coding_rankscore,
  Eigen_raw_coding_rankscore,
  Eigen_PC_raw_coding_rankscore,
  GenoCanyon_rankscore,
  integrated_fitCons_rankscore,
  LINSIGHT_rankscore
  ) 

x <- x %>% mutate(across(everything(), ~ as.numeric(.x) ))
y <- x %>% rowMeans(., na.rm = T )
tg.snp.summary$delRank_mean <- y
rm(x)

# Score splicing variants for tg genes:
tg.snp.summary$spliceScore_mean <- NA
x <- which(names(tg.snp.summary)=='dbscSNV_ADA_SCORE'|names(tg.snp.summary)=='dbscSNV_RF_SCORE')
tg.snp.summary$spliceScore_mean[which(
  sapply( 1:nrow(tg.snp.summary), function(i) grepl(tg.snp.summary$GeneName[i], 
                                                    as.character(tg.snp.summary$Gene_ensGene[i])))
  )] <-
  rowMeans( tg.snp.summary[which(
    sapply( 1:nrow(tg.snp.summary), function(i) grepl(tg.snp.summary$GeneName[i], 
                                                      as.character(tg.snp.summary$Gene_ensGene[i])))
  ), x], na.rm=T)
rm(x)

saveRDS(tg.snp.summary, 'results/tg_snp_summary.rds')
# write_csv(  tg.snp.summary, paste0(proj.dir,'/results/allENSG_tg_snp_summary.csv'))
cat(paste0(" - IN PROGRESS Target-snp summary file written to: results/allENSG_tg_snp_summary.csv"))

##  DeepSEA results ----
synapser::synLogin()
deepsea.dis <- read_tsv( synapser::synGet('syn26532726')$path )
deepsea.ev <- read_tsv( synapser::synGet('syn26533095')$path )

tg.snp.summary <- deepsea.dis %>% #data.frame() %>%
  dplyr::select(. , snp = sequence.name, ds_dis_score=value) %>%
  left_join(tg.snp.summary,.,by='snp')
 
tg.snp.summary <- deepsea.ev %>% #data.frame() %>%
  dplyr::select(. , snp = sequence.name, ds_eval=value) %>%
  left_join(tg.snp.summary,.,by='snp')

tg.snp.summary$ds_dis_score[which(tg.snp.summary$n_qtl == 0)] <- NA
tg.snp.summary$ds_eval[which(tg.snp.summary$n_qtl == 0)] <- NA

tg.snp.summary <- distinct(tg.snp.summary)

saveRDS(tg.snp.summary, 'results/tg_snp_summary.rds')
# write_csv(  tg.snp.summary, paste0(proj.dir,'/results/allENSG_tg_snp_summary.csv'))
cat(paste0(" - IN PROGRESS Target-snp summary file written to: results/allENSG_tg_snp_summary.csv"))

##  RegulomeDB results ----

synapser::synLogin()
regulome.db <- read_csv( synapser::synGet('syn26530109')$path )

tg.snp.summary <- regulome.db %>% 
  dplyr::select(., snp= rsids, regulomeDB_rank = ranking, regulomeDB_prob = probability) %>%
  left_join(tg.snp.summary, ., by='snp')

tg.snp.summary$regulomeDB_rank[which(tg.snp.summary$n_qtl == 0)] <- NA
tg.snp.summary$regulomeDB_prob[which(tg.snp.summary$n_qtl == 0)] <- NA

## Write Variant-level summary of variant function ----
saveRDS(tg.snp.summary, 'results/tg_snp_summary.rds')
# write_csv(  tg.snp.summary, paste0(proj.dir,'/results/allENSG_tg_snp_summary.csv'))
cat(paste0(" - FINAL Target-snp summary file written to: results/allENSG_tg_snp_summary.csv"))

foo <- synapser::synStore(
  synapser::File(paste0(here::here(),'/results/tg_snp_summary.rds'),
                 parent = 'syn26409209',
                 name = 'allENSG_tg_snp_summary.rds'))

##  Summarize variant info per tg ----

# Make table summarizing variants-by-target
tg.del.snps <- dplyr::group_by(tg.snp.summary, ENSG, GeneName) %>% 
  dplyr::select(., ENSG, GeneName) %>% data.frame() %>% distinct()

# Calculate variant properties
tg.del.snps <- tg.snp.summary %>% filter(ENSG %in% tg.del.snps$ENSG) %>% group_by(ENSG) %>% 
  summarise(.groups='keep', n_snps = n_distinct(snp)) %>%
  left_join(tg.del.snps, ., by='ENSG')
tg.del.snps <- tg.snp.summary %>% filter(ENSG %in% tg.del.snps$ENSG) %>% group_by(ENSG) %>% 
  filter(tg.coding==1) %>% summarise(.groups='keep', n_codingSNP = n_distinct(snp)) %>%
  left_join(tg.del.snps, ., by='ENSG')
tg.del.snps <- tg.snp.summary %>% filter(ENSG %in% tg.del.snps$ENSG) %>% group_by(ENSG) %>% 
  filter(tg.coding==1 & deleterious>2) %>% summarise(.groups='keep', n_codingDelSNP = n_distinct(snp)) %>%
  left_join(tg.del.snps, ., by='ENSG')
tg.del.snps$fx_codingDel <- NA
tg.del.snps <- tg.snp.summary %>% filter(ENSG %in% tg.del.snps$ENSG) %>% group_by(ENSG) %>% 
  filter(tg.coding==1) %>% summarise(.groups='keep', max_delRank = max(delRank_mean, na.rm=T)) %>%
  left_join(tg.del.snps, ., by='ENSG')
tg.del.snps <- tg.snp.summary %>% filter(ENSG %in% tg.del.snps$ENSG) %>% group_by(ENSG) %>% 
  summarise(.groups='keep', max_spliceScore = max(spliceScore_mean, na.rm=T)) %>%
  left_join(tg.del.snps, ., by='ENSG')
tg.del.snps <- tg.snp.summary %>% filter(ENSG %in% tg.del.snps$ENSG) %>% group_by(ENSG) %>% 
  filter(n_qtl > 0) %>% summarise(.groups='keep', 
                                  n_qtlSNP =  n_distinct(snp),
                                  min_regulomeRank =  first( sort(as.character(regulomeDB_rank), na.last=T) ),
                                  max_regulomeProb =  max(regulomeDB_prob, na.rm=T),
                                  max_regulomeProb =  max(regulomeDB_prob, na.rm=T),
                                  mean_regulomeProb = mean(regulomeDB_prob, na.rm=T),
                                  deepsea_disSig =  max(ds_dis_score, na.rm=T),
                                  deepsea_eval = max(ds_eval, na.rm=T),
                                  absv_qtl_beta = max(abs(as.numeric(unlist( str_split( qtl_beta, '\\|') ))), na.rm=T),
                                  qtl_beta_range = paste0(range(as.numeric(unlist( str_split(qtl_beta, '\\|') )), na.rm=T)[1],',',
                                                          range(as.numeric(unlist( str_split(qtl_beta, '\\|') )), na.rm=T)[2]),
                                  pos_beta_ratio = log2(
                                    as.numeric(unlist( str_split( qtl_beta, '\\|') )) %>% as.data.frame %>% filter(. > 0) %>% nrow() /
                                      as.numeric(unlist( str_split( qtl_beta, '\\|') )) %>% as.data.frame %>% filter(. < 0) %>% nrow()
                                  )  ) %>%  
  left_join(tg.del.snps, ., by='ENSG')

# Re-factor NA/Inf in some characteristics
tg.del.snps$n_codingSNP[which(is.na(tg.del.snps$n_codingSNP))] <- 0
tg.del.snps$n_codingDelSNP[which(is.na(tg.del.snps$n_codingDelSNP))] <- 0
tg.del.snps$fx_codingDel <- tg.del.snps$n_codingDelSNP/tg.del.snps$n_codingSNP
tg.del.snps$n_qtlSNP[which(is.na(tg.del.snps$n_qtlSNP))] <- 0
tg.del.snps$max_delRank[which(!(is.finite(tg.del.snps$max_delRank)))] <- NA
tg.del.snps$max_spliceScore[which(!(is.finite(tg.del.snps$max_spliceScore)))] <- NA

# tg.del.snps <- tg.del.snps %>% rename(ENSG=ensg)

write_csv( tg.del.snps, paste0(proj.dir,'/results/allENSG_tg_variantFunction_data.csv') )
cat(paste0(" - IN PROGRESS Target-snp summary file written to: results/allENSG_tg_variantFunction_data.csv"))

## GNOMAD pLoF scores ----

system('wget -q -O data/gnomad_lof_metrics_by_gene.txt.bgz https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')
gnomad <- read_tsv( gzfile(paste0(proj.dir,'/data/gnomad_lof_metrics_by_gene.txt.bgz'))) #, header=T)
# gnomad.filter <- dplyr::select(gnomad, gene, ENSG=gene_id, oe_lof_upper)
    # also: oe_lof_upper_rank, oe_lof_upper_bin
tg.del.snps <- dplyr::select(gnomad, ENSG=gene_id, gnomad_loeuf=oe_lof_upper) %>% 
  left_join(tg.del.snps, ., by='ENSG') %>% relocate(gnomad_loeuf, .before=n_snps)


## Write Target-level summary of variant function ----

write_csv( tg.del.snps, paste0(proj.dir,'/results/allENSG_tg_variantFunction_data.csv') )
cat(paste0(" - FINAL Target-snp summary file written to: results/allENSG_tg_variantFunction_data.csv"))

# synapser::synLogin()
foo <- synapser::synStore( synapser::File(
  paste0(here::here(),'/results/allENSG_tg_variantFunction_data.csv'),
  parent = 'syn26409209'
))

## EOF ####