## HEAD: Compute genetics scores ----

# setup -------------------------------------------------------------------------------------------------------

proj.dir <- here::here()

# Load libraries
suppressMessages({  
  # library(biomaRt)
  library(tidyverse)
})
cat(' - Packages Loaded: tidyverse \n')

synapser::synLogin()

# read score component files ----------------------------------------------------------------------------------

tg.list <- read_csv( synGet('syn26474902')$path )

results <- read_csv( synapser::synGet('syn25556226')$path )
tg.del.snps <- read_csv( synapser::synGet('syn25556077')$path )
hs.phenos <-  read_csv( synapser::synGet('syn26529912')$path )
ortho.phenos <-  read_csv( synapser::synGet('syn34162624')$path )
model.data <-  read_csv( synapser::synGet('syn26546173')$path )

# 11 GB!! # tg.snp.summary <-  read_csv( synapser::synGet('syn25556197')$path )

# results modifications -----------------------------------------

results$meanRank_gwasP <- results %>% 
  mutate_at( vars(contains('minPval')), 
             ~ case_when( n_gwas == 0 ~ NA_real_,
                          !is.finite(.x) ~ NA_real_,
                          .x == 0 ~ min(.x[.x>0], na.rm=T),
                          T ~ .x) ) %>% 
  mutate_at( vars(contains('minPval')), 
             ~ rank( -log10(.x), na.last = 'keep') / length(.x[!is.na(.x)]) ) %>% 
  select(contains('minPval')) %>% rowMeans(na.rm=T)
results$meanRank_gwasP[results$n_gwas == 0] <- NA
results$meanRank_gwasP <- (results$meanRank_gwasP - min(results$meanRank_gwasP, na.rm=T)+1e-6)/(max(results$meanRank_gwasP, na.rm=T)-min(results$meanRank_gwasP, na.rm=T))

results$meanRank_qtlFDR <- results %>% 
  mutate_at( vars(contains('qtl_FDR')), 
             ~ case_when( n_qtl == 0 ~ NA_real_,
                          !is.finite(.x) ~ NA_real_,
                          .x == 0 ~ min(.x[.x>0], na.rm=T),
                          T ~ .x) ) %>% 
  mutate_at( vars(contains('qtl_FDR')), 
             ~ rank( -log10(.x), na.last = 'keep') / length(.x[!is.na(.x)]) ) %>% 
  select(contains('qtl_FDR')) %>% rowSums(na.rm=T)
results$meanRank_qtlFDR[results$n_qtl == 0] <- NA
results$meanRank_qtlFDR <- (results$meanRank_qtlFDR - min(results$meanRank_qtlFDR, na.rm=T)+1e-6)/(max(results$meanRank_qtlFDR, na.rm=T)-min(results$meanRank_qtlFDR, na.rm=T))

# calculate coding / noncoding variant summary scores ---------------------------------------------------------

# coding
tg.del.snps <- left_join(tg.del.snps, results %>% select(ENSG, biotype), by = 'ENSG')

tg.del.snps = tg.del.snps %>% 
  mutate(gl = case_when(n_codingDelSNP==0 ~ NA_real_ , n_codingDelSNP > 0 ~ -log10(gnomad_loeuf) ),
         # gl = rank(rnk_cd, ties.method = 'min', na.last = 'keep'),
         rnk_cd = case_when(n_codingDelSNP==0 ~ NA_real_ , n_codingDelSNP > 0 ~ n_codingDelSNP),
         rnk_cd = rank(rnk_cd, ties.method = 'min', na.last = 'keep'),
         rnk_fcd = case_when(n_codingDelSNP==0 ~ NA_real_ , n_codingDelSNP > 0 ~ fx_codingDel),
         rnk_fcd = rank(rnk_fcd, ties.method = 'min', na.last = 'keep') 
         ) %>% 
  relocate(gl, rnk_cd, rnk_fcd, .after=gnomad_loeuf)

tg.del.snps$gl  <- (tg.del.snps$gl - min(tg.del.snps$gl, na.rm=T))/
  (max(tg.del.snps$gl, na.rm=T)-min(tg.del.snps$gl, na.rm=T))

tg.del.snps$rnk_cd  <- (tg.del.snps$rnk_cd - min(tg.del.snps$rnk_cd, na.rm=T)+1)/
  (max(tg.del.snps$rnk_cd, na.rm=T)-min(tg.del.snps$rnk_cd, na.rm=T)+1)

tg.del.snps$rnk_fcd  <- (tg.del.snps$rnk_fcd - min(tg.del.snps$rnk_fcd, na.rm=T)+1)/
  (max(tg.del.snps$rnk_fcd, na.rm=T)-min(tg.del.snps$rnk_fcd, na.rm=T)+1)

tg.del.snps <- tg.del.snps %>% 
  mutate(
  max_spliceScore = case_when( 
    !is.finite(max_spliceScore) ~ NA_real_,
    (n_codingDelSNP == 0 & is.na(max_delRank)) ~ NA_real_,
    T ~ max_spliceScore)
)

tg.del.snps$coding_variant_summary <- tg.del.snps %>% 
  select( gl, rnk_cd, rnk_fcd, max_delRank, max_spliceScore ) %>% 
  rowSums( na.rm=T )

tg.del.snps$coding_variant_summary <- 
  (tg.del.snps$coding_variant_summary - min(tg.del.snps$coding_variant_summary, na.rm=T)) /
  (max(tg.del.snps$coding_variant_summary, na.rm=T) - min(tg.del.snps$coding_variant_summary, na.rm=T) )

tg.del.snps$coding_variant_summary[which(tg.del.snps$coding_variant_summary == 0)] <- NA
tg.del.snps$coding_variant_summary[which(tg.del.snps$biotype != 'protein_coding')] <- NA

# noncoding

tg.del.snps = tg.del.snps %>% 
  mutate(rnk_qtl = case_when(n_qtlSNP==0 ~ NA_integer_ , 
                             n_qtlSNP > 0 ~ rank( n_qtlSNP, ties.method = 'min', na.last = 'keep') ),
         ds_ds = case_when(n_qtlSNP==0 ~ NA_real_ , 
                           !is.finite(deepsea_disSig) ~ NA_real_ , 
                           (n_qtlSNP > 0 & is.finite(deepsea_disSig)) ~ deepsea_disSig),
         ds_ev = case_when(n_qtlSNP==0 ~ NA_real_ , 
                           !is.finite(deepsea_eval) ~ NA_real_ , 
                           (n_qtlSNP > 0 & is.finite(deepsea_eval)) ~ deepsea_eval),
  ) %>% 
  mutate(
    rnk_qtl = (rnk_qtl - min(rnk_qtl, na.rm=T))/(max(rnk_qtl, na.rm=T) - min(rnk_qtl, na.rm=T))/2,
    ds_ds = (ds_ds - min(ds_ds[is.finite(ds_ds)], na.rm=T))/(max(ds_ds[is.finite(ds_ds)], na.rm=T) - min(ds_ds[is.finite(ds_ds)], na.rm=T)),
    ds_ev = (ds_ev - min(ds_ev[is.finite(ds_ev)], na.rm=T))/(max(ds_ev[is.finite(ds_ev)], na.rm=T) - min(ds_ev[is.finite(ds_ev)], na.rm=T))
  ) %>% 
  relocate(rnk_qtl, ds_ds, ds_ev, .after=n_qtlSNP)

tg.del.snps$noncoding_variant_summary <- tg.del.snps %>% 
  select( rnk_qtl, mean_regulomeProb, ds_ev ) %>% 
  rowSums( na.rm=T )

tg.del.snps$noncoding_variant_summary[which(tg.del.snps$noncoding_variant_summary == 0)] <- NA
tg.del.snps$noncoding_variant_summary <- 
  (tg.del.snps$noncoding_variant_summary - min(tg.del.snps$noncoding_variant_summary, na.rm=T)+0.1) /
  (max(tg.del.snps$noncoding_variant_summary, na.rm=T) - min(tg.del.snps$noncoding_variant_summary, na.rm=T)+0.1)

# Phenotype summary scores ------------------------------------------------

model.data$Mmus_phenos[which(model.data$Mmus_ADrelPheno==0)] <- NA

model.data$Drer_pheno_score[ model.data$Drer_n_relPheno == 0 ] <- NA

model.data$model_summary = model.data %>%
  dplyr::select( Mmus_pheno_score, Drer_pheno_score, Dmel_pheno_score ) %>%
  rowSums( na.rm = T )
model.data$model_summary[model.data$model_summary==0] <- NA

model.data$model_summary <- 
  (model.data$model_summary - min(model.data$model_summary, na.rm=T) +1e-6)/
  (max(model.data$model_summary, na.rm=T) - min(model.data$model_summary, na.rm=T))

model.data$MODELAD_strain[model.data$MODELAD_strain==TRUE] <- 1
model.data$MODELAD_corr[which(is.na(model.data$MODELAD_strain))] <- NA
model.data$MODELAD <- model.data$MODELAD_strain + model.data$MODELAD_corr

model.data$model_summary[model.data$model_summary==0] <- NA

model.data <- model.data %>% 
  arrange(desc(model_summary)) %>% 
  filter(!duplicated(ENSG))

# compile score component data --------------------------------------------------------------------------------

genetics.scores <- results %>% 
  dplyr::select( ENSG , GeneName , biotype, n_gwas, min_gwasP, meanRank_gwasP, n_qtl, min_qtlFDR, meanRank_qtlFDR)
genetics.scores <- tg.del.snps %>%
  dplyr::select(ENSG , coding_variant_summary, gnomad_loeuf, n_snps, n_codingSNP, n_codingDelSNP, max_delRank, max_spliceScore, 
                noncoding_variant_summary, n_qtlSNP, min_regulomeRank, mean_regulomeProb, deepsea_eval ) %>% 
  left_join( genetics.scores, ., by="ENSG" ) %>% distinct()
genetics.scores <- hs.phenos %>%
  dplyr::select(ENSG,  Hsap_pheno_score, n_hs_adRel, hs_ad_phenotypes, n_hs_demRel, hs_dementia_phenotypes, omim_titles) %>% 
  left_join( genetics.scores, ., by="ENSG" ) %>% distinct()
genetics.scores <- model.data %>%
  dplyr::select(ENSG, MODELAD) %>% 
  left_join( genetics.scores, ., by="ENSG" ) %>% distinct()
genetics.scores <- ortho.phenos %>%
  dplyr::select(ENSG, Ortho_pheno_score, n_ortho_adRel, ortho_ad_phenotypes, n_ortho_demRel, ortho_dementia_phenotypes) %>% 
  left_join( genetics.scores, ., by="ENSG" ) %>% distinct()


# make object for score calculation and rank score components -------------------------------------------------

# genetics.scores <- genetics.scores %>% select(-score)
gs <- genetics.scores %>% 
  dplyr::select(ENSG, GeneName, n_gwas, min_gwasP, meanRank_gwasP, n_qtl, min_qtlFDR, meanRank_qtlFDR, coding_variant_summary, noncoding_variant_summary, MODELAD, Hsap_pheno_score, Ortho_pheno_score)

## Correct for NA and INF val ==
gs$n_gwas[which(gs$n_gwas==0)] <- NA
gs$min_gwasP[which(!is.finite(gs$min_gwasP))] <- NA

gs$n_qtl[which(gs$n_qtl==0)] <- NA
gs$min_qtlFDR[which(!is.finite(gs$min_qtlFDR))] <- NA

gs$MODELAD[which(gs$MODELAD == 0)] <- NA

## Calculate ranks ==
gs$n_gwas = gs$n_gwas/max(gs$n_gwas, na.rm = T)
# gs$n_gwas = rank(gs$n_gwas, ties.method = 'average', na.last = 'keep')/length(which(!is.na(gs$n_gwas)))
gs$min_gwasP = rank(-log10(gs$min_gwasP), ties.method = 'max', na.last = 'keep')/length(which(!is.na(gs$min_gwasP)))
gs$n_qtl = gs$n_qtl/max(gs$n_qtl, na.rm = T)
# gs$n_qtl = rank(gs$n_qtl, ties.method = 'average', na.last = 'keep')/length(which(!is.na(gs$n_qtl)))
gs$min_qtlFDR = rank(-log10(gs$min_qtlFDR), ties.method = 'max', na.last = 'keep')/length(which(!is.na(gs$min_qtlFDR)))

# gs$coding_variant_summary <- rank(gs$coding_variant_summary, ties.method = 'max', na.last = 'keep')/length(which(!is.na(gs$coding_variant_summary)))
# gs$noncoding_variant_summary <- rank(gs$noncoding_variant_summary, ties.method = 'max', na.last = 'keep')/length(which(!is.na(gs$noncoding_variant_summary)))

gs$MODELAD = gs$MODELAD / max(gs$MODELAD, na.rm = T)
# gs$MODELAD = rank(gs$MODELAD, ties.method = 'max', na.last = 'keep')/length(which(!is.na(gs$MODELAD)))
# gs$Hsap_pheno_score = rank(gs$Hsap_pheno_score, ties.method = 'max', na.last = 'keep')/length(which(!is.na(gs$Hsap_pheno_score)))
# gs$Mmus_pheno_score = rank( gs$Mmus_pheno_score, ties.method = 'max', na.last = 'keep') / length(which(!is.na(gs$Mmus_pheno_score)))
# gs$model_summary = rank( gs$model_summary, ties.method = 'max', na.last = 'keep') / length(which(!is.na(gs$model_summary)))


# compute genetics score --------------------------------------------------------------------------------------

gs$meanScore <- 0
gs$sumScore <- 0

gs$meanScore <- gs %>% dplyr::select(.,
                                 n_gwas, min_gwasP, meanRank_gwasP, n_qtl, min_qtlFDR, meanRank_qtlFDR,
                                 coding_variant_summary, noncoding_variant_summary,
                                 MODELAD, Hsap_pheno_score, Ortho_pheno_score) %>%
  as.matrix() %>% rowMeans(na.rm=T) * 3


gs$sumScore <- gs %>% dplyr::select(., 
                                    n_gwas, min_gwasP, meanRank_gwasP, n_qtl, min_qtlFDR, meanRank_qtlFDR,
                                    coding_variant_summary, noncoding_variant_summary,
                                    MODELAD, Hsap_pheno_score, Ortho_pheno_score) %>%
  as.matrix() %>% rowSums(na.rm=T) 

gs$meanScore <- (gs$meanScore / max(gs$meanScore)) * 3
gs$sumScore <- (gs$sumScore / max(gs$sumScore)) * 3

# incorporate score into full genetics table ------------------------------------------------------------------

genetics.scores <- gs %>% dplyr::select(ENSG, sumScore, meanScore) %>% left_join(genetics.scores, ., by='ENSG') %>% distinct()

###
# FIX GENE NAMING ISSUE
###

if( all(genetics.scores$ENSG == tg.list$ENSG) ){
  print('all ENSG match')
  genetics.scores$GeneName <- tg.list$GeneName 
  genetics.scores$biotype <- tg.list$biotype } else { print( 'mismatched ENSG' )}

###

genetics.scores <- genetics.scores %>% relocate(sumScore, meanScore, .after=GeneName) %>% arrange(desc(sumScore))

genetics.scores$score_rank <- 1+length(which(!is.na(genetics.scores$sumScore))) - rank(genetics.scores$sumScore, ties.method = 'min')
genetics.scores <- relocate(genetics.scores, score_rank, .after=sumScore)

# compile final table ----------------------------------------------

genetics.scores <- genetics.scores %>% 
  dplyr::select( ENSG, GeneName, biotype, GeneticsScore=sumScore, score_rank,
                n_gwas, min_gwasP, meanRank_gwasP, n_qtl, min_qtlFDR, meanRank_qtlFDR, n_snps, 
                coding_variant_summary, gnomad_loeuf, n_codingSNP, n_codingDelSNP, max_delRank, max_spliceScore,
                noncoding_variant_summary, n_qtlSNP, min_regulomeRank, mean_regulomeProb, deepsea_eval, 
                Hsap_pheno_score, Hsap_n_adPheno = n_hs_adRel, Hsap_ad_phenos = hs_ad_phenotypes, 
                Hsap_n_demPheno = n_hs_demRel, Hsap_dem_phenotypes = hs_dementia_phenotypes, omim_titles, 
                MODELAD, Ortholog_pheno_score = Ortho_pheno_score, ortholog_n_adPheno = n_ortho_adRel, ortholog_ad_phenotypes = ortho_ad_phenotypes, 
                ortholog_n_demPheno = n_ortho_demRel, ortholog_dem_phenotypes  = ortho_dementia_phenotypes)


# intersect with omics ----------------------------------------------------

omics <- read_csv( synapser::synTableQuery("SELECT * FROM syn22758536")$filepath )

gen.scores <- genetics.scores %>% 
  mutate(
    GeneticsScore = case_when( 
      ( ENSG %in% omics$ENSG | 
         between(GeneticsScore, quantile(GeneticsScore,0.75), 3)) ~ GeneticsScore, 
      T ~ 0
      ) # 1.15
  ) %>%
  filter( GeneticsScore != 0) %>% 
  mutate(GeneticsScore = 3*( (GeneticsScore - min(GeneticsScore, na.rm=T)) / (max(GeneticsScore, na.rm=T)-min(GeneticsScore, na.rm=T)) ) ) %>%
  arrange(desc(GeneticsScore)) 

genetics.scores <- gen.scores

# # filter & remove duplicates ----------------------------------------------------------------------------------
# 
# dup <- genetics.scores %>% select(-ENSG) %>% distinct() %>% filter(!is.na(GeneName), duplicated(GeneName)) %>% pull(GeneName)
# x <- genetics.scores %>% filter(!(GeneName %in% dup))
# x1 <- genetics.scores %>% filter(GeneName %in% dup) %>% arrange(desc(GeneticsScore)) %>% filter(!duplicated(GeneName))
# y <- rbind(x,x1) %>% arrange(desc(GeneticsScore))
# y <- y %>% filter( ENSG != 'ENSG00000270231' ) # locus in genome assembly patch
# # y %>% filter( ENSG %in% (y %>% filter(duplicated(ENSG)) %>% pull(ENSG)) | GeneName %in% (y %>% filter(duplicated(GeneName))) ) %>% View()
# genetics.scores <- y

# write genetics score file -----------------------------------------------------------------------------------

write_csv(genetics.scores, 'results/TAD_genetics_scores_Aug2022.csv')#, quote = T, row.names = F)

# generate version for synapse table display (remove NAs) -----------------------------------------------------

gen_scoreTbl <- genetics.scores
gen_scoreTbl$min_gwasP[which(is.infinite(gen_scoreTbl$min_gwasP))] <- 1
gen_scoreTbl$meanRank_gwasP[which(is.na(gen_scoreTbl$meanRank_gwasP))] <- 0
gen_scoreTbl$min_qtlFDR[which(is.infinite(gen_scoreTbl$min_qtlFDR))] <- 1
gen_scoreTbl$meanRank_qtlFDR[which(is.na(gen_scoreTbl$meanRank_qtlFDR))] <- 0
gen_scoreTbl$coding_variant_summary[which(is.na(gen_scoreTbl$coding_variant_summary))] <- 0
gen_scoreTbl$noncoding_variant_summary[which(is.na(gen_scoreTbl$noncoding_variant_summary))] <- 0
gen_scoreTbl$MODELAD[which(is.na(gen_scoreTbl$MODELAD))] <- 0
gen_scoreTbl$Hsap_pheno_score[which(is.na(gen_scoreTbl$Hsap_pheno_score))] <- 0
gen_scoreTbl$Ortholog_pheno_score[which(is.na(gen_scoreTbl$Ortholog_pheno_score))] <- 0

write_csv(gen_scoreTbl, 'results/TAD_genetics_scores_Aug2022_table.csv')


# consolidate genetics score file versions --------------------------------

foo <- synStore( 
  File('results/TAD_genetics_scores_Aug2022.csv', parent = 'syn25556053', name = 'TAD_genetics_scores.csv'),
  versionLabel = 'Aug 2022')

# SynTable

# Pull the score table, add new columns:
synId = 'syn26844312'
schema <- synGet(synId)

cols = tibble(x = as.list(synGetColumns(synId)))
cols$col = map_chr(1:nrow(cols), ~ cols$x[[.x]]$name)

# to remove
rm_idx = which(cols$col %in% setdiff(cols$col, names(gen_scoreTbl)))

for(i in rm_idx){
  schema$removeColumn(cols$x[[i]]$id)
}

# to add
gen_scoreTbl <- rename(gen_scoreTbl, 
                       # Hsap_dem_phenos = Hsap_dem_phenotypes,
                       ortholog_ad_phenos = ortholog_ad_phenotypes,
                       ortholog_dem_phenos = ortholog_dem_phenotypes
                       )
setdiff( names(gen_scoreTbl), cols$col )

newColumn1 <- synStore(Column(name = "Hsap_n_adPheno", columnType = "INTEGER"))
schema$addColumn(newColumn1)

newColumn2 <- synStore(Column(name = "Hsap_ad_phenos", columnType = "LARGETEXT"))
schema$addColumn(newColumn2)

newColumn3 <- synStore(Column(name = "Hsap_n_demPheno", columnType = "INTEGER"))
schema$addColumn(newColumn3)

newColumn4 <- synStore(Column(name = "Hsap_dem_phenos", columnType = "LARGETEXT"))
schema$addColumn(newColumn4)

newColumn5 <- synStore(Column(name = "Ortholog_pheno_score", columnType = "DOUBLE"))
schema$addColumn(newColumn5)

newColumn6 <- synStore(Column(name = "ortholog_n_adPheno", columnType = "INTEGER"))
schema$addColumn(newColumn6)

newColumn7 <- synStore(Column(name = "ortholog_ad_phenos", columnType = "LARGETEXT"))
schema$addColumn(newColumn7)

newColumn8 <- synStore(Column(name = "ortholog_n_demPheno", columnType = "INTEGER"))
schema$addColumn(newColumn8)

newColumn9 <- synStore(Column(name = "ortholog_dem_phenos", columnType = "LARGETEXT"))
schema$addColumn(newColumn9)

schema <- synStore(schema)

# Change the entire Table Entity
syn_table <- synTableQuery(sprintf("select * from %s ", synId))
deleted <- synDelete(syn_table)
synStore(Table(syn_table$tableId, gen_scoreTbl),
         used = c(
           'syn26474902',
           'syn25556226',
           'syn25556077',
           'syn26546173',
           'syn26529912',
           'syn34162624'),
         executed = 'https://github.com/caryga/treatAD_genetics/blob/6274ad2d1a868d4d82055b8ce68bf915353eecf4/scripts/computing_genetics_score.R')

# Create a snapshot of the Table Entity
body_json <- rjson::toJSON(c(snapshotComment = "New Phenotype Scores", snapshotLabel = "Aug2022"))
snapshot <-  synRestPOST(paste0("/entity/", syn_table$tableId, "/table/snapshot"), body = body_json)

## EOF ####