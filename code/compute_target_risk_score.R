
# load libraries ----------------------------------------------------------

require('synapser')
require( 'tidyverse')
synLogin()

# Calculate Overall scores ------------------------------------------------------------------------------------

# read in score components ------------------------------------------------------------------------------------

# remove minimum values from genetics scores -- want to set these to zero
gs <- read_csv( synTableQuery('SELECT * FROM syn26844312')$filepath )  %>%
# gs <- read_csv( 'results/TAD_genetics_scores_Jan2022.csv' )  %>%   
  mutate(GeneticsScore = case_when(GeneticsScore!=min(GeneticsScore)~GeneticsScore)) %>%
  mutate(isScored_genetics = case_when( !is.na(GeneticsScore)~'Y', T~'N' )) %>% 
  filter(isScored_genetics == 'Y') %>% 
  select(ENSG, GeneName, GeneticsScore,isScored_genetics)

# omics scores
omics <- read_csv(synTableQuery("SELECT * FROM syn22758536")$filepath) %>% 
  select(ENSG, GeneName = GName, OmicsScore) %>% 
  mutate(isScored_omics = 'Y') %>% 
  arrange(desc(OmicsScore))

# literature scores
litScores <- read_csv(synTableQuery("SELECT * FROM syn23569348")$filepath) %>% 
  # mutate(isScored_lit = case_when( !is.na(litScore)~'Y', T~'N' )) %>% 
  # filter(isScored_lit == 'Y') %>%
  mutate(isScored_lit = 'Y') %>% 
  select(ENSG, GeneName = symbol, litScore, isScored_lit)

# Shulman lab neuropath score
flynp <- readxl::read_excel(synGet('syn26348115')$path,sheet='Sheet1') %>% 
  mutate(FlyNeuroPathScore = 2*(Rating / max(Rating)) ) %>% 
  arrange(desc(FlyNeuroPathScore)) %>% 
  filter(!duplicated(`Human Gene`)) %>%  
  mutate(isScored_neuropath = 'Y') %>%
  select(GeneName = `Human Gene`, FlyNeuroPathScore, isScored_neuropath)

# Combine omic, genetics, lit, flynp scores -------------------------------------------------------------------

tad_score <- full_join( 
  dplyr::select(gs, ENSG, GeneName, GeneticsScore, isScored_genetics), 
  dplyr::select(omics, ENSG, OmicsScore, isScored_omics), 
  by=c('ENSG'='ENSG')) %>% arrange(desc(OmicsScore)) 

tad_score <- litScores %>% dplyr::select(ENSG, GeneName, LiteratureScore = litScore, isScored_lit) %>% 
  left_join(tad_score,., by=c('ENSG'='ENSG','GeneName'='GeneName'))
  # 1 litScore omitted -- ENSG00000176782 DEFB104A

tad_score <- flynp %>% 
  left_join(tad_score,., by='GeneName')
  # 2 flynp scores omitted -- AL031587.1 & GI1

# read in and attach druggability info ------------------------------------------------------------------------

druggability <- data.table::fread(synGet('syn13363443')$path,data.table=F)
# geneTable <- utilityFunctions::convertEnsemblToHgnc(druggability$GeneID)
# druggability <- dplyr::left_join(druggability,geneTable,by=c('GeneID'='ensembl_gene_id'))
tad_score <- druggability %>% left_join(tad_score, ., by=c('ENSG'='GeneID'))
tad_score <- tad_score %>% relocate(starts_with('isScored'), 
                                    .after = 'tissue_engagement_bucket_definition')

# set isScored flags and zeros -------------------------------------------------

x <- tad_score %>% 
  mutate(
    isScored_genetics = case_when(!is.na(GeneticsScore)~'Y', T~isScored_genetics),
    isScored_omics = case_when(!is.na(OmicsScore)~'Y', T~isScored_omics),
    isScored_lit = case_when(!is.na(LiteratureScore)~'Y', T~isScored_lit),
    isScored_neuropath = case_when(!is.na(FlyNeuroPathScore)~'Y', T~isScored_lit),
  ) %>% 
  mutate(
    isScored_genetics = case_when(!is.na(isScored_genetics)~isScored_genetics,T~'N'),
    isScored_omics = case_when(!is.na(isScored_omics)~isScored_omics,T~'N'),
    isScored_lit = case_when(!is.na(isScored_lit)~isScored_lit,T~'N'),
    isScored_neuropath = case_when(!is.na(isScored_neuropath)~isScored_neuropath,T~'N')
  ) %>% 
  mutate(
    GeneticsScore = case_when(is.na(GeneticsScore)~0,T~GeneticsScore),
    OmicsScore = case_when(is.na(OmicsScore)~0,T~OmicsScore),
    LiteratureScore = case_when(is.na(LiteratureScore)~0,T~LiteratureScore),
    FlyNeuroPathScore = case_when(is.na(FlyNeuroPathScore)~0,T~FlyNeuroPathScore),
  )
tad_score  <-  x

# calculate Overall Risk score -------------------------------------------------

tad_score$Overall <- tad_score %>% 
  dplyr::select(GeneticsScore, OmicsScore, LiteratureScore, FlyNeuroPathScore) %>%  
  rowSums(., na.rm=T) 

tad_score$Overall_rank <- 1+length(which(!is.na(tad_score$Overall))) - rank(tad_score$Overall, ties.method = 'min')

tad_score <- tad_score %>%  arrange(desc(Overall)) %>% relocate( Overall, Overall_rank, .after=GeneName )

# write tad score file ---------------------------------------------------------

write_csv(tad_score, paste0(here::here(), '/results/TAD_overall_scores_Feb2022.csv'))#, quote=T, row.names = F)

# consolidate overall score syn file versions --------------------------------

# SynFiles
foo <- synStore(
  File(
    paste0(here::here(),'/results/TAD_overall_scores_Feb2022.csv'), 
    parent = 'syn25556053', 
    name = 'TAD_overall_scores.csv'),
  versionLabel = 'Feb 2022 n2')

# update synapse table with new columns -----------------------------------

scores <- tad_score
names(scores) <- names(scores) %>% str_replace_all(., ' ', '_')

# Pull Version 3 of the scores, add new columns:
synId = 'syn25575156'

# Change the entire Table Entity
tad_table <- synTableQuery(sprintf("select * from %s ", synId))
deleted <- synDelete(tad_table)
synStore(Table(tad_table$tableId, scores),
         used = c('syn26844312','syn22758536','syn23569348','syn26348115','syn13363443'),
         executed = 'https://github.com/caryga/treatAD_overall/blob/9656ca69f03e9fdbd5cf236f17c8ce3c111ec76e/scripts/computing_scores_NEW2.R')

# Create a snapshot of the Table Entity
body_json <- rjson::toJSON(c(snapshotComment = "Version 9 Scores Feb 2022; updated Omics", snapshotLabel = "v9"))
snapshot <-  synRestPOST(paste0("/entity/", tad_table$tableId, "/table/snapshot"), body = body_json)


# #syns used:
# syn26844312
# syn22758536
# syn23569348
# syn26348115
# syn13363443

## EOF ####