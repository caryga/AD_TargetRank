# Submit jobs to the DeepSEA API interface

# Load required packages ----

# Working directory
proj.dir <- here::here()

# Load libraries
suppressMessages({ 
  library(tidyverse)
  library(httr)
  library(jsonlite)
  })

# Declare path to the API interface ----
api_path <- "https://hb.flatironinstitute.org/api/deepsea/jobs/"

  # # Test API interface; example commands can be found at: https://hb.flatironinstitute.org/api
  # # e.g. for GET /api/deepsea/
  # fromJSON(httr::content(GET(path), as='text', encoding = 'UTF-8'),flatten=T)

# ID the VCF files for scoring ----
vcf.files <- list.files(paste0(proj.dir, '/results/deepsea'))
# system(" grep -v '-' results/deepsea/",vcf.files[i]," > ",vcf.files[i],".vcf ")

# Post all vcf jobs ====
file.jobs <- tibble( file = vcf.files, jobid = NA)
# for( f in 1 ){
for(f in 1:nrow(file.jobs)){
  
  id=file.jobs$file[f]
  
  # file = paste0("@/Users/caryg/projects/treatAD_genetics/results/deepsea/", id) # MLG
  file = paste0("@",proj.dir,"/results/deepsea/", id) 
  
  title = str_remove(id,".vcf")
  
  # POST the VCF file to the DeepSEA server for scoring ----
  # must be GRCh37/hg19 coordinates 
  # each VCF should be limited to 10,000 variants per submission (break up larger files)
  message(paste0('Posting VCF: ', id))
  # cmd <- paste0("curl -F upload_file=\"",file,"\" -F input_type=vcf -F ds_model=beluga -F title=", title, " ", api_path)
  cmd <- paste0("curl -F upload_file=\"",file,"\" -F input_type=vcf -F ds_model=sei -F genome_assembly=hg38 -F title=", title, " ", api_path)
  post <- system( cmd, intern = T )
  
  file.jobs$jobid[f] <- fromJSON(post, flatten = T)$job_id
  
  write_tsv(file.jobs, paste0(proj.dir, '/results/',
                              Sys.Date() %>% str_replace_all(.,'-','_'),'deepsea_file_jobs.txt'))
}
print(file.jobs, n = nrow(file.jobs))
write_tsv(file.jobs, paste0(proj.dir, '/results/',
                            Sys.Date() %>% str_replace_all(.,'-','_'),'deepsea_file_jobs.txt'))

# Wait while jobs run ====

file.jobs$status <- map_chr(
  1:nrow(file.jobs),
  ~ fromJSON( httr::content( GET(paste0(api_path,file.jobs$jobid[.x])) , as='text', encoding = 'UTF-8'),  flatten=T)$status
)

t=900
# message(paste0('running, please wait... ', t/60,'min'))

Sys.sleep(900)

repeat{  
  t=t+900
  file.jobs$status <- map_chr(
    1:nrow(file.jobs),
    ~ fromJSON( httr::content( GET(paste0(api_path,file.jobs$jobid[.x])) , as='text', encoding = 'UTF-8'),  flatten=T)$status
  )
  if(any(file.jobs$status!='completed')){
    # message(paste0('running, please wait... ', t,'s'))
    Sys.sleep(900)
  } else{ 
    message(paste0('all deepsea jobs completed after ', t/60,' minutes!'))
    break}}  

# Retrieve scores ====

seqclass <- map_dfr(
  1,#:nrow(file.jobs),
  ~fromJSON( httr::content( GET(paste0(api_path,file.jobs$jobid[1],'/variant_scores/?score_type=seqclass' )), as='text', encoding = 'UTF-8'), flatten=T)[['scores']] 
)
  
diffs <- map_dfr(
  1:nrow(file.jobs),
  ~
)

# dis_scores <- map_dfr(
#   1:nrow(file.jobs),
#   ~ fromJSON( httr::content( GET(paste0(api_path,file.jobs$jobid[.x],'/variant_scores/?score_type=dis')), as='text', encoding = 'UTF-8'), flatten=T)[['scores']]
# )
# 
# mneval_scores <- map_dfr(
#   1:nrow(file.jobs),
#   ~ fromJSON( httr::content( GET(paste0(api_path,file.jobs$jobid[.x],'/variant_scores/?score_type=mean-evalue')), as='text', encoding = 'UTF-8'), flatten=T)[['scores']]
# )

      # dis_scores <- c()
      # mneval_scores <- c()
      # for(i in 1:nrow(file.jobs)){
      #   # Retrieve scores ====
      #   d.scores <- fromJSON( httr::content( GET(paste0(api_path,file.jobs$jobid[i],'/variant_scores/?score_type=dis')), as='text', encoding = 'UTF-8'), flatten=T)[['scores']]
      #   dis_scores <- rbind(dis_scores, d.scores)
      #   
      #   e.scores <- fromJSON( httr::content( GET(paste0(api_path,file.jobs$jobid[i],'/variant_scores/?score_type=mean-evalue')), as='text', encoding = 'UTF-8'), flatten=T)[['scores']]
      #   mneval_scores <- rbind(mneval_scores, e.scores)
      # }

# Write output ====
synapser::synLogin()
# 
# write_tsv(dis_scores, paste0(proj.dir,'/results/deepsea/allENSG_deepsea_dis.txt'))
# foo <- synapser::synStore(
#   synapser::File( paste0(proj.dir,'/results/deepsea/allENSG_deepsea_dis.txt'),
#                   parent = 'syn26409209')
# )

write_tsv(mneval_scores, paste0(proj.dir,'/results/deepsea/allENSG_deepsea_mean-evalue.txt'))
foo <- synapser::synStore(
  synapser::File( paste0(proj.dir,'/results/deepsea/allENSG_deepsea_mean-evalue.txt'),
                  parent = 'syn26409209')
)


# Variant scores
# Disease Impact Score (DIS): DIS is calculated by training a logistic regression model that prioritizes likely 
# disease-associated mutations on the basis of the predicted transcriptional or post-transcriptional regulatory effects 
# of these mutations (See Zhou et. al, 2019). The predicted DIS probabilities are then converted into ‘DIS e-values’, 
# computed based on the empirical distributions of predicted effects for gnomAD variants. The final DIS score is:
#   
# −𝑙𝑜𝑔10(𝐷𝐼𝑆𝑒𝑣𝑎𝑙𝑢𝑒𝑓𝑒𝑎𝑡𝑢𝑟𝑒)
# 
# Mean -log e-value (MLE): For each predicted regulatory feature effect (𝑎𝑏𝑠(𝑝𝑎𝑙𝑡−𝑝𝑟𝑒𝑓)) of a variant, w
# e calculate a ‘feature e-value’ based on the empirical distribution of that feature’s effects among gnomAD variants 
# (see above Regulatory feature scores: e-value). The MLE score of a variant is
# 
# ∑−𝑙𝑜𝑔10(𝑒𝑣𝑎𝑙𝑢𝑒𝑓𝑒𝑎𝑡𝑢𝑟𝑒)/𝑁

## EOF ####