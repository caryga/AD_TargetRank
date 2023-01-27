#!/bin/bash

# set necessary paths (project, annovar, singularity imgs)
project_folder=</proj/folder/path>
annovar_path=</annovar/directory/path>
sif=</path/to/container> 

cd $project_folder

# instantiate log file
log="logs/$set-log_$(date +"%Y.%m.%d").txt" && touch $log
printf "\nTHIS IS THE LOG ()__) \t $(date +"%a %Y.%m.%d %T") \n" >> $log
printf "\n - Processing on: $(hostname)" >> $log

set = 'allENSG'

# Target list & summary --------------------------------------------------------

# Generate TARGETS file
printf "\n\n#### GENERATING THE TARGET SET FILE ####\n" >> $log
date +"%T" >> $log
singularity exec $sif/treatAD.sif \
  Rscript scripts/01_target_list_generation.R $log >> $log
printf "\n####\n" >> $log
echo $? >> $log

# run the target-level summary script
printf "\n\n#### RUNNING target level summary ####\n" >> $log
date +"%T" >> $log
singularity exec $sif/treatAD.sif \
  Rscript scripts/02_target_level_summary.R $log >> $log
printf "\n####\n" >> $log
echo $? >> $log

# Phenotypes -------------------------------------------------------------------

# run the human phenotype script
printf "\n\n#### SUMMARIZING human phenotype DATA ####\n" >> $log
date +"%T" >> $log
singularity exec $sif/treatAD.sif \
 Rscript scripts/03_human_phenotype_summary.R $log >> $log
printf "\n####\n" >> $log
echo $? >> $log

# run the animal model data script
printf "\n\n#### SUMMARIZING animal model DATA ####\n" >> $log
date +"%T" >> $log
singularity exec $sif/treatAD.sif \
 Rscript scripts/04_ortholog_phenotype_summary.R $log >> $log
printf "\n####\n" >> $log
echo $? >> $log

# Run RegulomeDB screen --------------------------------------------------------
  cd $project_folder
printf "\n\n#### DOWNLOADING AND COLATING RegulomeDB SCORES ... ####\n" >> $log
date +"%T" >> $log
singularity exec $sif/treatAD.sif \
 Rscript scripts/05_regulomedb_results.R $log >> $log
rm ${set}_tmp_regulomedb.tsv
printf "\n####\n" >> $log
echo $? >> $log

# Use ANNOVAR for snp ID conversion --------------------------------------------

# convert RSID to avinput for annovar input 
cd ${annovar_path}
./annotate_variation.pl -downdb -webfrom annovar -buildver hg38 avsnp150 humandb/

${annovar_path}/convert2annovar.pl -format rsid \
  ${project_folder}/results/new_rsid.txt \
  -avsnpfile ${annovar_path}/humandb/hg38_avsnp150.txt \
  > ${project_folder}/results/new_rsid_snplist.avinput
 
printf " - Using hg38 coordinates \n" 
lines=$(wc -l ${project_folder}/results/new_rsid.txt | cut -f1 -d't')
printf " - Processed $lines variants \n" 
printf " - Produced file under results/new_rsid_snplist.avinput \n" 
printf "\n####\n" 
echo $? 

# Convert RSID outputs for QTL snps to VCF for DeepSEA input (hg19 only)

wget --no-check-certificate http://www.openbioinformatics.org/annovar/download/hg19_avsnp150.txt.gz
wget --no-check-certificate http://www.openbioinformatics.org/annovar/download/hg19_avsnp150.txt.idx.gz

gunzip hg19_avsnp150.txt.gz
gunzip hg19_avsnp150.txt.idx.gz

cd ${annovar_path}
printf "\n\n#### Convert QTL RSIDs to VCF for DeepSEA input ####\n" 
${annovar_path}/convert2annovar.pl -format rsid \
  ${project_folder}/results/allENSG_qtl_rsid.txt \
  -avsnpfile ${annovar_path}/humandb/hg19_avsnp150.txt \
  | awk '{print $1"\t"$2"\t"$6"\t"$4"\t"$5}' - > ${project_folder}/results/allENSG_qtl_snplist.vcf
 
printf " - Using hg19 coordinates \n"
lines=$(wc -l ${project_folder}/results/allENSG_qtl_rsid.txt | cut -f1 -d't')
printf " - Processed $lines variants \n" 
printf " - Produced file under results/allENSG_qtl_snplist.vcf \n" 

gzip ${annovar_path}/humandb/hg19_avsnp150.txt
 
# Split VCF into files that are 10k lines long for processing
printf " - Splitting up for DeepSEA API input \n" 
cd $project_folder/results
split -a 3 -l 10000 allENSG_qtl_snplist.vcf allENSG-
vcffiles=(allENSG-*)
rm -fr deepsea && mkdir deepsea
for f in ${vcffiles[@]}; do mv ${f} deepsea/${f}.vcf ; done
cd $project_folder
 
lines=$(ls -1 results/deepsea/allENSG-* | wc -l | cut -f1 -d'r')  
printf " - Produced $lines files under results/deepsea/allENSG- \n" 
printf "\n####\n"
echo $? 

# Run DeepSEA scoring ----------------------------------------------------------

printf "\n\n#### Running DeepSEA API... ####\n" >> $log
date +"%T" >> $log
singularity exec $sif/treatAD.sif \
 Rscript scripts/deepsea_api_results.R $log >> $log
printf "\n####\n" >> $log
echo $? >> $log

# Run ANNOVAR scoring ----------------------------------------------------------

cd $annovar_path
./annotate_variation.pl -downdb -webfrom annovar -buildver hg38 ensGene humandb/
./annotate_variation.pl -downdb -webfrom annovar -buildver hg38 dbnsfp42a humandb/
./annotate_variation.pl -downdb -webfrom annovar -buildver hg38 dbscsnv11 humandb/

printf "\n\n#### Running ANNOVAR variant scoring... ####\n"

perl $annovar_path/table_annovar.pl \
 --protocol ensGene,dbnsfp42a,dbscsnv11 \
 --operation g,f,f \
 --buildver hg38 \
 --dot2underline --otherinfo --thread 10 \
 --outfile $project_folder/results/annovar/new_rsid \
  $project_folder/results/new_rsid_snplist.avinput \
  $annovar_path/humandb/ 
printf "\n####\n" 
echo $? 

printf " - Processed variants with ANNOVAR \n" 
printf "\n####\n" 
echo $? 

# cleanup
rm -fr ${annovar_path}/humandb/hg38_dbnsfp42a.txt ${annovar_path}/humandb/hg38_dbscsnv11.txt

# run the variant-level summary script -----------------------------------------
printf "\n\n#### SUMMARIZING VARIANT FUNCTIONAL CLASSIFICATION ####\n" >> $log
date +"%T" >> $log
singularity exec $sif/treatAD.sif \
 Rscript scripts/variant_level_summary.R $log >> $log
printf "\n####\n" >> $log
echo $? >> $log

# Calculate Genetics Score -----------------------------------------------------

# run the score computation script
printf "\n\n#### COMPUTING SCORES ####\n" >> $log
date +"%T" >> $log
singularity exec $sif/treatAD.sif \
 Rscript scripts/08_computing_genetics_score.R $log >> $log
printf "\n####\n" >> $log
echo $? >> $log

# wrap up
printf "\nHERE'S THE OTHER END OF THE LOG (__() \t $(date +"%T")" >> $log

## EOF ####