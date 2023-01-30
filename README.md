# AD_TargetRank
Rank and prioritize AD targets based on genetic and Genomic Evidence


## Repository Structure and Design:

All code is stored in ```code/```

The scripts to compile and compute the genetics score are stored under ```code/genetics_score/``` and the scripts to compute the multi-omics score are stored under ```code/multi-omics_score/```. Once these two scores have been computed, the overall Target Risk Score can be computed with ```code/compute_target_risk_score.R```.

##Biological domain scripts:

The ```annotate_biodomains.R``` script is used to download from BioMart the gene identifiers specific to each biological domain term. There are two scripts to calculate and process biological domain GO term enrichments: (1) ```biodom_GSEA.R``` is a wrapper around ```clusterProfiler::gseGO``` that will annotate the resulting enrichments with biological domain information, (2) ```tally_significantly_enriched_biodomains.R``` is used to generate a high-level summary of the terms enriched within each biological domain. 

##To run the genetics pipeline:

The script ```code/genetics_score/00_pipeline.sh``` is a bash wrapper around a set of R scripts. Each R script can also be run individually and this may provide more flexibility for long-running processes. The pipeline will generate a list of gene targets along with their coordinates from ```Homo_sapiens.GRCh38.104.gff3```. The target level summary script will overlap pre-filtered and cleaned genetics datasets downloaded from Synapse storage and generate summaries of minimum significance and number of significant variants for each gene and each study. There are two scripts that pull relevant phenotype data from the Monarch Initiative for human genes and also from orthologs. Variant severity metrics are calculated using Annovar scripts (for coding variants) and API calls to RegulomeDB and DeepSEA (for non-coding variants). The variant RefSNP IDs (RSIDs) need to be converted to ```.avinput``` for Annovar and ```.vcf``` files for DeepSEA prior to running these scripts. The processing of these files is detailed in the ```00_pipeline.sh``` script. 

##To run the multi-omics pipeline:

The script to run the pipeline is ```code/Main.R``` it will draw-on and run the induvidual module scripts stored in ```code/modules/```

### Pipeline input:

Input to the pipeline is in YAML format. This format will specify the target list to load and run the pipeline on, as well as which modules scripts to run. It will also specify the synapse parentID and the folder name to create a synapse folder containing the results. Induvidual config files will be stored in ```configs/```. An example YAML input file will be displayed shortly.

### Output:

The compiled output will be a Wiki on the target synapse folder. This folder will also contain the output plots in induvidual files from the script. These will be stored in ```figures/``` and that entire folder is also specified in ```.gitignore``` to prevent this repo from growing too large.

### Setting up the pipeline run environment:

All dependencies will be contained within the docker image (idealy...) and can be run in an interactive RStudio session from an AWS instance utalizing the following commands:
```
sudo yum update -y
sudo amazon-linux-extras install docker
sudo service docker start
sudo usermod -aG docker <USER ID>
sudo docker run -v "/home/jgockley/AD_TargetRank:/home/jgockley/AD_TargetRank" -e USER=$(id -u ${USER}) -e PASSWORD=<PASSWORD> -d -p 8787:8787 <DockerImageID>

#Inside the docker image make the followning file
nano ~/.synapseConfig

#[authentication]
#username = <Synapse ID>
#password = <Synapse Password>

```

### Running the pipeline:

Sourcing the Inlitiazer with your yaml configuration file copies the path into the RMarkdown file and creates the Run.Rmd file. It then sources that file and knits the markdown file into the wiki of the synapse folder specified in the configuration script and copies the figures and tables into sub folders within the destination synapse folder specified in the configuration yaml file.

### Configuration file fields:

