#!/bin/bash

# activate synapse conda env and download synapse data
conda activate synapse
synapse get syn16984409 --downloadLocation data # DLPFC_ROSMAP_cis_eQTL_release.csv
synapse get syn16984410 --downloadLocation data # TCX_Mayo_cis_eQTL_release
synapse get syn21213340 --downloadLocation data # pQTLresults.csv

# Filter QTL data for FDR < 0.05
awk -F "\"*,\"*" '$8  <= 0.05 {print $0}' data/DLPFC_ROSMAP_cis_eQTL_release.csv > data/DLPFC_ROSMAP_cis_eQTL_release_fdr05.csv
awk -F "\"*,\"*" '$8  <= 0.05 {print $0}' data/TCX_Mayo_cis_eQTL_release.csv > data/TCX_Mayo_cis_eQTL_release_fdr05.csv
awk -F "\"*,\"*" '$12 <= 0.05 {print $0}' data/pQTLresults.csv > data/pQTLresults_fdr05.csv

#EOF
