#!/bin/bash
#$ -S /bin/bash
#$ -m eas -M galen.guo@canada.ca
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -N biomobpan

source ~/miniconda3/bin/activate
conda activate anvio7


GALEN=/isilon/ottawa-rdc/users/shared/chenw_lab/galen
WORK=$GALEN/BioMob/acinetobacter
RAW=$WORK/raw_data
PAN=$WORK/pangenome

anvi-gen-genomes-storage -e $PAN/acinetobacter_genome_short.txt -o $PAN/acinetobacter-GENOMES_short.db

anvi-pan-genome -g $PAN/acinetobacter-GENOMES_short.db -n biomob_acinetobacter_short -T 24 --align-with famsa -o $PAN/biomob_acinetobacter_short --enforce-hierarchical-clustering
