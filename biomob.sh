#!/bin/bash
#$ -S /bin/bash



source ~/miniconda3/bin/activate
conda activate anvio7


## CHANGE SAMPLE NAME HERE
## CHANGE SAMPLE NAME HERE
## CHANGE SAMPLE NAME HERE
## CHANGE SAMPLE NAME HERE

GALEN=/isilon/ottawa-rdc/users/shared/chenw_lab/galen
WORK=$GALEN/BioMob/acinetobacter
RAW=$GALEN/BioMob/raw_data

for SAMPLE in `awk '{print $1}' $RAW/sample_name`;
do
echo $SAMPLE

#directories


TRIM=$WORK/$SAMPLE/trimming_clean
ASSEMBLY=$WORK/$SAMPLE/assembly
ANVIO=$WORK/$SAMPLE/anvio/
MAPPING=$WORK/$SAMPLE/mapping/



####################################################################
### QC of fastqc pre-trim
####################################################################
mkdir $WORK/$SAMPLE/
mkdir $WORK/$SAMPLE/fastqc
fastqc $RAW/$SAMPLE* -t 10 -o $WORK/$SAMPLE/fastqc/ -f fastq


###################################################################
### trimming with trimmomatic, removal of adapter using bbduk or fastp
###################################################################

echo Trimming

mkdir $WORK/$SAMPLE/trimming_clean
TRIM=$WORK/$SAMPLE/trimming_clean

#fastp

fastp \
--in1 $RAW/${SAMPLE}_R1.fastq.gz \
--in2 $RAW/${SAMPLE}_R2.fastq.gz \
--out1 $TRIM/${SAMPLE}_clean1.fq.gz \
--out2 $TRIM/${SAMPLE}_clean2.fq.gz \
--html $TRIM/${SAMPLE}.html \
--json $TRIM/${SAMPLE}.json \
-f 10 \
-F 10 \
-c 


####################################################################
### co-assembly megahit
### 
####################################################################


ASSEMBLY=$WORK/$SAMPLE/assembly
megahit -1 $TRIM/${SAMPLE}_clean1.fq.gz -2 $TRIM/${SAMPLE}_clean2.fq.gz -o $ASSEMBLY -t 80 --min-contig-len 500


####################################################################
### loading into anvio ### 
####################################################################

mkdir $WORK/$SAMPLE/anvio/
ANVIO=$WORK/$SAMPLE/anvio/
anvi-script-reformat-fasta $ASSEMBLY/final.contigs.fa -o $ANVIO/$SAMPLE.fa  -l 500 --simplify-names -r $ANVIO/$SAMPLE_report.txt


###create anvio database

#Gene calling, taxonomical, functional assignment to contigs and generate summary


anvi-gen-contigs-database -f $ANVIO/$SAMPLE.fa -o $ANVIO/$SAMPLE.db -n $SAMPLE -T 80
anvi-run-hmms -c $ANVIO/$SAMPLE.db -T 80
anvi-run-kegg-kofams -c $ANVIO/$SAMPLE.db -T 80 --hmmer-program hmmsearch --just-do-it
anvi-run-ncbi-cogs -c $ANVIO/$SAMPLE.db --cog-data-dir $GALEN/COG --num-threads 80 --search-with diamond
anvi-run-pfams -c $ANVIO/$SAMPLE.db --pfam-data-dir $GALEN/database/pfam --hmmer-program hmmscan --num-threads 80
anvi-run-scg-taxonomy -c $ANVIO/$SAMPLE.db -T 80 
anvi-scan-trnas  -c $ANVIO/$SAMPLE.db --log-file $ANVIO/log_trna.txt --trna-hits-file trna_hit.txt
anvi-estimate-scg-taxonomy -c $ANVIO/$SAMPLE.db -o $ANVIO/${SAMPLE}_SCG.txt --just-do-it

anvi-display-contigs-stats $ANVIO/$SAMPLE.db  --report-as-text  --output-file $ANVIO/${SAMPLE}_contigs_summary.txt
anvi-get-sequences-for-gene-calls -c $ANVIO/$SAMPLE.db --get-aa-sequences -o $ANVIO/${SAMPLE}_aa.fa


####################################################################
### mapping with bowtie2
### 
####################################################################

echo mapping

## create mapping directory

mkdir $WORK/$SAMPLE/mapping/
MAPPING=$WORK/$SAMPLE/mapping/




bowtie2-build $ANVIO/$SAMPLE.fa  $MAPPING/$SAMPLE --threads 100

## transfering sample name file to mapping folder.

### mapping read to contigs

bowtie2 -x $MAPPING/$SAMPLE \
-1 $TRIM/${SAMPLE}_clean1.fq.gz \
-2 $TRIM/${SAMPLE}_clean2.fq.gz \
-S $MAPPING/$SAMPLE.sam \
--threads 100



echo convert sam to bam
conda deactivate
samtools view -S -b $MAPPING/${SAMPLE}.sam > $MAPPING/${SAMPLE}.bam


## bam to anvio_bam profiling
conda activate anvio7
anvi-init-bam $MAPPING/${SAMPLE}.bam -o $MAPPING/${SAMPLE}_anvi.bam

rm $MAPPING/*.bt2
rm $MAPPING/*.sam

echo mapping_done

####################################################################
### Profiling  bam --> anvio to ensure clean genomes
####################################################################

anvi-profile -i $MAPPING/${SAMPLE}_anvi.bam -c $ANVIO/$SAMPLE.db  \
--output-dir $ANVIO/profile --sample-name $SAMPLE -T 100 

#add collection name
anvi-script-add-default-collection -c $ANVIO/$SAMPLE.db -p $ANVIO/profile/PROFILE.db -C $SAMPLE 

#generate new fasta file with contigs associated only with taxa of choice.
anvi-summarize -c $ANVIO/$SAMPLE.db -p $ANVIO/profile/PROFILE.db -C $SAMPLE -o $ANVIO/$SAMPLE


done




conda activate /isilon/ottawa-rdc/users/shared/chenw_lab/galen/gtdbtk

GTDBTK_DATA_PATH="/isilon/common/reference/databases/gtdb/release202"


for SAMPLE in `awk '{print $1}' $RAW/sample_name`;
do
echo $SAMPLE

#directories


TRIM=$WORK/$SAMPLE/trimming_clean
ASSEMBLY=$WORK/$SAMPLE/assembly
ANVIO=$WORK/$SAMPLE/anvio/
MAPPING=$WORK/$SAMPLE/mapping/

gtdbtk ani_rep --genome_dir $ANVIO/$SAMPLE/bin_by_bin/EVERYTHING/ --out_dir $ANVIO/$SAMPLE/gtdb_output -x fa --cpus 20 --prefix $SAMPLE 

done


gtdbtk ani_rep --genome_dir $WORK/all_fas/ --out_dir $WORK/all_fas/gtdb_output -x fa --cpus 20 

### pangenomic?

#create external file

name	contigs_db_path
Name_01	/path/to/contigs-01.db
Name_02	/path/to/contigs-02.db
Name_03	/path/to/contigs-03.db
(…)	(…)

ls -d /path_to_db > acinetobacter_genome.txt



# load to anvio
mkdir $WORK/pangenome
PAN=$WORK/pangenome





anvi-gen-genomes-storage -e $PAN/acinetobacter_genome.txt -o $PAN/acinetobacter-GENOMES.db

# run pangenomic analysis

anvi-pan-genome -g $PAN/acinetobacter-GENOMES.db -n biomob_acinetobacter -T 24 --align-with famsa -o $PAN/biomob_acinetobacter --enforce-hierarchical-clustering



anvi-compute-genome-similarity -e $PAN/acinetobacter_genome.txt \
-p $PAN/biomob_acinetobacter/biomob_acinetobacter-PAN.db  \
                               -o $PAN/genome-similarity \
                               --program fastANI



anvi-display-pan -p $PAN/biomob_acinetobacter/biomob_acinetobacter-PAN.db -g $PAN/acinetobacter-GENOMES.db

