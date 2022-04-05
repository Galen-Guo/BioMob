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
RAW=$WORK/raw_data

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

####################################################################
### checkm
####################################################################

for SAMPLE in `awk '{print $1}' $RAW/sample_name`;
do
echo $SAMPLE

#directories


TRIM=$WORK/$SAMPLE/trimming_clean
ASSEMBLY=$WORK/$SAMPLE/assembly
ANVIO=$WORK/$SAMPLE/anvio/
MAPPING=$WORK/$SAMPLE/mapping/


checkm lineage_wf \
$ANVIO/$SAMPLE/bin_by_bin/EVERYTHING/ \
$ANVIO/$SAMPLE/checkm \
-x fa -t 16
done

## concatenate all "bin_summary" file into one. 


for SAMPLE in `awk '{print $1}' $RAW/sample_name`;
do
echo $SAMPLE

#directories


TRIM=$WORK/$SAMPLE/trimming_clean
ASSEMBLY=$WORK/$SAMPLE/assembly
ANVIO=$WORK/$SAMPLE/anvio/
MAPPING=$WORK/$SAMPLE/mapping/

echo $SAMPLE >> all_acinetobacter_summary.txt
cat $ANVIO/$SAMPLE/bins_summary.txt >> all_acinetobacter_summary.txt

done


####################################################################
### loading into anvio ### 
####################################################################

find $RAW/*.fasta  -printf "%f\n" > sample_name
sed -e s/.fasta//g -i $RAW/sample_name

for i in `awk '{print $1}' $RAW/sample_name`;
do
echo $i
anvi-script-reformat-fasta $RAW/${i}.fasta -o $ANVIO/${i}.fasta  -l 5000 --simplify-names -r $ANVIO/${i}_report.txt
done





###create anvio database

mkdir $ANVIO

#Gene calling, taxonomical, functional assignment to contigs and generate summary

for i in `awk '{print $1}' $RAW/sample_name`;
do

echo $i
mkdir $ANVIO/$i
cd $ANVIO/$i
anvi-gen-contigs-database -f $ANVIO/$i/${i}.fasta -o $ANVIO/$i/${i}.db -n GRDI -T 80
anvi-run-hmms -c $ANVIO/$i/${i}.db -T 80
anvi-run-kegg-kofams -c $ANVIO/$i/${i}.db -T 80 --hmmer-program hmmsearch --keep-all-hits --just-do-it
anvi-run-ncbi-cogs -c $ANVIO/$i/${i}.db --cog-data-dir $ANVIO/$i/COG --num-threads 80 --search-with diamond
anvi-get-sequences-for-gene-calls -c $ANVIO/$i/${i}.db --get-aa-sequences -o -c $ANVIO/$i/${i}_aa.fa
anvi-run-scg-taxonomy -c $ANVIO/$i/${i}.db -T 80 
anvi-estimate-scg-taxonomy -c $ANVIO/$i/${i}.db --output-file $ANVIO/$i/${i}_scg_est_output.txt
anvi-display-contigs-stats $ANVIO/$i/${i}.db  --report-as-text  --output-file $ANVIO/$i/${i}_contigs_summary.txt
done




### pangenomic?

#create external file

name	contigs_db_path
Name_01	/path/to/contigs-01.db
Name_02	/path/to/contigs-02.db
Name_03	/path/to/contigs-03.db
(…)	(…)


# load to anvio
anvi-gen-genomes-storage -e $ANVIO/external_genome.txt -o $ANVIO/aqua_pathos-GENOMES.db

# run pangenomic analysis
anvi-pan-genome -g $ANVIO/aqua_pathos-GENOMES.db -n BioMob

anvi-display-pan -p $ANVIO/BioMob/BioMob-PAN.db -g $ANVIO/aqua_pathos-GENOMES.db

