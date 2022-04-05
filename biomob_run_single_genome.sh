#!/bin/bash
#$ -S /bin/bash


source ~/miniconda3/bin/activate
conda activate anvio7


## CHANGE SAMPLE NAME HERE
## CHANGE SAMPLE NAME HERE
## CHANGE SAMPLE NAME HERE
## CHANGE SAMPLE NAME HERE


GALEN=/isilon/ottawa-rdc/users/shared/chenw_lab/galen
WORK=$GALEN/BioMob
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

#this need to be done on local PC, won't tunnel properly on cluster
#create bins based on taxonomy. LOOK AT PHYLOTREE!

anvi-merge $ANVIO/profile/PROFILE.db -c $ANVIO/CaF467b.db  -C $SAMPLE

#generate new fasta file with contigs associated only with taxa of choice.
anvi-summarize -c $ANVIO/$SAMPLE.db -p $ANVIO/profile/PROFILE.db -C $SAMPLE -o $ANVIO/$SAMPLE

checkm lineage_wf 
$ANVIO/$SAMPLE/bin_by_bin/Helicobacter/Helicobacter-contigs.fa \
$ANVIO/$SAMPLE/checkm \
-x fa -t 20


#extract read for the helicobacter bin
anvi-get-short-reads-from-bam 
-c $ANVIO/$SAMPLE.db \
-p $ANVIO/profile/PROFILE.db \
-C $SAMPLE \
-b Helicobacter \
$MAPPING/${SAMPLE}_anvi.bam \
--gzip-output \
--split-R1-and-R2 \
-O $ANVIO/helicobacter

####################################################################
### RAGTAG
####################################################################
#download best reference genome. For this genome, Helicobacter canadensisn is best

# correct a query assembly
ragtag.py correct $WORK/$SAMPLE/reference_fasta/helicobacter_canadensis/chr.fna \
$ANVIO/$SAMPLE/bin_by_bin/Helicobacter/Helicobacter-contigs.fa \
-o $ANVIO/ragtag_output/

# scaffold a query assembly
ragtag.py scaffold $WORK/$SAMPLE/reference_fasta/helicobacter_canadensis/chr.fna \
$ANVIO/ragtag_output/ragtag.correct.fasta \
-o $ANVIO/ragtag_output/

# make joins and fill gaps in target.fa using sequences from query.fa ## NOT SURE ABOUT THIS STEP.
ragtag.py patch $ANVIO/$SAMPLE/bin_by_bin/Helicobacter/Helicobacter-contigs.fa \
$ANVIO/ragtag_output/ragtag.scaffold.fasta \
-o $ANVIO/ragtag_output/

####################################################################
### phylogeny
####################################################################
conda deactivate
conda activate /isilon/ottawa-rdc/users/shared/chenw_lab/galen/gtdbtk

#create copy of fasta file and change extension so gtdb-tk will only use fasta with new extension

cp $ANVIO/ragtag_output/ragtag.scaffold.fasta $ANVIO/ragtag_output/ragtag.scaffold.fa

gtdbtk de_novo_wf --genome_dir $ANVIO/ragtag_output/ \
--bacteria \
--out_dir $ANVIO/$SAMPLE/gtdbtk_final/ \
--taxa_filter  g__Helicobacter \
--outgroup_taxon f__Arcobacteraceae \
-x fa --cpu 20 


gtdbtk classify_wf --genome_dir $ANVIO/ragtag_output/ \
--out_dir $ANVIO/$SAMPLE/gtdbtk_final \
-x fa --cpu 20 


gtdbtk ani_rep --genome_dir $ANVIO/ragtag_output/ \
--out_dir $ANVIO/$SAMPLE/gtdbtk_final \
-x fa --cpu 20 


####################################################################
### run abricate for AMR and Pathogenic elements
####################################################################

conda activate anvio7

abricate  $ANVIO/ragtag_final/ragtag.scaffold.fa  -db ncbi > $ANVIO/abricate/abricate_ncbi.txt
abricate $ANVIO/ragtag_final/ragtag.scaffold.fa  -db card > $ANVIO/abricate/abricate_card.txt
abricate $ANVIO/ragtag_final/ragtag.scaffold.fa  -db vfdb > $ANVIO/abricate/abricate_vfdb.txt
abricate --summary $ANVIO/abricate/abricate_ncbi.txt $ANVIO/abricate/abricate_card.txt $ANVIO/abricate/abricate_vfdb.txt > $ANVIO/abricate/abricate_summary.txt 

####################################################################
### run genome similarity on anvio
####################################################################

#download on https://www.ncbi.nlm.nih.gov/assembly ref genome 
#search organism (species), then look for full genome only and representative only
#click on big download dataset > genomic sequence (fasta)
#upload to reference_fasta folder in $SAMPLE
mv all fasta to $WORK/$SAMPLE/reference_fasta/

# generate anvi database of reference fasta
mkdir $WORK/$SAMPLE/reference_fasta/anvi_ref/


####for ref in `find $WORK/$SAMPLE/reference_fasta/ncbi-genomes/*.fna  -printf "%f\n"`; 


for ref in `awk '{print $1}' $WORK/$SAMPLE/reference_fasta/ncbi-genomes/ncbi_ref.txt`;

do
echo $ref;
anvi-script-reformat-fasta $WORK/$SAMPLE/reference_fasta/ncbi-genomes/${ref}.fna -o $WORK/$SAMPLE/reference_fasta/anvi_ref/${ref}.fa --simplify-names; 
anvi-gen-contigs-database -f $WORK/$SAMPLE/reference_fasta/anvi_ref/${ref}.fa -o $WORK/$SAMPLE/reference_fasta/anvi_ref/${ref}.db; 
anvi-run-hmms -c $WORK/$SAMPLE/reference_fasta/anvi_ref/${ref}.db -T 20
done


#create file for external genome: needs header (name	contigs_db_path)

ls -d $WORK/$SAMPLE/reference_fasta/anvi_ref/*.db > $WORK/$SAMPLE/reference_fasta/anvi_ref/ref_list_db.txt
#upload  ref_list_db.txt on excel and add taxa names

#create internal genome with this as header (put in the correct info)
#name	bin_id	collection_id	profile_db_path	contigs_db_path

anvi-gen-genomes-storage -e $WORK/$SAMPLE/EXTERNAL-GENOMES.txt -i $WORK/$SAMPLE/internal-genomes.txt -o $WORK/$SAMPLE/HELICOBACTER-GENOMES.db

anvi-pan-genome -g $WORK/$SAMPLE/HELICOBACTER-GENOMES.db -n helicobacter --output-dir $WORK/$SAMPLE/helicobacter

anvi-compute-genome-similarity -e $WORK/$SAMPLE/external-genomes.txt -i $WORK/$SAMPLE/internal-genomes.txt --program pyANI \
                               --output-dir $WORK/$SAMPLE/helicobacter/ANI \
                               --num-threads 20 \
                               --pan-db $WORK/$SAMPLE/helicobacter/helicobacter-PAN.db


anvi-get-sequences-for-hmm-hits -e $WORK/$SAMPLE/external-genomes.txt  \
-i $WORK/$SAMPLE/internal-genomes.txt \
                                   --hmm-source Bacteria_71 \
                                   --list-available-gene-names


#create phylogenomic tree
anvi-get-sequences-for-hmm-hits -e $WORK/$SAMPLE/external-genomes.txt \
-i $WORK/$SAMPLE/internal-genomes.txt \
 -o $WORK/$SAMPLE/helicobacter/phylogeny/concatenated-genes.fa \
--hmm-sources Bacteria_71 \
--concatenate \
--align-with famsa \
--get-aa-sequences \
--return-best-hit \
--just-do-it 


trimal -in $WORK/$SAMPLE/helicobacter/phylogeny/concatenated-align-genes-fasta \
       -out $WORK/$SAMPLE/helicobacter/phylogeny/align-genes-fasta_REMOVED.fa \
       -gt 0.50

anvi-gen-phylogenomic-tree -f $WORK/$SAMPLE/helicobacter/phylogeny/concatenated-align-genes-fasta -o $WORK/$SAMPLE/helicobacter/phylogeny/phylogeny



#### trimming not used

adapter=/home/AAFC-AAC/guog/miniconda3/envs/anvio7/opt/bbmap-38.18/resources

bbduk.sh -Xmx1g \
in1=$RAW/${SAMPLE}_*_R1_001.fastq.gz \
in2=$RAW/${SAMPLE}_*_R2_001.fastq.gz \
out1=$TRIM/${SAMPLE}.clean1.fq \
out2=$TRIM/${SAMPLE}.clean2.fq \
ref=$adapter/adapters.fa k=31 hdist=1 \
stats=$TRIM/${SAMPLE}_stats.txt
	
java -jar $GALEN/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-threads 20 $TRIM/${SAMPLE}.clean1.fq $TRIM/${SAMPLE}.clean2.fq \
$TRIM/${SAMPLE}.pair1.fq.gz $TRIM/${SAMPLE}.unpair1.fq.gz \
$TRIM/${SAMPLE}.pair2.fq.gz $TRIM/${SAMPLE}.unpair2.fq.gz \
ILLUMINACLIP:$GALEN/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 


