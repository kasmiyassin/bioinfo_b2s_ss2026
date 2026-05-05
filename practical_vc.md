
* info ref genome at ``` /courses/master_b2s/references/GRCh/hg19/hg19.chr5_12_17.fa.gz``` 
* working data at ```/courses/master_b2s/raw_data/s02_vc/cancer```


1. mkdir for working 

```bash
cd 
mkdir -p session02vc/{scripts,results,figures,reports,alignments,vcf_files} 
tree -d session02vc/

cd session02vc
```

2. fastqc

```bash
## version and availiblity

fastqc -v
## create working direcotr
mkdir -p ./reports/fastqc/

# fastqc
fastqc -o ./reports/fastqc/ \
        --threads 1 \
        /courses/master_b2s/raw_data/s02_vc/cancer/*.gz

# see file
ls -l ./reports/fastqc/
```
now open WinSCP and download and open html in your browser

3. MultiQC
```bash
## version and availiblity
multiqc -v
## create working direcotr
mkdir -p ./reports/multiqc_step01/
# cd ./reports/multiqc_step01/
# run analyse
multiqc ./reports/fastqc/ -o ./reports/multiqc_step01/

# see file
ls -l ./reports/multiqc_step01/
```

4. trimming
```bash
# 1. Create the output directory
mkdir -p /courses/master_b2s/raw_data/s02_vc/cancer_trim

/courses/master_b2s/raw_data/s02_vc/cancer
# Run Trimmomatic on the Normal (N) sample
java -jar /courses/master_b2s/software/Trimmomatic/trimmomatic-0.40.jar PE -threads 4 \
    /courses/master_b2s/raw_data/s02_vc/cancer/normal_r1.fastq.gz \
    /courses/master_b2s/raw_data/s02_vc/cancer/normal_r2.fastq.gz \
    /courses/master_b2s/raw_data/s02_vc/cancer_trim/normal_paired_R1.fastq.gz \
    /courses/master_b2s/raw_data/s02_vc/cancer_trim/normal_unpaired_R1.fastq.gz \
    /courses/master_b2s/raw_data/s02_vc/cancer_trim/normal_paired_R2.fastq.gz \
    /courses/master_b2s/raw_data/s02_vc/cancer_trim/normal_unpaired_R2.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads \
    HEADCROP:3 \
    TRAILING:10 \
    MINLEN:25



# Run Trimmomatic on the Tumor (T) sample
java -jar /courses/master_b2s/software/Trimmomatic/trimmomatic-0.40.jar PE -threads 4 \
    /courses/master_b2s/raw_data/s02_vc/cancer/tumor_r1.fastq.gz \
    /courses/master_b2s/raw_data/s02_vc/cancer/tumor_r2.fastq.gz \
    /courses/master_b2s/raw_data/s02_vc/cancer_trim/tumor_paired_R1.fastq.gz /courses/master_b2s/raw_data/s02_vc/cancer_trim/tumor_unpaired_R1.fastq.gz \
    /courses/master_b2s/raw_data/s02_vc/cancer_trim/tumor_paired_R2.fastq.gz /courses/master_b2s/raw_data/s02_vc/cancer_trim/tumor_unpaired_R2.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads \
    HEADCROP:3 \
    TRAILING:10 \
    MINLEN:25
```
5. repeate QC

```bash
# fastqc
cd 
cd session02vc
fastqc -o ./reports/fastqc/ \
        --threads 1 \
        /courses/master_b2s/raw_data/s02_vc/cancer_trim/*.gz

# see file
ls -l ./reports/fastqc/

# now open WinSCP and download and open html in your browser
#MultiQC
cd 
cd session02vc
mkdir ./reports/multiqc_step02/
# run analyse
multiqc ./reports/fastqc/ -o ./reports/multiqc_step02/ -t 10

# see file
ls -l ./reports/multiqc_step02/
```




```bash
mkdir -p fastqc/original fastqc/trimmed trimm_out multiqc

for R1 in *_r1.fastq.gz; do

    R2=${R1/_R1/_R2}

    BASE=$(basename $R1 _r1.fastq.gz)

    echo "Processing $BASE..."

    java -jar /courses/master_b2s/software/Trimmomatic/trimmomatic-0.40.jar PE -threads 4 $R1 $R2 trimm_out/${BASE}_R1_paired.fastq.gz trimm_out/${BASE}_R1_unpaired.fastq.gz trimm_out/${BASE}_R2_paired.fastq.gz trimm_out/${BASE}_R2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> trimm_out/${BASE}.log

done



fastqc *_r1.fastq.gz *_r2.fastq.gz -o fastqc/original/

fastqc trimm_out/*_paired.fastq.gz -o fastqc/trimmed/



multiqc . -o multiqc/

```




--- 

6. index reference genome

```bash
mkdir /courses/master_b2s/references/GRCh/hg19
cd /courses/master_b2s/references/GRCh/hg19
bwa index hg19.chr5_12_17.fa.gz

```

7. alignemnt bam

```bash

REFERENCE_SEQUENCE=/courses/master_b2s/references/GRCh/hg19/hg19.chr5_12_17.fa.gz
LEFT_READS=/courses/master_b2s/raw_data/s02_vc/cancer/normal_r1.fastq.gz
RIGHT_READS=/courses/master_b2s/raw_data/s02_vc/cancer/normal_r2.fastq.gz
SAMPLE=`basename $LEFT_READS _r1.fastq.gz`
#SAMPLE=normal
SAM_FILE=~/session02vc/alignments/${SAMPLE}_hg19.sam

# Align reads with bwa
bwa mem \
    -M \
    -t 1 \
    -R "@RG\tID:$SAMPLE\tPL:illumina\tPU:$SAMPLE\tSM:$SAMPLE" \
    $REFERENCE_SEQUENCE \
    $LEFT_READS \
    $RIGHT_READS \
    -o $SAM_FILE

```

for tumor too


```bash
REFERENCE_SEQUENCE=/courses/master_b2s/references/GRCh/hg19/hg19.chr5_12_17.fa.gz
LEFT_READS=/courses/master_b2s/raw_data/s02_vc/cancer/tumor_r1.fastq.gz
RIGHT_READS=/courses/master_b2s/raw_data/s02_vc/cancer/tumor_r2.fastq.gz
#SAMPLE=tumor
SAMPLE=`basename $LEFT_READS _r1.fastq.gz`
SAM_FILE=~/session02vc/alignments/${SAMPLE}_hg19.sam

# Align reads with bwa
bwa mem \
    -M \
    -t 1 \
    -R "@RG\tID:$SAMPLE\tPL:illumina\tPU:$SAMPLE\tSM:$SAMPLE" \
    $REFERENCE_SEQUENCE \
    $LEFT_READS \
    $RIGHT_READS \
    -o $SAM_FILE
```

8. Alignment file processing


```bash

# Assign file paths to variables
SAMPLE_NAME=normal
SAM_FILE=~/session02vc/alignments/${SAMPLE_NAME}_hg19.sam
REPORTS_DIRECTORY=~/session02vc/reports/picard/${SAMPLE_NAME}/
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
REMOVE_DUPLICATES_BAM_FILE=`echo ${SAM_FILE%sam}remove_duplicates.bam`
METRICS_FILE=${REPORTS_DIRECTORY}/${SAMPLE_NAME}.remove_duplicates_metrics.txt
COORDINATE_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}coordinate_sorted.bam`

# Make reports directory
mkdir -p $REPORTS_DIRECTORY

# Query-sort alginment file and convert to BAM
#java -jar $PICARD/picard.jar SortSam \
picard SortSam \
  --INPUT $SAM_FILE \
  --OUTPUT $QUERY_SORTED_BAM_FILE \
  --SORT_ORDER queryname

# Mark and remove duplicates
# java -jar $PICARD/picard.jar
picard MarkDuplicates \
  --INPUT $QUERY_SORTED_BAM_FILE \
  --OUTPUT $REMOVE_DUPLICATES_BAM_FILE \
  --METRICS_FILE $METRICS_FILE \
  --REMOVE_DUPLICATES true

# Coordinate-sort BAM file and create BAM index file
#java -jar $PICARD/picard.jar 
picard SortSam \
  --INPUT $REMOVE_DUPLICATES_BAM_FILE \
  --OUTPUT $COORDINATE_SORTED_BAM_FILE \
  --SORT_ORDER coordinate \
  --CREATE_INDEX true

```
we can add also read group

```bash
 # Add or replace read group information
 #java -jar $PICARD/picard.jar 
 picard AddOrReplaceReadGroups \
   --INPUT $SAM_FILE \
   --OUTPUT $BAM_FILE \
   --RGID $READ_GROUP_ID \
   --RGLB $READ_GROUP_LIBRARY \
   --RGPL $READ_GROUP_PLATFORM \
   --RGPU $READ_GROUP_PLATFORM_UNIT \
   --RGSM $READ_GROUP_SAMPLE
```

We do for tumoral

```bash
# Assign file paths to variables
SAMPLE_NAME=tumor
SAM_FILE=~/session02vc/alignments/${SAMPLE_NAME}_hg19.sam
REPORTS_DIRECTORY=~/session02vc/reports/picard/${SAMPLE_NAME}/
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
REMOVE_DUPLICATES_BAM_FILE=`echo ${SAM_FILE%sam}remove_duplicates.bam`
METRICS_FILE=${REPORTS_DIRECTORY}/${SAMPLE_NAME}.remove_duplicates_metrics.txt
COORDINATE_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}coordinate_sorted.bam`
# Make reports directory
mkdir -p $REPORTS_DIRECTORY

# Query-sort alginment file and convert to BAM
#java -jar $PICARD/picard.jar SortSam \
picard SortSam \
  --INPUT $SAM_FILE \
  --OUTPUT $QUERY_SORTED_BAM_FILE \
  --SORT_ORDER queryname

# Mark and remove duplicates
# java -jar $PICARD/picard.jar
picard MarkDuplicates \
  --INPUT $QUERY_SORTED_BAM_FILE \
  --OUTPUT $REMOVE_DUPLICATES_BAM_FILE \
  --METRICS_FILE $METRICS_FILE \
  --REMOVE_DUPLICATES true

# Coordinate-sort BAM file and create BAM index file
#java -jar $PICARD/picard.jar 
picard SortSam \
  --INPUT $REMOVE_DUPLICATES_BAM_FILE \
  --OUTPUT $COORDINATE_SORTED_BAM_FILE \
  --SORT_ORDER coordinate \
  --CREATE_INDEX true

```



10. alignment_QC.m

if for some reason you could not do this work. You can copie the files

```bash
 cp /courses/master_b2s/session02vc/alignments/*_hg19.coordinate_sorted.bam session2vc/alignments/

```


```bash
# Assign variables
INPUT_BAM=~/session02vc/alignments/normal_hg19.coordinate_sorted.bam
REFERENCE=/courses/master_b2s/references/GRCh/hg19/hg19.chr5_12_17.fa.gz
OUTPUT_METRICS_FILE=~/session02vc/reports/picard/normal/normal_hg19.CollectAlignmentSummaryMetrics.txt

# Run Picard CollectAlignmentSummaryMetrics
#java -jar $PICARD/picard.jar 
picard CollectAlignmentSummaryMetrics \
  --INPUT $INPUT_BAM \
  --REFERENCE_SEQUENCE $REFERENCE \
  --OUTPUT $OUTPUT_METRICS_FILE


# Assign variables

INPUT_BAM=~/session02vc/alignments/tumor_hg19.coordinate_sorted.bam
REFERENCE=/courses/master_b2s/references/GRCh/hg19/hg19.chr5_12_17.fa.gz
OUTPUT_METRICS_FILE=~/session02vc/reports/picard/tumor/tumor_hg19.CollectAlignmentSummaryMetrics.txt

# Run Picard CollectAlignmentSummaryMetrics
#java -jar $PICARD/picard.jar 
picard CollectAlignmentSummaryMetrics \
  --INPUT $INPUT_BAM \
  --REFERENCE_SEQUENCE $REFERENCE \
  --OUTPUT $OUTPUT_METRICS_FILE

```

11. Running MultiQC Aggregating QC


```bash
less ~/session02vc/reports/picard/normal/normal_hg19.CollectAlignmentSummaryMetrics.txt
```


```bash
tree ~/session02vc/reports/
PLACEHOLDER FOR COMMAND

```


```bash

# Assign variables
REPORTS_DIRECTORY=~/session02vc/reports/
NORMAL_SAMPLE_NAME=normal
TUMOR_SAMPLE_NAME=tumor
REFERENCE=hg19
NORMAL_PICARD_METRICS=${REPORTS_DIRECTORY}picard/${NORMAL_SAMPLE_NAME}/${NORMAL_SAMPLE_NAME}_${REFERENCE}.CollectAlignmentSummaryMetrics.txt
TUMOR_PICARD_METRICS=${REPORTS_DIRECTORY}picard/${TUMOR_SAMPLE_NAME}/${TUMOR_SAMPLE_NAME}_${REFERENCE}.CollectAlignmentSummaryMetrics.txt
NORMAL_FASTQC_1=${REPORTS_DIRECTORY}fastqc/${NORMAL_SAMPLE_NAME}_r1_fastqc.zip
NORMAL_FASTQC_2=${REPORTS_DIRECTORY}fastqc/${NORMAL_SAMPLE_NAME}_r2_fastqc.zip
TUMOR_FASTQC_1=${REPORTS_DIRECTORY}fastqc/${TUMOR_SAMPLE_NAME}_r1_fastqc.zip
TUMOR_FASTQC_2=${REPORTS_DIRECTORY}fastqc/${TUMOR_SAMPLE_NAME}_r2_fastqc.zip
OUTPUT_DIRECTORY=${REPORTS_DIRECTORY}/multiqc/

# Create directory for output
mkdir -p $OUTPUT_DIRECTORY

# Run MultiQC
multiqc \
  $NORMAL_PICARD_METRICS \
  $TUMOR_PICARD_METRICS \
  $NORMAL_FASTQC_1 \
  $NORMAL_FASTQC_2 \
  $TUMOR_FASTQC_1 \
  $TUMOR_FASTQC_2 \
  --outdir $OUTPUT_DIRECTORY

```

to solve issue here you need to change name in zip to correct name **normal_1_fastqc**

12. variant_calling


```bash
conda activate variant-analysis
# Verify it now says version 17 or 21
java -version

```
create diction GATK of ref genome
```bash
cd /courses/master_b2s/references/GRCh/hg19/
gatk CreateSequenceDictionary -R hg19.chr5_12_17.fa

# YOU DON'T NEED TO RUN THIS
samtools faidx \
  reference_sequence.fa

# YOU DON'T NEED TO RUN THIS
 java -jar $PICARD/picard.jar CreateSequenceDictionary \
 --REFERENCE /n/groups/hbctraining/variant_calling/reference/GRCh37.p13_genomic.fa
 --OUTPUT /n/groups/hbctraining/variant_calling/reference/GRCh37.p13_genomic.dict


```

```bash
# Assign variables
REFERENCE_SEQUENCE=/courses/master_b2s/references/GRCh/hg19/hg19.chr5_12_17.fa
REFERENCE_DICTIONARY=`echo ${REFERENCE_SEQUENCE%fa}dict`
NORMAL_SAMPLE_NAME=normal
NORMAL_BAM_FILE=~/session02vc/alignments/${NORMAL_SAMPLE_NAME}_hg19.coordinate_sorted.bam
TUMOR_SAMPLE_NAME=tumor
TUMOR_BAM_FILE=~/session02vc/alignments/${TUMOR_SAMPLE_NAME}_hg19.coordinate_sorted.bam
VCF_OUTPUT_FILE=~/session02vc/vcf_files/mutect2_${NORMAL_SAMPLE_NAME}_${TUMOR_SAMPLE_NAME}_hg19-raw.vcf

# Run MuTect2 toooooooooooo lang
gatk Mutect2 \
  --sequence-dictionary $REFERENCE_DICTIONARY \
  --reference $REFERENCE_SEQUENCE \
  --input $NORMAL_BAM_FILE \
  --normal-sample $NORMAL_SAMPLE_NAME \
  --input $TUMOR_BAM_FILE \
  --tumor-sample $TUMOR_SAMPLE_NAME \
  --annotation ClippingRankSumTest --annotation DepthPerSampleHC --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation FisherStrand --annotation MappingQuality --annotation DepthPerAlleleBySample --annotation Coverage \
  --output $VCF_OUTPUT_FILE


## fast 
REFERENCE_SEQUENCE=/courses/master_b2s/references/GRCh/hg19/hg19.chr5_12_17.fa
REFERENCE_DICTIONARY=/courses/master_b2s/references/GRCh/hg19/hg19.chr5_12_17.dict

NORMAL_SAMPLE_NAME=normal
NORMAL_BAM_FILE=~/session02vc/alignments/${NORMAL_SAMPLE_NAME}_hg19.coordinate_sorted.bam

TUMOR_SAMPLE_NAME=tumor
TUMOR_BAM_FILE=~/session02vc/alignments/${TUMOR_SAMPLE_NAME}_hg19.coordinate_sorted.bam

VCF_OUTPUT_FILE=~/session02vc/vcf_files/mutect2_somatic_turbo.vcf

gatk Mutect2 \
  -R $REFERENCE_SEQUENCE \
  -L chr5 -L chr12 -L chr17 \
  -I $NORMAL_BAM_FILE -normal $NORMAL_SAMPLE_NAME \
  -I $TUMOR_BAM_FILE -tumor $TUMOR_SAMPLE_NAME \
  --native-pair-hmm-threads 4 \
  --annotation ClippingRankSumTest --annotation DepthPerSampleHC \
  --annotation MappingQualityRankSumTest --annotation MappingQualityZero \
  --annotation QualByDepth --annotation ReadPosRankSumTest \
  --annotation RMSMappingQuality --annotation FisherStrand \
  --annotation MappingQuality --annotation DepthPerAlleleBySample \
  --annotation Coverage \
  -O $VCF_OUTPUT_FILE

```


13. variant_filtering

```bash 

# Assign variables
REFERENCE_SEQUENCE=/courses/master_b2s/references/GRCh/hg19/hg19.chr5_12_17.fa
RAW_VCF_FILE=~/session02vc/vcf_files/mutect2_normal_tumor_hg19-raw.vcf
LCR_FILE=/courses/master_b2s/references/GRCh/hg19/LCR-hs37_with_chr.bed
MUTECT_FILTERED_VCF=${RAW_VCF_FILE%raw.vcf}filt.vcf
PASSING_FILTER_VCF=${RAW_VCF_FILE%raw.vcf}pass-filt.vcf
LCR_FILTERED_VCF=${RAW_VCF_FILE%raw.vcf}pass-filt-LCR.vcf

# Filter Mutect Calls
gatk FilterMutectCalls \
  --reference $REFERENCE_SEQUENCE \
  --variant $RAW_VCF_FILE \
  --output $MUTECT_FILTERED_VCF

# Filter for only SNPs with PASS in the FILTER field
/courses/master_b2s/software/miniconda3/envs/snpeff-env/bin/java -jar /courses/master_b2s/software/snpEff/SnpSift.jar filter \
  -noLog \
  "( FILTER = 'PASS' )" \
  $MUTECT_FILTERED_VCF > $PASSING_FILTER_VCF

# Filter LCR
/courses/master_b2s/software/miniconda3/envs/snpeff-env/bin/java -jar /courses/master_b2s/software/snpEff/SnpSift.jar intervals \
  -noLog \
  -x \
  -i $PASSING_FILTER_VCF \
  $LCR_FILE > $LCR_FILTERED_VCF
  

###########

less ~/session02vc/vcf_files/mutect2_normal_tumor_hg19-pass-filt-LCR.vcf

grep -v "^#" ~/session02vc/vcf_files/mutect2_normal_tumor_hg19-pass-filt-LCR.vcf | wc -l

```


14. variant annotation
 
```bash
cd ~/session02vc/scripts/
nano ~/session02vc/scripts/normal_tumor_pedigree_header.txt
##PEDIGREE=<Derived=tumor,Original=normal>

# Assign variables
REPORTS_DIRECTORY=~/session02vc/reports/snpeff/
SAMPLE_NAME=mutect2_normal_tumor
REFERENCE_SEQUENCE_NAME=hg19
CSV_STATS=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.csv"`
HTML_REPORT=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.html"`
#REFERENCE_DATABASE=hg19
REFERENCE_DATABASE=GRCh37.p13.RefSeq
DATADIR=/courses/master_b2s/references/GRCh/hg38/snpeff/data/
FILTERED_VCF_FILE=~/session02vc/vcf_files/mutect2_normal_tumor_hg19-pass-filt-LCR.vcf
PEDIGREE_HEADER_FILE=~/session02vc/scripts/normal_tumor_pedigree_header.txt
FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER=${FILTERED_VCF_FILE%.vcf}.pedigree_header.vcf
SNPEFF_ANNOTATED_VCF_FILE=${FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER%.vcf}.snpeff.vcf
DBSNP_DATABASE=/courses/master_b2s/references/GRCh/hg19/GRCh37.p13.dbSNP.vcf.gz
DBSNP_ANNOTATED_VCF_FILE=${SNPEFF_ANNOTATED_VCF_FILE%.vcf}.dbSNP.vcf

# Create reports directory
mkdir -p $REPORTS_DIRECTORY

# Append Header
bcftools annotate \
  --header-lines $PEDIGREE_HEADER_FILE \
  $FILTERED_VCF_FILE \
  > $FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER

# Run SnpEff
/courses/master_b2s/software/miniconda3/envs/snpeff-env/bin/java -jar /courses/master_b2s/software/snpEff/snpEff.jar databases | less

/courses/master_b2s/software/miniconda3/envs/snpeff-env/bin/java -jar -Xmx4g /courses/master_b2s/software/snpEff/snpEff.jar  eff \
  -dataDir $DATADIR \
  -cancer \
  -noLog \
  -csvStats $CSV_STATS \
  -s $HTML_REPORT \
  $REFERENCE_DATABASE \
  $FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER \
  > $SNPEFF_ANNOTATED_VCF_FILE

# Use dbSNP VCF to annotate our VCF
/courses/master_b2s/software/miniconda3/envs/snpeff-env/bin/java -jar /courses/master_b2s/software/snpEff/SnpSift.jar annotate \
  $DBSNP_DATABASE \
  -tabix \
  -noLog \
  $SNPEFF_ANNOTATED_VCF_FILE \
  > $DBSNP_ANNOTATED_VCF_FILE


less  ~/session02vc/vcf_files/mutect2_normal_tumor_hg19-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf 

##SnpEffCmd="SnpEff  -cancer -csvStats /root/session02vc/reports/snpeff/annotation_mutect2_normal_tumor_hg19-effects-stats.csv -s /ro>
```

## install database

```bash

curl -o GRCh37.p13.dbSNP.vcf.gz -L https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz

tabix GRCh37.p13.dbSNP.vcf.gz

# Use dbSNP VCF to annotate our VCF
java -jar $SNPEFF/SnpSift.jar annotate \
  $DBSNP_DATABASE \
  -tabix \
  -noLog \
  $SNPEFF_ANNOTATED_VCF_FILE \
  > $DBSNP_ANNOTATED_VCF_FILE
```



15. Prioritizing Variants

```bash

cd ~/session02vc/vcf_files/

/courses/master_b2s/software/miniconda3/envs/snpeff-env/bin/java -jar /courses/master_b2s/software/snpEff/SnpSift.jar filter \
  -noLog \
  "( CHROM = '5' )" \
  mutect2_normal_tumor_hg19-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf  | less
# in our data we need to use chr5 instead 5

/courses/master_b2s/software/miniconda3/envs/snpeff-env/bin/java -jar /courses/master_b2s/software/snpEff/SnpSift.jar filter \
  "( CHROM = 'chr1' ) | ( CHROM = 'chr5' )" \
  mutect2_normal_tumor_hg19-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf  | less


# Alternatively, we could be interested in variants on Chromosome 1 between positions 1000000 and 2000000. This command would look like:

/courses/master_b2s/software/miniconda3/envs/snpeff-env/bin/java -jar /courses/master_b2s/software/snpEff/SnpSift.jar filter \
  -noLog \
  "( CHROM = '1' ) & ( POS > 1000000 ) & ( POS < 2000000 )" \
  mutect2_normal_tumor_hg19-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf   | less


# If you are interested in all of the variants corresponding to a single gene of interest, you can filter by the gene name in this case CPSF3L:

/courses/master_b2s/software/miniconda3/envs/snpeff-env/bin/java -jar /courses/master_b2s/software/snpEff/SnpSift.jar filter \
  -noLog \
  "( ANN[*].GENE = 'SDHAP3' )" mutect2_normal_tumor_hg19-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf | less



# You can also filter on the transcript ID which in this case is the NCBI accession number.

/courses/master_b2s/software/miniconda3/envs/snpeff-env/bin/java -jar /courses/master_b2s/software/snpEff/SnpSift.jar filter \
  -noLog \
  "( ANN[*].TRID = 'XM_017001557.1' )" mutect2_normal_tumor_hg19-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf | less


# If you want to filter your output by the effects the variants have on the annotated gene models, the syntax for this is quite similar to the example for genes:

/courses/master_b2s/software/miniconda3/envs/snpeff-env/bin/java -jar /courses/master_b2s/software/snpEff/SnpSift.jar filter \
  -noLog \
  "( ANN[*].EFFECT has 'missense_variant' )" \
  mutect2_normal_tumor_hg19-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf  | less

```


```bash

#Let's go ahead and select out all of our HIGH impact muations:

/courses/master_b2s/software/miniconda3/envs/snpeff-env/bin/java -jar /courses/master_b2s/software/snpEff/SnpSift.jar filter \
  -noLog \
  "( ANN[*].IMPACT has 'HIGH' )" \
  mutect2_normal_tumor_hg19-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf  | less


#Let's go ahead and redirect the output of these "high-impact" mutations to a new VCF file:

/courses/master_b2s/software/miniconda3/envs/snpeff-env/bin/java -jar /courses/master_b2s/software/snpEff/SnpSift.jar filter \
  -noLog \
  "( ANN[*].IMPACT has 'HIGH' )"  \
  mutect2_normal_tumor_hg19-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf  > mutect2_normal_tumor_hg19-pass-filt-LCR.pedigree_header.snpeff.dbSNP.high_impact.vcf 



# A useful tool within the SnpSift toolkit is the perl script named vcfEffOnePerLine.pl. This script allows the user to separate each effect onto its own line instead of having them lumped into a single line. In order to utilize this script we need to pipe the output of our filter command into $SNPEFF/scripts/vcfEffOnePerLine.pl. We can use it on our previous example to demonstrate:
# one time chmod +x /courses/software/vcfInfoOnePerLine.pl
/courses/master_b2s/software/miniconda3/envs/snpeff-env/bin/java -jar /courses/master_b2s/software/snpEff/SnpSift.jar filter \
  -noLog \
  "( ANN[*].IMPACT has 'HIGH' )"  \
  mutect2_normal_tumor_hg19-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf | \
  /courses/master_b2s/software/snpEff/scripts/vcfEffOnePerLine.pl | less

```

```bash
1	6471577	.	A	ACTCACGTGCAAGCATCACACCGGCACGC	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=49,52|6,4;ClippingRankSum=-1.498;DP=119;ECNT=1;FS=2.779;GERMQ=93;MBQ=32,32;MFRL=341,338;MMQ=60,60;MPOS=20;MQ=60;MQ0=0;MQRankSum=0;NALOD=1.82;NLOD=19.17;POPAF=6;ReadPosRankSum=-1.659;TLOD=33.37;LOF=(PLEKHG5|PLEKHG5|8|1.00);NMD=(PLEKHG5|PLEKHG5|8|1.00);ANN=ACTCACGTGCAAGCATCACACCGGCACGC|frameshift_variant&stop_gained|HIGH|PLEKHG5|PLEKHG5|transcript|NM_001265592.1|protein_coding|13/22|c.1428_1429insGCGTGCCGGTGTGATGCTTGCACGTGAG|p.Trp477fs|1493/4794|1428/3258|476/1085||	GT:AD:AF:DP:F1R2:F2R1:SB	0/0:63,0:0.015:63:37,0:25,0:33,30,0,0	0/1:38,10:0.216:48:19,5:17,4:16,22,6,4
1	6471577	.	A	ACTCACGTGCAAGCATCACACCGGCACGC	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=49,52|6,4;ClippingRankSum=-1.498;DP=119;ECNT=1;FS=2.779;GERMQ=93;MBQ=32,32;MFRL=341,338;MMQ=60,60;MPOS=20;MQ=60;MQ0=0;MQRankSum=0;NALOD=1.82;NLOD=19.17;POPAF=6;ReadPosRankSum=-1.659;TLOD=33.37;LOF=(PLEKHG5|PLEKHG5|8|1.00);NMD=(PLEKHG5|PLEKHG5|8|1.00);ANN=ACTCACGTGCAAGCATCACACCGGCACGC|frameshift_variant&stop_gained|HIGH|PLEKHG5|PLEKHG5|transcript|NM_001265593.1|protein_coding|12/21|c.1398_1399insGCGTGCCGGTGTGATGCTTGCACGTGAG|p.Trp467fs|1424/4725|1398/3228|466/1075||	GT:AD:AF:DP:F1R2:F2R1:SB	0/0:63,0:0.015:63:37,0:25,0:33,30,0,0	0/1:38,10:0.216:48:19,5:17,4:16,22,6,4
1	6471577	.	A	ACTCACGTGCAAGCATCACACCGGCACGC	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=49,52|6,4;ClippingRankSum=-1.498;DP=119;ECNT=1;FS=2.779;GERMQ=93;MBQ=32,32;MFRL=341,338;MMQ=60,60;MPOS=20;MQ=60;MQ0=0;MQRankSum=0;NALOD=1.82;NLOD=19.17;POPAF=6;ReadPosRankSum=-1.659;TLOD=33.37;LOF=(PLEKHG5|PLEKHG5|8|1.00);NMD=(PLEKHG5|PLEKHG5|8|1.00);ANN=ACTCACGTGCAAGCATCACACCGGCACGC|frameshift_variant&stop_gained|HIGH|PLEKHG5|PLEKHG5|transcript|NM_001042665.1|protein_coding|12/21|c.1191_1192insGCGTGCCGGTGTGATGCTTGCACGTGAG|p.Trp398fs|1397/4698|1191/3021|397/1006||	GT:AD:AF:DP:F1R2:F2R1:SB	0/0:63,0:0.015:63:37,0:25,0:33,30,0,0	0/1:38,10:0.216:48:19,5:17,4:16,22,6,4
1	6471577	.	A	ACTCACGTGCAAGCATCACACCGGCACGC	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=49,52|6,4;ClippingRankSum=-1.498;DP=119;ECNT=1;FS=2.779;GERMQ=93;MBQ=32,32;MFRL=341,338;MMQ=60,60;MPOS=20;MQ=60;MQ0=0;MQRankSum=0;NALOD=1.82;NLOD=19.17;POPAF=6;ReadPosRankSum=-1.659;TLOD=33.37;LOF=(PLEKHG5|PLEKHG5|8|1.00);NMD=(PLEKHG5|PLEKHG5|8|1.00);ANN=ACTCACGTGCAAGCATCACACCGGCACGC|frameshift_variant&stop_gained|HIGH|PLEKHG5|PLEKHG5|transcript|NM_001042664.1|protein_coding|12/21|c.1191_1192insGCGTGCCGGTGTGATGCTTGCACGTGAG|p.Trp398fs|1414/4715|1191/3021|397/1006||	GT:AD:AF:DP:F1R2:F2R1:SB	0/0:63,0:0.015:63:37,0:25,0:33,30,0,0	0/1:38,10:0.216:48:19,5:17,4:16,22,6,4
1	6471577	.	A	ACTCACGTGCAAGCATCACACCGGCACGC	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=49,52|6,4;ClippingRankSum=-1.498;DP=119;ECNT=1;FS=2.779;GERMQ=93;MBQ=32,32;MFRL=341,338;MMQ=60,60;MPOS=20;MQ=60;MQ0=0;MQRankSum=0;NALOD=1.82;NLOD=19.17;POPAF=6;ReadPosRankSum=-1.659;TLOD=33.37;LOF=(PLEKHG5|PLEKHG5|8|1.00);NMD=(PLEKHG5|PLEKHG5|8|1.00);ANN=ACTCACGTGCAAGCATCACACCGGCACGC|frameshift_variant&stop_gained|HIGH|PLEKHG5|PLEKHG5|transcript|NM_001265594.1|protein_coding|12/22|c.1191_1192insGCGTGCCGGTGTGATGCTTGCACGTGAG|p.Trp398fs|1428/4529|1191/2793|397/930||	GT:AD:AF:DP:F1R2:F2R1:SB	0/0:63,0:0.015:63:37,0:25,0:33,30,0,0	0/1:38,10:0.216:48:19,5:17,4:16,22,6,4
1	6471577	.	A	ACTCACGTGCAAGCATCACACCGGCACGC	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=49,52|6,4;ClippingRankSum=-1.498;DP=119;ECNT=1;FS=2.779;GERMQ=93;MBQ=32,32;MFRL=341,338;MMQ=60,60;MPOS=20;MQ=60;MQ0=0;MQRankSum=0;NALOD=1.82;NLOD=19.17;POPAF=6;ReadPosRankSum=-1.659;TLOD=33.37;LOF=(PLEKHG5|PLEKHG5|8|1.00);NMD=(PLEKHG5|PLEKHG5|8|1.00);ANN=ACTCACGTGCAAGCATCACACCGGCACGC|frameshift_variant&stop_gained|HIGH|PLEKHG5|PLEKHG5|transcript|NM_020631.4|protein_coding|12/21|c.1191_1192insGCGTGCCGGTGTGATGCTTGCACGTGAG|p.Trp398fs|1343/4644|1191/3021|397/1006||	GT:AD:AF:DP:F1R2:F2R1:SB	0/0:63,0:0.015:63:37,0:25,0:33,30,0,0	0/1:38,10:0.216:48:19,5:17,4:16,22,6,4
1	6471577	.	A	ACTCACGTGCAAGCATCACACCGGCACGC	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=49,52|6,4;ClippingRankSum=-1.498;DP=119;ECNT=1;FS=2.779;GERMQ=93;MBQ=32,32;MFRL=341,338;MMQ=60,60;MPOS=20;MQ=60;MQ0=0;MQRankSum=0;NALOD=1.82;NLOD=19.17;POPAF=6;ReadPosRankSum=-1.659;TLOD=33.37;LOF=(PLEKHG5|PLEKHG5|8|1.00);NMD=(PLEKHG5|PLEKHG5|8|1.00);ANN=ACTCACGTGCAAGCATCACACCGGCACGC|frameshift_variant&stop_gained|HIGH|PLEKHG5|PLEKHG5|transcript|NM_001042663.1|protein_coding|13/22|c.1359_1360insGCGTGCCGGTGTGATGCTTGCACGTGAG|p.Trp454fs|1460/4761|1359/3189|453/1062||	GT:AD:AF:DP:F1R2:F2R1:SB	0/0:63,0:0.015:63:37,0:25,0:33,30,0,0	0/1:38,10:0.216:48:19,5:17,4:16,22,6,4
1	6471577	.	A	ACTCACGTGCAAGCATCACACCGGCACGC	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=49,52|6,4;ClippingRankSum=-1.498;DP=119;ECNT=1;FS=2.779;GERMQ=93;MBQ=32,32;MFRL=341,338;MMQ=60,60;MPOS=20;MQ=60;MQ0=0;MQRankSum=0;NALOD=1.82;NLOD=19.17;POPAF=6;ReadPosRankSum=-1.659;TLOD=33.37;LOF=(PLEKHG5|PLEKHG5|8|1.00);NMD=(PLEKHG5|PLEKHG5|8|1.00);ANN=ACTCACGTGCAAGCATCACACCGGCACGC|frameshift_variant&stop_gained|HIGH|PLEKHG5|PLEKHG5|transcript|NM_198681.3|protein_coding|13/22|c.1422_1423insGCGTGCCGGTGTGATGCTTGCACGTGAG|p.Trp475fs|1972/5273|1422/3252|474/1083||	GT:AD:AF:DP:F1R2:F2R1:SB	0/0:63,0:0.015:63:37,0:25,0:33,30,0,0	0/1:38,10:0.216:48:19,5:17,4:16,22,6,4
```



## to determine where is missing position
```bash
SnpSift filter \
  -noLog \
  "( ANN[*].EFFECT has 'missense_variant' )"  \
  mutect2_normal_tumor_hg19-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf  | \
  /courses/software/vcfInfoOnePerLine.pl | \
  SnpSift extractFields \
  - \
  "CHROM" "POS" "ANN[*].GENE" "ANN[*].TRID" "EFF[*].HGVS_P" "ANN[*].HGVS_C" "ANN[*].EFFECT" | less

  ```
    mutect2_normal_tumor_hg19-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf  | less



## to determine 
```bash
SnpSift filter \
  -noLog \
  "( ANN[*].EFFECT has 'missense_variant' )"  \
  mutect2_normal_tumor_hg19-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf  | \
  /courses/software/vcfInfoOnePerLine.pl | \
  SnpSift extractFields \
  - \
  "CHROM" "POS" "ANN[*].GENE" "ANN[*].TRID" "EFF[*].HGVS_P" "ANN[*].HGVS_C" "ANN[*].EFFECT" | \
  grep 'missense_variant' | less

```
