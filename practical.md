
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
        --threads 10 \
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
cd ./reports/multiqc_step01/
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
    /courses/master_b2s/raw_data/s02_vc/cancer/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz \
    /courses/master_b2s/raw_data/s02_vc/cancer/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz \
    /courses/master_b2s/raw_data/s02_vc/cancer_trim/SLGFSK-N_paired_R1.fastq.gz /courses/master_b2s/raw_data/s02_vc/cancer_trim/SLGFSK-N_unpaired_R1.fastq.gz \
    /courses/master_b2s/raw_data/s02_vc/cancer_trim/SLGFSK-N_paired_R2.fastq.gz /courses/master_b2s/raw_data/s02_vc/cancer_trim/SLGFSK-N_unpaired_R2.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads \
    HEADCROP:3 \
    TRAILING:10 \
    MINLEN:25



# Run Trimmomatic on the Tumor (T) sample
java -jar /courses/master_b2s/software/Trimmomatic/trimmomatic-0.40.jar PE -threads 4 \
    /courses/master_b2s/raw_data/s02_vc/cancer/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz \
    /courses/master_b2s/raw_data/s02_vc/cancer/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz \
    /courses/master_b2s/raw_data/s02_vc/cancer_trim/SLGFSK-T_paired_R1.fastq.gz /courses/master_b2s/raw_data/s02_vc/cancer_trim/SLGFSK-T_unpaired_R1.fastq.gz \
    /courses/master_b2s/raw_data/s02_vc/cancer_trim/SLGFSK-T_paired_R2.fastq.gz /courses/master_b2s/raw_data/s02_vc/cancer_trim/SLGFSK-T_unpaired_R2.fastq.gz \
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
cd ./reports/multiqc_step01/
# run analyse
multiqc ./reports/fastqc/ -o ./reports/multiqc_step01/ -t 10

# see file
ls -l ./reports/multiqc_step01/
```




```bash
mkdir -p fastqc/original fastqc/trimmed trimm_out multiqc

for R1 in *_R1_001.pe.fq.gz; do

    R2=${R1/_R1_/_R2_}

    BASE=$(basename $R1 _R1_001.pe.fq.gz)

    echo "Processing $BASE..."

    java -jar /courses/master_b2s/software/Trimmomatic/trimmomatic-0.40.jar PE -threads 4 $R1 $R2 trimm_out/${BASE}_R1_paired.fq.gz trimm_out/${BASE}_R1_unpaired.fq.gz trimm_out/${BASE}_R2_paired.fq.gz trimm_out/${BASE}_R2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> trimm_out/${BASE}.log

done



fastqc *_R1_001.pe.fq.gz *_R2_001.pe.fq.gz -o fastqc/original/

fastqc trimm_out/*_paired.fq.gz -o fastqc/trimmed/



multiqc . -o multiqc/

```




--- 

6. index reference genome

```bash
mkdir /courses/master_b2s/references/GRCh/hg19#
bwa index hg19.chr5_12_17.fa.gz

```

7. alignemnt bam

```bash

REFERENCE_SEQUENCE=/courses/master_b2s/references/GRCh/hg19/hg19.chr5_12_17.fa.gz
LEFT_READS=/courses/master_b2s/raw_data/s02_vc/cancer/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
RIGHT_READS=/courses/master_b2s/raw_data/s02_vc/cancer/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
#SAMPLE=`basename $LEFT_READS _1.fq.gz`
SAMPLE=syn3_normal
SAM_FILE=./session02vc/alignments/${SAMPLE}_GRCh38.p7.sam

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
LEFT_READS=/courses/master_b2s/raw_data/s02_vc/cancer/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
RIGHT_READS=/courses/master_b2s/raw_data/s02_vc/cancer/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz
SAMPLE=syn3_tumor
SAM_FILE=./session02vc/alignments/${SAMPLE}_GRCh38.p7.sam

# Align reads with bwa
bwa mem \
    -M \
    -t 10 \
    -R "@RG\tID:$SAMPLE\tPL:illumina\tPU:$SAMPLE\tSM:$SAMPLE" \
    $REFERENCE_SEQUENCE \
    $LEFT_READS \
    $RIGHT_READS \
    -o $SAM_FILE
```

8. Alignment file processing


```bash

# Assign file paths to variables
SAMPLE_NAME=syn3_normal
SAM_FILE=~/session02vc/alignments/${SAMPLE_NAME}_GRCh38.p7.sam
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
SAMPLE_NAME=syn3_tumor
SAM_FILE=~/session02vc/alignments/${SAMPLE_NAME}_GRCh38.p7.sam
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
 cp /courses/master_b2s/session02vc/alignments/syn3_*_GRCh38.p7.coordinate_sorted.ba? session2vc/alignments/

```


```bash
# Assign variables
INPUT_BAM=~/session02vc/alignments/syn3_normal_GRCh38.p7.coordinate_sorted.bam
REFERENCE=/courses/master_b2s/references/GRCh/hg19/hg19.chr5_12_17.fa.gz
OUTPUT_METRICS_FILE=~/session02vc/reports/picard/syn3_normal/syn3_normal_GRCh38.p7.CollectAlignmentSummaryMetrics.txt

# Run Picard CollectAlignmentSummaryMetrics
#java -jar $PICARD/picard.jar 
picard CollectAlignmentSummaryMetrics \
  --INPUT $INPUT_BAM \
  --REFERENCE_SEQUENCE $REFERENCE \
  --OUTPUT $OUTPUT_METRICS_FILE


# Assign variables

INPUT_BAM=~/session02vc/alignments/syn3_tumor_GRCh38.p7.coordinate_sorted.bam
REFERENCE=/courses/master_b2s/references/GRCh/hg19/hg19.chr5_12_17.fa.gz
OUTPUT_METRICS_FILE=~/session02vc/reports/picard/syn3_tumor/syn3_tumor_GRCh38.p7.CollectAlignmentSummaryMetrics.txt

# Run Picard CollectAlignmentSummaryMetrics
#java -jar $PICARD/picard.jar 
picard CollectAlignmentSummaryMetrics \
  --INPUT $INPUT_BAM \
  --REFERENCE_SEQUENCE $REFERENCE \
  --OUTPUT $OUTPUT_METRICS_FILE

```

11. Running MultiQC Aggregating QC


```bash
less ~/session02vc/reports/picard/syn3_normal/syn3_normal_GRCh38.p7.CollectAlignmentSummaryMetrics.txt
```


```bash
tree ~/session02vc/reports/
PLACEHOLDER FOR COMMAND

```


```bash

# Assign variables
REPORTS_DIRECTORY=~/session02vc/reports/
NORMAL_SAMPLE_NAME=syn3_normal
TUMOR_SAMPLE_NAME=syn3_tumor
REFERENCE=GRCh38.p7
NORMAL_PICARD_METRICS=${REPORTS_DIRECTORY}picard/${NORMAL_SAMPLE_NAME}/${NORMAL_SAMPLE_NAME}_${REFERENCE}.CollectAlignmentSummaryMetrics.txt
TUMOR_PICARD_METRICS=${REPORTS_DIRECTORY}picard/${TUMOR_SAMPLE_NAME}/${TUMOR_SAMPLE_NAME}_${REFERENCE}.CollectAlignmentSummaryMetrics.txt
NORMAL_FASTQC_1=${REPORTS_DIRECTORY}fastqc/${NORMAL_SAMPLE_NAME}_1_fastqc.zip
NORMAL_FASTQC_2=${REPORTS_DIRECTORY}fastqc/${NORMAL_SAMPLE_NAME}_2_fastqc.zip
TUMOR_FASTQC_1=${REPORTS_DIRECTORY}fastqc/${TUMOR_SAMPLE_NAME}_1_fastqc.zip
TUMOR_FASTQC_2=${REPORTS_DIRECTORY}fastqc/${TUMOR_SAMPLE_NAME}_2_fastqc.zip
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

to solve issue here you need to change name in zip to correct name **syn3_normal_1_fastqc**

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
 --REFERENCE /n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
 --OUTPUT /n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.dict


```

```bash
# Assign variables
REFERENCE_SEQUENCE=/courses/master_b2s/references/GRCh/hg19/hg19.chr5_12_17.fa
REFERENCE_DICTIONARY=`echo ${REFERENCE_SEQUENCE%fa}dict`
NORMAL_SAMPLE_NAME=syn3_normal
NORMAL_BAM_FILE=~/session02vc/alignments/${NORMAL_SAMPLE_NAME}_GRCh38.p7.coordinate_sorted.bam
TUMOR_SAMPLE_NAME=syn3_tumor
TUMOR_BAM_FILE=~/session02vc/alignments/${TUMOR_SAMPLE_NAME}_GRCh38.p7.coordinate_sorted.bam
VCF_OUTPUT_FILE=~/session02vc/vcf_files/mutect2_${NORMAL_SAMPLE_NAME}_${TUMOR_SAMPLE_NAME}_GRCh38.p7-raw.vcf

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

NORMAL_SAMPLE_NAME=syn3_normal
NORMAL_BAM_FILE=~/session02vc/alignments/${NORMAL_SAMPLE_NAME}_GRCh38.p7.coordinate_sorted.bam

TUMOR_SAMPLE_NAME=syn3_tumor
TUMOR_BAM_FILE=~/session02vc/alignments/${TUMOR_SAMPLE_NAME}_GRCh38.p7.coordinate_sorted.bam

VCF_OUTPUT_FILE=~/session02vc/vcf_files/mutect2_somatic_turbo.vcf1

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
RAW_VCF_FILE=~/session02vc/vcf_files/mutect2_syn3_normal_syn3_tumor_GRCh38.p7-raw.vcf
LCR_FILE=/courses/master_b2s/references/GRCh/hg19/LCR-hs38_with_chr.bed
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

less session02vc/vcf_files/mutect2_syn3_normal_syn3_tumor_GRCh38.p7-pass-filt-LCR.vcf

grep -v "^#" session02vc/vcf_files/mutect2_syn3_normal_syn3_tumor_GRCh38.p7-pass-filt-LCR.vcf | wc -l

```


