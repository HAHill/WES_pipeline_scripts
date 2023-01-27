#!/usr/bin

################################################### 
Prep files (TO BE 
GENERATED ONLY ONCE) 
##########################################################

# download reference files
wget -P ~/Desktop/demo/supporting_files/hg38/ 
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip ~/Desktop/demo/supporting_files/hg38/hg38.fa.gz

# index ref - .fai file before running haplotype caller
samtools faidx ~/Desktop/demo/supporting_files/hg38/hg38.fa

# ref dict - .dict file before running haplotype caller
gatk CreateSequenceDictionary 
R=~/Desktop/demo/supporting_files/hg38/hg38.fa 
O=~/Desktop/demo/supporting_files/hg38/hg38.dict



###################################################### VARIANT CALLING STEPS ####################################################################


# directories
ref="/rsrch4/scratch/bcb/hahill1/files/support_files/Homo_sapiens_assembly38.fasta.fai"
known_sites="/rsrch4/scratch/bcb/hahill1/files/support_files/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/rsrch4/scratch/bcb/hahill1/files/yliu/aligned"
reads="/rsrch4/scratch/bcb/hahill1/files/yliu/YL-MCL-WES_HC320"
results= "/rsrch4/scratch/bcb/hahill1/files/results"
data= "/rsrch4/scratch/bcb/hahill1/files/yliu/data"

# ------------------
# Alt STEP 1: QC and preprocessing - Run fastp # ---------------------------

reads= "/rsrch4/scratch/bcb/hahill1/files/yliu/YL-MCL-WES_HC320"

module load fastp

echo "STEP 1: QC and preprocessing - Run fastp" 

cat ids.txt | parallel echo fastp -I ${reads}/{}.fq.gz -O ${reads}/
cat ids.txt | parallel echo fastp -i ${reads}/{}.fq.gz -o ${reads}/

#Fastp automatically filters bad reads and trims if adapters are detected. Will still need to do base quality recalibration.

# -------------------
# STEP 1: QC - Run fastqc # -------------------

echo "STEP 1: QC - Run fastqc"

fastqc ${reads}/1-MCH1_S*.fq.gz -o ${reads}/
fastqc ${reads}/1-MCH1_S1_L001_R2_001.fq.gz -o ${reads}/

#If quality looks ok, no trimming required.


# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

#!/usr/bin

echo "STEP 2: Map to reference using BWA-MEM"

# BWA index reference and paths to storage folders
ref="/rsrch4/scratch/bcb/hahill1/files/support_files/Homo_sapiens_assembly38.fasta"
aligned_reads="/rsrch4/scratch/bcb/hahill1/files/yliu/aligned"
reads="/rsrch4/scratch/bcb/hahill1/files/yliu/results"


bwa index ${ref}

# BWA alignment

bwa mem -t 4 -R "@RG\tID:1-MCH1_S1_L001_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S1R1L1.fq.gz ${reads}/S1R3L1.fq.gz > ${aligned_reads}/1-MCH1_S1_L001.sam
bwa mem -t 4 -R "@RG\tID:1-MCH1_S1_L002_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S1R1L2.fq.gz ${reads}/S1R3L2.fq.gz > ${aligned_reads}/1-MCH1_S1_L002.sam
bwa mem -t 4 -R "@RG\tID:1-MCH2_S2_L001_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S2R1L1.fq.gz ${reads}/S2R3L1.fq.gz > ${aligned_reads}/1-MCH2_S2_L001.sam
bwa mem -t 4 -R "@RG\tID:1-MCH2_S2_L002_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S2R1L2.fq.gz ${reads}/S2R3L2.fq.gz > ${aligned_reads}/1-MCH2_S2_L002.sam
bwa mem -t 4 -R "@RG\tID:1-MCH3_S3_L001_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S3R1L1.fq.gz ${reads}/S3R3L1.fq.gz > ${aligned_reads}/1-MCH3_S3_L001.sam
bwa mem -t 4 -R "@RG\tID:1-MCH3_S3_L002_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S3R1L2.fq.gz ${reads}/S3R3L2.fq.gz > ${aligned_reads}/1-MCH3_S3_L002.sam
bwa mem -t 4 -R "@RG\tID:1-MCH4_S4_L001_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S4R1L1.fq.gz ${reads}/S4R3L1.fq.gz > ${aligned_reads}/1-MCH4_S4_L001.sam
bwa mem -t 4 -R "@RG\tID:1-MCH4_S4_L002_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S4R1L2.fq.gz ${reads}/S4R3L2.fq.gz > ${aligned_reads}/1-MCH4_S4_L002.sam
bwa mem -t 4 -R "@RG\tID:1-MCH5_S5_L001_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S5R1L1.fq.gz ${reads}/S5R3L1.fq.gz > ${aligned_reads}/1-MCH5_S5_L001.sam
bwa mem -t 4 -R "@RG\tID:1-MCH5_S5_L002_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S5R1L2.fq.gz ${reads}/S5R3L2.fq.gz > ${aligned_reads}/1-MCH5_S5_L002.sam
bwa mem -t 4 -R "@RG\tID:1-MCH1_S6G_L001_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S6GR1L1.fq.gz ${reads}/S6GR3L1.fq.gz > ${aligned_reads}/1-MCH1_S6G_L001.sam
bwa mem -t 4 -R "@RG\tID:1-MCH1_S6G_L002_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S6GR1L2.fq.gz ${reads}/S6GR3L2.fq.gz > ${aligned_reads}/1-MCH1_S6G_L002.sam
bwa mem -t 4 -R "@RG\tID:1-MCH2_S7G_L001_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S7GR1L1.fq.gz ${reads}/S7GR3L1.fq.gz > ${aligned_reads}/1-MCH2_S7G_L001.sam
bwa mem -t 4 -R "@RG\tID:1-MCH2_S7G_L002_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S7GR1L2.fq.gz ${reads}/S7GR3L2.fq.gz > ${aligned_reads}/1-MCH2_S7G_L002.sam
bwa mem -t 4 -R "@RG\tID:1-MCH3_S8G_L001_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S8GR1L1.fq.gz ${reads}/S8GR3L1.fq.gz > ${aligned_reads}/1-MCH3_S8G_L001.sam
bwa mem -t 4 -R "@RG\tID:1-MCH3_S8G_L002_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S8GR1L2.fq.gz ${reads}/S8GR3L2.fq.gz > ${aligned_reads}/1-MCH3_S8G_L002.sam
bwa mem -t 4 -R "@RG\tID:1-MCH4_S9G_L001_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S9GR1L1.fq.gz ${reads}/S9GR3L1.fq.gz > ${aligned_reads}/1-MCH4_S9G_L001.sam
bwa mem -t 4 -R "@RG\tID:1-MCH4_S9G_L002_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S9GR1L2.fq.gz ${reads}/S9GR3L2.fq.gz > ${aligned_reads}/1-MCH4_S9G_L002.sam
bwa mem -t 4 -R "@RG\tID:1-MCH5_S10G_L001_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S10GR1L1.fq.gz ${reads}/S10GR3L1.fq.gz > ${aligned_reads}/1-MCH5_S10G_L001.sam
bwa mem -t 4 -R "@RG\tID:1-MCH5_S10G_L002_001\tPL:ILLUMINA\tSM:HC320" ${ref} ${reads}/S10GR1L2.fq.gz ${reads}/S10GR3L2.fq.gz > ${aligned_reads}/1-MCH5_S10G_L002.sam


# -----------------------------------------
# STEP 3:Merge Sam Files,  Mark Duplicates and Sort - GATK4
# -----------------------------------------

module load picard

aligned_reads="/rsrch4/scratch/bcb/hahill1/files/yliu/aligned"
merged="/rsrch4/scratch/bcb/hahill1/files/yliu/aligned/merged"
PICARDHOME="/risapps/noarch/picard/2.23.8/build/libs"

java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH1_S1_L001.sam  -I ${aligned_reads}/1-MCH1_S1_L002.sam -O ${merged}1-MCH1_S1_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH2_S2_L001.sam  -I ${aligned_reads}/1-MCH2_S2_L002.sam  -O ${merged}1-MCH2_S2_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH3_S3_L001.sam  -I ${aligned_reads}/1-MCH3_S3_L002.sam  -O ${merged}1-MCH3_S3_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH4_S4_L001.sam  -I ${aligned_reads}/1-MCH4_S4_L002.sam  -O ${merged}1-MCH4_S4_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH5_S5_L001.sam  -I ${aligned_reads}/1-MCH5_S5_L002.sam  -O ${merged}1-MCH5_S5_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH1_S6G_L001.sam  -I ${aligned_reads}/1-MCH1_S6G_L002.sam  -O ${merged}1-MCH1_S6G_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH2_S7G_L001.sam  -I ${aligned_reads}/1-MCH2_S7G_L002.sam  -O ${merged}1-MCH2_S7G_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH3_S8G_L001.sam  -I ${aligned_reads}/1-MCH3_S8G_L002.sam  -O ${merged}1-MCH3_S8G_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH4_S9G_L001.sam  -I ${aligned_reads}/1-MCH4_S9G_L002.sam  -O  ${merged}1-MCH4_S9G_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH5_S10G_L001.sam  -I ${aligned_reads}/1-MCH5_S10G_L002.sam  -O ${merged}1-MCH5_S10G_merged.sam

echo "STEP 3: Mark Duplicates and Sort - GATK4"

aligned_reads="/rsrch4/scratch/bcb/hahill1/files/yliu/aligned"
merged="/rsrch4/scratch/bcb/hahill1/files/yliu/aligned/merged"

java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${merged}/1-MCH1_S1_L001.merged.sam  -I ${aligned_reads}/1-MCH1_S1_L002.sam -O ${merged}1-MCH1_S1_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH2_S2_L001.sam  -I ${aligned_reads}/1-MCH2_S2_L002.sam  -O ${merged}1-MCH2_S2_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH3_S3_L001.sam  -I ${aligned_reads}/1-MCH3_S3_L002.sam  -O ${merged}1-MCH3_S3_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH4_S4_L001.sam  -I ${aligned_reads}/1-MCH4_S4_L002.sam  -O ${merged}1-MCH4_S4_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH5_S5_L001.sam  -I ${aligned_reads}/1-MCH5_S5_L002.sam  -O ${merged}1-MCH5_S5_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH1_S6G_L001.sam  -I ${aligned_reads}/1-MCH1_S6G_L002.sam  -O ${merged}1-MCH1_S6G_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH2_S7G_L001.sam  -I ${aligned_reads}/1-MCH2_S7G_L002.sam  -O ${merged}1-MCH2_S7G_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH3_S8G_L001.sam  -I ${aligned_reads}/1-MCH3_S8G_L002.sam  -O ${merged}1-MCH3_S8G_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH4_S9G_L001.sam  -I ${aligned_reads}/1-MCH4_S9G_L002.sam  -O  ${merged}1-MCH4_S9G_merged.sam
java -jar $PICARDHOME/picard.jar MergeSamFiles  -I ${aligned_reads}/1-MCH5_S10G_L001.sam  -I ${aligned_reads}/1-MCH5_S10G_L002.sam  -O ${merged}1-MCH5_S10G_merged.sam

gatk MarkDuplicatesSpark -I ${aligned_reads}/1-MCH1_S1_L001.merged.sam -O 
${aligned_reads}/1-MCH1_S1_L001_sorted_dedup_reads.bam

# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------

echo "STEP 4: Base quality recalibration"

# 1. build the model
gatk BaseRecalibrator -I ${aligned_reads}/1-MCH1_S1_L001_sorted_dedup_reads.bam 
-R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table


# 2. Apply the model to adjust the base quality scores
gatk ApplyBQSR -I ${aligned_reads}/1-MCH1_S1_L001_sorted_dedup_reads.bam -R 
${ref} --bqsr-recal-file {$data}/recal_data.table -O 
${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 


# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------

echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics 
R=${ref} 
-I ${aligned_reads}/1-MCH1_S1_L001_sorted_dedup_bqsr_reads.bam 
O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics 
INPUT=${aligned_reads}/1-MCH1_S1_L001_sorted_dedup_bqsr_reads.bam 
OUTPUT=${aligned_reads}/insert_size_metrics.txt 
HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf


# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ----------------------------------------------

echo "STEP 6: Call Variants - gatk haplotype caller"

gatk HaplotypeCaller
 -R ${ref}
 -I ${aligned_reads}/1-MCH1_S1_L001_sorted_dedup_bqsr_reads.bam
 -O ${results}/raw_variants.vcf

# ----------------------------------------------
# Alt STEP 6: Call Variants - Mutect2
# ----------------------------------------------

echo "STEP 7: Call Variants - gatk Mutect2"

gatk Mutect2 
  -R ${ref} 
  -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 
  -I ${aligned_reads}/1-MCH1_S1_L00G1_sorted_dedup_bqsr_reads.bam
  -tumor 1-MCH1_S1_tumor
  -normal 1-MCH1_S1_normal
  --germline-resource af-only-gnomad.vcf.gz 
  --panel-of-normals pon.vcf.gz 
  -O ${results}/somatic.vcf.gz




# extract SNPs & INDELS

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type 
SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type 
INDEL -O ${results}/raw_indels.vcf





