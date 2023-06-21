# WES_pipeline_scripts
Basic scripts for conducting WES analysis on germline and somatic DNA samples 
** "wes_script.sh" puts all pieces of the pipeline together, but is usually run in pieces for QC assessment.

STEPS:
1. Run fastp for QC and moderate trimming on demultiplexed data (FASTQs). Check QC metrics for all files. Can process UMI if needed.
2. Map (align) to reference genome using BWA-MEM algorithm.
3. Merge SAM or BAM Files (if dual end reads) with Picard.
4. Mark duplicates with GATK Spark and sort– resulting in sorted, deduped BAM files.
5. Do base quality recalibration with GATK (first, build a model with data tables – second, apply the model to adjust the base quality scores.
6. Checking QC step only: collect alignment and insert size metrics with GATK. Can generate a histogram of insert sizes.
7. Call Variants – use GATK haplotype caller on sorted, deduped, base-quality recalibrated BAM files. Output = VCF
8. Call Somatic Variants – use Mutect2. Use germline if available – if not utilize panel of normal or population databases. Output = Somatic VCF
9. Use Funcotator (or preferred annotator) to annotate and filter somatic variants. Output = VCF or MAF. 
10. Optional, call CNV utilize ExomeCNV or preferred software
11. Visualization and analyses.

Special thanks to https://github.com/kpatel427 for her amazing YouTube bioinformatics tutorials and super useful script.
