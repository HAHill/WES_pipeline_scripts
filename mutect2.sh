#!/usr/bin

ref="/rsrch4/scratch/bcb/hahill1/files/support_files/Homo_sapiens_assembly38.fasta"
gnomad="/rsrch4/files/support_files/af-only-gnomad.hg38.vcf"
pon="rsrch4/files/support_files/somatic-hg38-1000g_pon.hg38.vcf.gz.vcf"
aligned_reads="/rsrch4/scratch/bcb/hahill1/files/yliu/aligned"
results="/rsrch4/scratch/bcb/hahill1/files/yliu/aligned/mutect2"

echo "STEP 7: Call Variants - gatk Mutect2"

gatk Mutect2 
  -R ${ref} -I ${aligned_reads}/recalibrated1-MCH1_S1_merged_d_r.sam -I ${aligned_reads}/1-MCH1_S6G_merged_d_r.sam -normal 1-MCH1_S6G --germline-resource {gnomad} --panel-of-normals {pon} -O ${results}/MCH1_somatic.vcf.gz
  -R ${ref} -I ${aligned_reads}/recalibrated1-MCH2_S2_merged_d_r.sam -I ${aligned_reads}/1-MCH2_S7G_merged_d_r.sam-normal 1-MCH2_S7G --germline-resource af-only-gnomad.vcf.gz --panel-of-normals pon.vcf.gz -O ${results}/MCH2_somatic.vcf.gz
  -R ${ref} -I ${aligned_reads}/recalibrated1-MCH3_S3_merged_d_r.sam -I ${aligned_reads}/1-MCH3_S8G_merged_d_r.sam -normal 1-MCH3_S8G --germline-resource af-only-gnomad.vcf.gz --panel-of-normals pon.vcf.gz -O ${results}/MCH3_somatic.vcf.gz
  -R ${ref} -I ${aligned_reads}/recalibrated1-MCH4_S4_merged_d_r.sam -I ${aligned_reads}/1-MCH4_S9G_merged_d_r.sam -normal 1-MCH4_S9G --germline-resource af-only-gnomad.vcf.gz --panel-of-normals pon.vcf.gz -O ${results}/MCH4_somatic.vcf.gz
  -R ${ref} -I ${aligned_reads}/recalibrated1-MCH5_S5_merged_d_r.sam -I ${aligned_reads}/1-MCH5_S10G_merged_d_r.sam -normal 1-MCH5_S10G --germline-resource af-only-gnomad.vcf.gz --panel-of-normals pon.vcf.gz -O ${results}/MCH5_somatic.vcf.gz

