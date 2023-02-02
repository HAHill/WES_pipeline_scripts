

#!/usr/bin

echo "STEP 4: Base quality recalibration"

ref="/rsrch4/scratch/bcb/hahill1/files/support_files/Homo_sapiens_assembly38.fasta"
data="/rsrch4/scratch/bcb/hahill1/files/yliu/data"
ref="/rsrch4/scratch/bcb/hahill1/files/support_files/Homo_sapiens_assembly38.dbsnp138.vcf"

#Merged, marked reads
deduped="/rsrch4/scratch/bcb/hahill1/files/yliu/aligned/deduped"
recalibrated="/rsrch4/scratch/bcb/hahill1/files/yliu/aligned/recalibrated/"



gatk ApplyBQSR -I ${deduped}/1-MCH1_S1_merged_d.sam  -R ${ref} --bqsr-recal-file ${data}/S1_recal_data.table -O ${recalibrated}1-MCH1_S1_merged_d_r.sam 
gatk ApplyBQSR -I ${deduped}/1-MCH2_S2_merged_d.sam  -R ${ref} --bqsr-recal-file ${data}/S2_recal_data.table -O ${recalibrated}1-MCH2_S2_merged_d_r.sam 
gatk ApplyBQSR -I ${deduped}/1-MCH3_S3_merged_d.sam  -R ${ref} --bqsr-recal-file ${data}/S3_recal_data.table -O ${recalibrated}1-MCH3_S3_merged_d_r.sam
gatk ApplyBQSR -I ${deduped}/1-MCH4_S4_merged_d.sam  -R ${ref} --bqsr-recal-file ${data}/S4_recal_data.table -O ${recalibrated}1-MCH4_S4_merged_d_r.sam
gatk ApplyBQSR -I ${deduped}/1-MCH5_S5_merged_d.sam  -R ${ref} --bqsr-recal-file ${data}/S5_recal_data.table -O ${recalibrated}1-MCH5_S5_merged_d_r.sam
gatk ApplyBQSR -I ${deduped}/1-MCH1_S6G_merged_d.sam  -R ${ref} --bqsr-recal-file ${data}/S6_recal_data.table -O ${recalibrated}1-MCH1_S6G_merged_d_r.sam
gatk ApplyBQSR -I ${deduped}/1-MCH2_S7G_merged_d.sam  -R ${ref} --bqsr-recal-file ${data}/S7_recal_data.table -O ${recalibrated}1-MCH2_S7G_merged_d_r.sam
gatk ApplyBQSR -I ${deduped}/1-MCH1_S8G_merged_d.sam  -R ${ref} --bqsr-recal-file ${data}/S8_recal_data.table -O ${recalibrated}1-MCH3_S8G_merged_d_r.sam
gatk ApplyBQSR -I ${deduped}/1-MCH1_S9G_merged_d.sam  -R ${ref} --bqsr-recal-file ${data}/S9_recal_data.table -O ${recalibrated}1-MCH4_S9G_merged_d_r.sam
gatk ApplyBQSR -I ${deduped}/1-MCH1_S10G_merged_d.sam  -R ${ref} --bqsr-recal-file ${data}/S10_recal_data.table -O ${recalibrated}1-MCH5_S10G_merged_d_r.sam



