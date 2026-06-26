# AD-PBMCs-HERVs
=======
# AD APOE Multi-Omics Pipeline

This repository contains preprocessing and analysis scripts for RNA-seq and scATAC-seq data used in the AD APOE project.

## Preprocessing

### merge_by_GSM.sh
This script merges SRR-level FASTQ files into GSM-level samples using a mapping file.

### Input
- SRR_GSM_map.csv
- Raw FASTQ files (SRR level)

### Output
- Merged FASTQ files per GSM sample:
  - GSMxxx_R1.fastq
  - GSMxxx_R2.fastq
 
