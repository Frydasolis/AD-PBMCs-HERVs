#!/bin/bash
#SBATCH --job-name=STAR_AD_Herpes
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=/storage/lemus_g/roldan/AD/alignment_%j.out

# ------
source /home/roldan/anaconda3/etc/profile.d/conda.sh || source /home/roldan/miniconda3/etc/profile.d/conda.sh || source /etc/profile.d/conda.sh
conda activate star_env

# ------
FASTQ_DIR="/storage/lemus_g/roldan/AD/fastqs_finales"
INDEX="/home/roldan/AD/resources/index_human_herpes"
WHITELIST="/home/roldan/AD/resources/737K-august-2016.txt"
PLAN="/home/roldan/AD/resources/merge_plan.csv"
OUT_DIR="/storage/lemus_g/roldan/AD/results"

mkdir -p $OUT_DIR

# ------
DONORS=$(tail -n +2 "$PLAN" | cut -d',' -f1 | sort -u)

for GSM in $DONORS; do
    echo "Processing Donor: $GSM"
    SRRS=$(awk -F',' -v donor="$GSM" '$1==donor {print $3}' "$PLAN" | tr -d '\r')
    
    R4_FILES=""
    R3_FILES=""

    for SRR in $SRRS; do
        F4="$FASTQ_DIR/${SRR}_4.fastq"
        F3="$FASTQ_DIR/${SRR}_3.fastq"
        if [[ -f "$F4" && -f "$F3" ]]; then
            R4_FILES+="$F4,"
            R3_FILES+="$F3,"
        fi
    done

    R4_FILES=${R4_FILES%,}
    R3_FILES=${R3_FILES%,}

    if [[ -n "$R4_FILES" ]]; then
        
        /home/roldan/.conda/envs/star_env/bin/STAR --runThreadN 16 \
             --genomeDir "$INDEX" \
             --readFilesIn "$R4_FILES" "$R3_FILES" \
             --soloType CB_UMI_Simple \
             --soloCBstart 1 --soloCBlen 16 \
             --soloUMIstart 17 --soloUMIlen 10 \
             --soloBarcodeReadLength 0 \
             --soloCBwhitelist "$WHITELIST" \
             --soloStrand Reverse \
             --soloFeatures Gene GeneFull \
             --soloMultiMappers EM \
             --outFilterMultimapNmax 500 \
             --outFilterMultimapScoreRange 5 \
             --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
             --soloUMIfiltering MultiGeneUMI_CR \
             --soloUMIdedup 1MM_CR \
             --clipAdapterType CellRanger4 \
             --outFilterScoreMin 30 \
             --limitOutSJcollapsed 5000000 \
             --outSAMattributes NH HI AS NM nM MD CR CY UR UY CB UB GX GN sS sQ sM \
             --outSAMunmapped Within \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix "${OUT_DIR}/${GSM}_"
    fi
done
