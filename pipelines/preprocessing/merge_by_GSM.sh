#!/bin/bash

FASTQ_DIR="/storage/lemus_g/roldan/AD/fastq"
MAP="/home/roldan/AD/resources/SRR_GSM_map.csv"
OUT="/storage/lemus_g/roldan/AD/merged_fastq"

mkdir -p $OUT

echo "Starting merge..."

tail -n +2 "$MAP" | tr -d '\r' | while IFS=',' read SM Assay Disease SRR
do

    R1="$FASTQ_DIR/${SRR}_3.fastq"
    R2="$FASTQ_DIR/${SRR}_4.fastq"

    if [[ -s "$R1" && -s "$R2" ]]; then

        echo "Merging $SRR -> $SM"

        cat "$R1" >> "$OUT/${SM}_R1.fastq"
        cat "$R2" >> "$OUT/${SM}_R2.fastq"

    else

        echo "WARNING: missing $SRR"

    fi

done

echo "Merge completed."
