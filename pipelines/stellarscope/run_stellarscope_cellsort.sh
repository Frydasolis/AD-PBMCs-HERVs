#!/bin/bash
#SBATCH --job-name=stellar_sort
#SBATCH --output=stellar_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=24:00:00
# Nota: He eliminado --partition para evitar el error. 
# Slurm usará la partición por defecto del servidor borromeo.

# 1. Conda
source /software/anaconda/etc/profile.d/conda.sh
conda activate stellarscope_protocol_env

# 
INPUT_DIR="/storage/lemus_g/roldan/AD/results"
OUTDIR="/storage/lemus_g/roldan/AD/Stellarscope/cellsort"
TMPDIR="/tmp"

mkdir -p "$OUTDIR"

# 3. Proccessing
echo "Stellarscope working..."

for INPUT_BAM in "$INPUT_DIR"/*.bam; do
    FILENAME=$(basename "$INPUT_BAM" .bam)
    SAMPLE_ID=$(echo "$FILENAME" | cut -d'_' -f1)
    BARCODES="$INPUT_DIR/${SAMPLE_ID}_Solo.out/Gene/filtered/barcodes.tsv"

    if [ -f "$BARCODES" ]; then
        echo "------------------------------------------------"
        echo "Proccessing: $SAMPLE_ID"
        stellarscope cellsort \
          --nproc $SLURM_CPUS_PER_TASK \
          --tempdir "$TMPDIR" \
          --outfile "$OUTDIR/${SAMPLE_ID}.Aligned.sortedByCB.bam" \
          "$INPUT_BAM" \
          "$BARCODES"
    else
        echo "SALTANDO $SAMPLE_ID: There is not bardcodes  $BARCODES"
    fi
done
