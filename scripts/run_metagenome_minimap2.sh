#!/bin/bash

#########################################################################
# Author: Alexander Graf  (alexander.graf@lmu.de)
#
# Run Nextflow metagenome workflow using Minimap2 on an input fastq folder.
#
##########################################################################

INPUT_DIR=$1
OUTPUT_DIR=$2
PROFILE=$3         # must be "default" or "full"
THREADS=${4:-8}     # optional, default: 8

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$PROFILE" ]; then
    echo ""
    echo "Run alignment-based metagenomics pipeline using Minimap2"
    echo ""
    echo "Usage:"
    echo "  $(basename "$0") <FASTQ_INPUT_FOLDER> <OUTPUT_DIR> <PROFILE> [THREADS]"
    echo ""
    echo "Arguments:"
    echo "  FASTQ_INPUT_FOLDER   Path to input folder with barcode FASTQs"
    echo "  OUTPUT_DIR           Path to output directory"
    echo "  PROFILE              'default' (3 genomes) or 'full' (7 genomes)"
    echo "  THREADS              (Optional) Number of CPU threads [default: 8]"
    echo ""
    echo "Examples:"
    echo "  $(basename "$0") /data/fastqs results default"
    echo "  $(basename "$0") /data/fastqs results full 16"
    echo ""
    exit 1
fi

if [[ "$PROFILE" != "default" && "$PROFILE" != "full" ]]; then
  echo "Invalid profile '$PROFILE'. Must be 'default' or 'full'."
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TAXONOMY_PATH="$SCRIPT_DIR/../referenceGenome/kraken2_db/taxonomy/new_taxdump.tar.gz"

eval "$(conda shell.bash hook)"
conda activate zugspitze_metagenome || {
  echo "❌ Failed to activate conda environment 'zugspitze_metagenome'"
  exit 6
}

# Determine reference set
if [[ "$PROFILE" == "full" ]]; then
  REF_DIR="$SCRIPT_DIR/../referenceGenome/all_references"
elif [[ "$PROFILE" == "default" ]]; then
  REF_DIR="$SCRIPT_DIR/../referenceGenome/combined_genomes"
else
  echo "Unknown profile '$PROFILE'. Use 'default' or 'full'."
  exit 2
fi

# Check that required files exist
if [[ ! -f "$REF_DIR/$(basename "$REF_DIR").mmi" ]]; then
  echo "Reference index not found at: $REF_DIR/$(basename "$REF_DIR").mmi"
  exit 3
fi

if [[ ! -f "$REF_DIR/ref2taxid.targloci.tsv" ]]; then
  echo "ref2taxid file not found at: $REF_DIR/ref2taxid.targloci.tsv"
  exit 4
fi

if [[ ! -f "$TAXONOMY_PATH" ]]; then
  echo "Taxonomy file not found at: $TAXONOMY_PATH"
  exit 5
fi

echo "✅ Running pipeline with profile: $PROFILE"
echo "📁 Reference directory: $REF_DIR"
echo "🚀 Threads: $THREADS"

nextflow run epi2me-labs/wf-metagenomics \
        --fastq $INPUT_DIR \
        --classifier minimap2 \
        --reference "$REF_DIR/$(basename "$REF_DIR").mmi" \
        --ref2taxid "$REF_DIR/ref2taxid.targloci.tsv" \
        --out_dir $OUTPUT_DIR \
        --keep_bam \
        --threads $THREADS \
        --taxonomy "$TAXONOMY_PATH" \
        -process.maxForks=1

