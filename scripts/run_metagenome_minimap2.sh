#!/bin/bash

#########################################################################
# Author: Alexander Graf  (alexander.graf@lmu.de)
#
# Run Nextflow metagenome workflow using Minimap2 on an input fastq folder.
#
##########################################################################

INPUT_DIR=$1
OUTPUT_DIR=$2
THREADS=${3:-8}

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] ; then
    echo ""
    echo "Run alignment-based metagenomics pipeline"
    echo ""
    echo "Usage:"
    echo "  $(basename "$0") <FASTQ_INPUT_FOLDER> <OUTPUT_DIR> [THREADS]"
    echo ""
    echo "Arguments:"
    echo "  FASTQ_INPUT_FOLDER   Path to input folder containing barcode subfolders with FASTQ files"
    echo "  OUTPUT_DIR           Path to output folder for pipeline results"
    echo "  THREADS              (Optional) Number of CPU threads to use [default: 8]"
    echo ""
    echo "Example:"
    echo "  $(basename "$0") /data/fastq_barcodes results_metagenomics 8"
    echo ""
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

nextflow run epi2me-labs/wf-metagenomics \
        --fastq $INPUT_DIR \
        --classifier minimap2 \
        --reference $SCRIPT_DIR/../referenceGenome/combined_genomes/combined_genomes.mmi \
        --ref2taxid $SCRIPT_DIR/../referenceGenome/combined_genomes/ref2taxid.targloci.tsv \
        --out_dir $OUTPUT_DIR \
        --keep_bam \
        --threads $THREADS \
        --taxonomy $SCRIPT_DIR/../referenceGenome/taxdump/new_taxdump.tar.gz \
        -process.maxForks=1

