#!/bin/bash

#########################################################################
# Author: Alexander Graf  (alexander.graf@lmu.de)
#
# Run Nextflow metagenome workflow using Minimap2 on an input fastq folder.
# Create BED coverage files for each containing barcode using mosdepth
# and generates a HTML report for genome coverage for selected species
#
# --species options are: "Lagopus muta", "Lepus timidus", "Lyrurus tetrix"
##########################################################################


INPUT_DIR=$1
OUTPUT_DIR=$2
SPECIES=$3
BIN_SIZE=${4:-1000}  # default: 1000 if not provided

GENOME_COVERAGE_BED="/home/alex/zugspitze/referenceGenome/mappings/combined_mapping.tsv"
GENOME_COVERAGE_PYTHON="/home/alex/PycharmProjects/zugspitze/genome_coverage.py"

VALID_SPECIES=("Lagopus muta" "Lepus timidus" "Lyrurus tetrix")

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$SPECIES" ]; then
  echo "Usage: $0 <FASTQ_INPUT_FOLDER> <OUTPUT_DIR> <SPECIES> [BIN_SIZE]"
  echo "Allowed species: Lagopus muta, Lepus timidus, Lyrurus tetrix"
  exit 1
fi

is_valid_species=false
for s in "${VALID_SPECIES[@]}"; do
  if [ "$SPECIES" == "$s" ]; then
    is_valid_species=true
    break
  fi
done

if [ "$is_valid_species" = false ]; then
  echo "Error: Invalid species '$SPECIES'"
  echo "Valid options: Lagopus muta, Lepus timidus, Lyrurus tetrix"
  exit 1
fi

nextflow run epi2me-labs/wf-metagenomics --fastq $INPUT_DIR --classifier minimap2 --reference referenceGenome/combined_genomes/LeTim_LagMut_LyrTet.genome.mmi --ref2taxid referenceGenome/combined_genomes/ref2taxid.targloci.tsv --out_dir $OUTPUT_DIR --keep_bam --taxonomy referenceGenome/taxdump/new_taxdump.tar.gz -process.maxForks=1

mkdir -p ${OUTPUT_DIR}/genome_coverage_analysis
for barcode_bam in ${OUTPUT_DIR}/bams/*.bam; do
        bc=$(basename "$barcode_bam" .reference.bam)
        out_folder="${OUTPUT_DIR}/genome_coverage_analysis/$bc"
        mkdir -p $out_folder
        echo "Calculate coverage from BAM file for sample: $bc"
        mosdepth -t 12 $out_folder/$bc $barcode_bam
        echo "Genome coverage analysis for: $bc"
        python3 $GENOME_COVERAGE_PYTHON \
                --mapping $GENOME_COVERAGE_BED \
                --sample $bc \
                --bed ${OUTPUT_DIR}/genome_coverage_analysis/$bc/${bc}.per-base.bed.gz \
                --bam $barcode_bam \
                --output_dir ${OUTPUT_DIR}/genome_coverage_analysis/$bc \
                --species "Lagopus muta" \
                --bin_size 1000
done