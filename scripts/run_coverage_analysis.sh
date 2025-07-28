#!/bin/bash

#########################################################################
# Author: Alexander Graf  (alexander.graf@lmu.de)
#
# Create BED coverage files for each containing barcode using mosdepth
# and generates a HTML report for genome coverage for selected species
#
# --species options are: "Lagopus muta", "Lepus timidus", "Lyrurus tetrix"
##########################################################################
set -e

eval "$(mamba shell hook --shell bash)"
mamba activate zugspitze_metagenome

RED=$(tput setaf 1)
BOLD=$(tput bold)
RESET=$(tput sgr0)

INPUT_DIR=$1
OUTPUT_DIR=$2
SPECIES=$3
BIN_SIZE=${4:-10000}  # default: 10000 if not provided

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(realpath "${SCRIPT_DIR}/..")"

GENOME_COVERAGE_BED="${PROJECT_ROOT}/referenceGenome/mappings/combined_mapping.tsv"
GENOME_COVERAGE_PYTHON="${PROJECT_ROOT}/genome_coverage.py"

VALID_SPECIES=("Lagopus muta" "Lepus timidus" "Lyrurus tetrix")

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$SPECIES" ]; then
  echo ""
  echo "${RED}[ERROR] Missing required arguments.${RESET}"
  echo ""
  echo "${BOLD}Genome Coverage Analysis Tool${RESET}"
  echo "------------------------------------"
  echo ""
  echo "${BOLD}Usage:${RESET}"
  echo "  $SCRIPT_NAME <BAM_FOLDER> <OUTPUT_DIR> <SPECIES> [BIN_SIZE]"
  echo ""
  echo "${BOLD}Arguments:${RESET}"
  echo "  BAM_FOLDER     Directory containing input BAM files"
  echo "  OUTPUT_DIR     Directory where HTML report file will be saved"
  echo "  SPECIES        Target species name (in quotes)"
  echo "  BIN_SIZE       Optional: bin size in base pairs (default: 10000)"
  echo ""
  echo "${BOLD}Supported species:${RESET}"
  echo "  - Lagopus muta"
  echo "  - Lepus timidus"
  echo "  - Lyrurus tetrix"
  echo ""
  echo "${BOLD}Example:${RESET}"
  echo "  $SCRIPT_NAME ./bam_files ./output 'Lagopus muta' 5000"
  echo ""
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

echo "Running coverage analysis for all BAM files"
echo ""
echo "Species: $SPECIES"
echo ""

mkdir -p ${OUTPUT_DIR}

for barcode_bam in ${INPUT_DIR}/*.bam; do
	bc=$(basename "$barcode_bam" .reference.bam)
	out_folder="${OUTPUT_DIR}/$bc"
	mkdir -p $out_folder
	echo ""
	echo "##################################################"
	echo "Analyzing sample: $bc"
	echo "##################################################"
	echo ""
	echo -e "[mosdepth]\tProcessing sample: $bc"
       	echo -e "   Generating BED coverage from BAM file"
	mosdepth -t 6 $out_folder/$bc $barcode_bam
	echo "Starting genome coverage analysis for: $bc"
	python3 $GENOME_COVERAGE_PYTHON \
		--mapping $GENOME_COVERAGE_BED \
		--sample $bc \
		--bed ${OUTPUT_DIR}/$bc/${bc}.per-base.bed.gz \
		--bam $barcode_bam \
		--output_dir ${OUTPUT_DIR}/$bc \
		--species "$SPECIES" \
		--bin_size $BIN_SIZE
done

