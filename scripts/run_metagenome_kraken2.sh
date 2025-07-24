#!/bin/bash

#########################################################################
# Author: Alexander Graf  (alexander.graf@lmu.de)
#
# Run Nextflow metagenome workflow using kraken2 on an input fastq folder.
#
##########################################################################


INPUT_DIR=$1
OUTPUT_DIR=$2
KRAKEN2_DB="referenceGenome/kraken2_db/k2_pluspf_16gb/"
TAXONOMY="referenceGenome/kraken2_db/taxonomy/"

if [ -z "$1" ] || [ -z "$2" ]; then
  echo "Usage: $0 <FASTQ_INPUT_FOLDER> <OUTPUT_DIR>"
  exit 1
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

nextflow run epi2me-labs/wf-metagenomics --fastq $INPUT_DIR --database $KRAKEN2_DB --taxonomy $TAXONOMY --out_dir $OUTPUT_DIR