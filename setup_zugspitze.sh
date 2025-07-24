#!/bin/bash
set -e

ENV_NAME="zugspitze_metagenome"
REQUIRED_PACKAGES=("nextflow" "seqkit" "kraken2" "pysam" "pandas" "plotly" "mosdepth" "matplotlib" "minimap2")
REFERENCE_DIR="referenceGenome"
COMBINED_FASTA="${REFERENCE_DIR}/combined_genomes/combined_genomes.fasta.gz"
MMI_INDEX="${REFERENCE_DIR}/combined_genomes/combined_genomes.mmi"
REF2TAXID="${REFERENCE_DIR}/combined_genomes/ref2taxid.targloci.tsv"

echo "🧬 Setting up Zugspitze Metagenomics Environment"
echo "🔍 Checking for Conda environment: $ENV_NAME"

# Function: Check if Conda environment exists
env_exists() {
  conda env list | grep -qE "^$ENV_NAME\s"
}

# Function: Check if all required packages are installed
env_has_all_packages() {
  for pkg in "${REQUIRED_PACKAGES[@]}"; do
    if ! conda run -n "$ENV_NAME" which "$pkg" &>/dev/null; then
      return 1
    fi
  done
  return 0
}

# Check for mamba
if ! command -v mamba &>/dev/null; then
  echo "❌ Mamba not found. Please install Mamba or Miniconda first."
  exit 1
fi

# Check and (if needed) create environment
if env_exists; then
  echo "✅ Conda environment '$ENV_NAME' already exists."

  if env_has_all_packages; then
    echo "🎉 All required packages are already installed."
  else
    echo "📦 Some packages are missing – installing..."
    mamba install -n "$ENV_NAME" -c conda-forge -c bioconda "${REQUIRED_PACKAGES[@]}" --yes
  fi
else
  echo "📦 Creating new environment '$ENV_NAME'..."
  mamba create -n "$ENV_NAME" -c conda-forge -c bioconda "${REQUIRED_PACKAGES[@]}" --yes
fi

echo "🐍 Activating conda environment '$ENV_NAME'..."
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

# Docker check and install
echo "🐳 Checking for Docker..."
if ! command -v docker &>/dev/null; then
  echo "⚠️ Docker is not installed – attempting installation (Debian/Ubuntu only)..."
  sudo apt update
  sudo apt install -y docker.io
  sudo systemctl enable --now docker
  sudo usermod -aG docker $USER
  echo "👉 Please run 'newgrp docker' or log out and back in to apply group changes."
else
  echo "✅ Docker is already installed."
fi


# Reference Genomes
echo "📥 Preparing reference genomes..."

declare -A GENOMES=(
  ["Lepus_timidus"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/893/245/GCA_040893245.2_mLepTim1.1_pri/GCA_040893245.2_mLepTim1.1_pri_genomic.fna.gz"
  ["Lagopus_muta"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/343/835/GCF_023343835.1_bLagMut1_primary/GCF_023343835.1_bLagMut1_primary_genomic.fna.gz"
  ["Lyrurus_tetrix"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/043/882/375/GCA_043882375.1_ASM4388237v1/GCA_043882375.1_ASM4388237v1_genomic.fna.gz"
)

mkdir -p "${REFERENCE_DIR}/combined_genomes"

for SPECIES in "${!GENOMES[@]}"; do
  DIR="${REFERENCE_DIR}/${SPECIES}"
  URL="${GENOMES[$SPECIES]}"
  mkdir -p "$DIR"
  FILE="$DIR/$(basename "$URL")"
  if [[ ! -f "$FILE" ]]; then
    echo "🔽 Downloading $SPECIES genome..."
    wget -q -O "$FILE" "$URL"
  else
    echo "✅ Genome for $SPECIES already exists."
  fi
done


# Combine genomes
if [[ ! -f "$COMBINED_FASTA" ]]; then
  echo "🧬 Combining genomes..."
  cat \
    ${REFERENCE_DIR}/Lepus_timidus/*.fna.gz \
    ${REFERENCE_DIR}/Lagopus_muta/*.fna.gz \
    ${REFERENCE_DIR}/Lyrurus_tetrix/*.fna.gz \
    > "$COMBINED_FASTA"
else
  echo "✅ Combined genome file already exists."
fi

# Minimap2 Index
if [[ ! -f "$MMI_INDEX" ]]; then
  echo "🔧 Building Minimap2 index..."
  conda run -n "$ENV_NAME" minimap2 -t 8 -x map-ont -d "$MMI_INDEX" "$COMBINED_FASTA"
else
  echo "✅ Minimap2 index already exists."
fi


# ref2taxid
if [[ ! -f "$REF2TAXID" ]]; then
  echo "🔗 Creating ref2taxid file..."
  cd "${REFERENCE_DIR}/combined_genomes"
  zgrep "^>" ../Lepus_timidus/*.fna.gz | cut -d' ' -f1 | sed 's/^>//' | awk '{print $0 "\t62621"}' > lepus.tsv
  zgrep "^>" ../Lagopus_muta/*.fna.gz | cut -d' ' -f1 | sed 's/^>//' | awk '{print $0 "\t64668"}' > lagopus.tsv
  zgrep "^>" ../Lyrurus_tetrix/*.fna.gz | cut -d' ' -f1 | sed 's/^>//' | awk '{print $0 "\t1233216"}' > lyrurus.tsv
  cat lepus.tsv lagopus.tsv lyrurus.tsv > ref2taxid.targloci.tsv
  echo "✅ ref2taxid file created."
  cd - >/dev/null
else
  echo "✅ ref2taxid file already exists."
fi

# -----------------------------------------------
# 🗺️ Generate mapping file for genome_coverage.py
# -----------------------------------------------

echo "🗺️ Generating mapping file for genome_coverage.py..."

MAPPING_DIR="${REFERENCE_DIR}/mappings"
mkdir -p "$MAPPING_DIR"

# Lagopus muta
echo "📍 Lagopus muta..."
zcat "${REFERENCE_DIR}/Lagopus_muta/Lagopus_muta.fna.gz" \
  | grep "^>" \
  | sed 's/^>//' \
  | awk '{
      chrom = "unknown"
      for (i=1; i<NF; i++) {
          if ($i == "chromosome") {
              chrom = $i " " $(i+1)
              if ($(i+2) == "unlocalized") {
                  chrom = chrom " unlocalized"
              }
              break
          }
      }
      gsub(/,/, "", chrom)
      print $1 "\tLagopus muta\t" chrom
  }' > "$MAPPING_DIR/lagopus_mapping.tsv"

# Lepus timidus
echo "📍 Lepus timidus..."
zcat "${REFERENCE_DIR}/Lepus_timidus/Lepus_timidus.fna.gz" \
  | grep "^>" \
  | sed 's/^>//' \
  | awk '{
      chrom = "unknown"
      for (i=1; i<NF; i++) {
          if ($i == "chromosome") {
              chrom = $i " " $(i+1)
              if ($(i+2) == "unlocalized") {
                  chrom = chrom " unlocalized"
              }
              break
          }
      }
      gsub(/,/, "", chrom)
      print $1 "\tLepus timidus\t" chrom
  }' > "$MAPPING_DIR/lepus_mapping.tsv"

# Lyrurus tetrix
echo "📍 Lyrurus tetrix..."
zcat "${REFERENCE_DIR}/Lyrurus_tetrix/Lyrurus_tetrix.fna.gz" \
  | grep "^>" \
  | sed 's/^>//' \
  | awk '{
      scaf = "unknown"
      for (i=1; i<=NF; i++) {
          if ($i ~ /HRSCAF_[0-9]+/) {
              scaf = $i
              gsub(/.*HRSCAF_/, "HRSCAF_", scaf)
              gsub(/,/, "", scaf)
              break
          }
      }
      print $1 "\tLyrurus tetrix\t" scaf
  }' > "$MAPPING_DIR/lyrurus_mapping.tsv"

# Combine
cat "$MAPPING_DIR"/{lagopus_mapping.tsv,lepus_mapping.tsv,lyrurus_mapping.tsv} > "$MAPPING_DIR/combined_mapping.tsv"

echo "✅ Mapping file created at: $MAPPING_DIR/combined_mapping.tsv"

echo "🎯 Setup complete!"