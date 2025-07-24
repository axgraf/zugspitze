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
    if ! conda list -n "$ENV_NAME" "$pkg" | grep -q "^$pkg"; then
      return 1
    fi
  done
  return 0
}

prepare_minimap2_index() {
  local name=$1         # "combined_genomes" oder "all_references"
  local species_list=("${!2}")  # Array mit Artennamen
  local taxid_map=$3    # assoziatives Array mit Art → TaxID

  local OUT_DIR="${REFERENCE_DIR}/${name}"
  mkdir -p "$OUT_DIR"
  local FASTA="${OUT_DIR}/${name}.fasta.gz"
  local MMI="${OUT_DIR}/${name}.mmi"
  local REF2TAXID="${OUT_DIR}/ref2taxid.targloci.tsv"
  local MAPPING_DIR="${OUT_DIR}/mappings"
  local COMBINED_MAPPING="${MAPPING_DIR}/combined_mapping.tsv"

  echo "🧬 Building FASTA for index: $name"
  if [[ ! -f "$FASTA" ]]; then
    for SPECIES in "${species_list[@]}"; do
      cat "${REFERENCE_DIR}/${SPECIES}"/*.fna.gz
    done > "$FASTA"
  fi

  echo "🔧 Building minimap2 index for $name..."
  if [[ ! -f "$MMI" ]]; then
    minimap2 -t 8 -x map-ont -d "$MMI" "$FASTA"
  fi

  echo "🔗 Creating ref2taxid for $name..."
  if [[ ! -f "$REF2TAXID" ]]; then
    > "$REF2TAXID"
    for SPECIES in "${species_list[@]}"; do
      TAXID=${taxid_map[$SPECIES]}
      zgrep "^>" "${REFERENCE_DIR}/${SPECIES}"/*.fna.gz | cut -d' ' -f1 | sed 's/^>//' | awk -v taxid="$TAXID" '{print $0 "\t" taxid}' >> "$REF2TAXID"
    done
  fi

  echo "🗺️ Generating mapping file for $name..."
  mkdir -p "$MAPPING_DIR"
  if [[ ! -f "$COMBINED_MAPPING" ]]; then
    > "$COMBINED_MAPPING"
    for SPECIES in "${species_list[@]}"; do
      case $SPECIES in
        Lagopus_muta|Lepus_timidus)
          zcat "${REFERENCE_DIR}/${SPECIES}"/*.fna.gz \
            | grep "^>" \
            | sed 's/^>//' \
            | awk -v sname="${SPECIES//_/ }" '{
                chrom = "unknown"
                for (i=1; i<NF; i++) {
                  if ($i == "chromosome") {
                    chrom = $i " " $(i+1)
                    if ($(i+2) == "unlocalized") chrom = chrom " unlocalized"
                    break
                  }
                }
                gsub(/,/, "", chrom)
                print $1 "\t" sname "\t" chrom
            }' >> "$COMBINED_MAPPING"
          ;;
        Lyrurus_tetrix)
          zcat "${REFERENCE_DIR}/${SPECIES}"/*.fna.gz \
            | grep "^>" \
            | sed 's/^>//' \
            | awk -v sname="Lyrurus tetrix" '{
                scaf = "unknown"
                for (i=1; i<=NF; i++) {
                  if ($i ~ /HRSCAF_[0-9]+/) {
                    scaf = $i
                    gsub(/.*HRSCAF_/, "HRSCAF_", scaf)
                    gsub(/,/, "", scaf)
                    break
                  }
                }
                print $1 "\t" sname "\t" scaf
            }' >> "$COMBINED_MAPPING"
          ;;
        *)
          # Für andere: dummy mapping (optional später anpassen)
          zcat "${REFERENCE_DIR}/${SPECIES}"/*.fna.gz \
            | grep "^>" \
            | sed 's/^>//' \
            | awk -v sname="${SPECIES//_/ }" '{print $1 "\t" sname "\tunknown"}' \
            >> "$COMBINED_MAPPING"
          ;;
      esac
    done
  fi

  echo "✅ [$name] Setup done."
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

# 3-Arten-Index
DEFAULT_SPECIES=("Lepus_timidus" "Lagopus_muta" "Lyrurus_tetrix")
declare -A TAXIDS_DEFAULT=(
  ["Lepus_timidus"]=62621
  ["Lagopus_muta"]=64668
  ["Lyrurus_tetrix"]=1233216
)
prepare_minimap2_index "combined_genomes" DEFAULT_SPECIES[@] TAXIDS_DEFAULT

# 7-Arten-Index (optional)
ALL_SPECIES=("Lepus_timidus" "Lagopus_muta" "Lyrurus_tetrix" "Tetrao_urogallus" "Mustela_erminea" "Rupicapra_rupicapra" "Ovis_aries")
declare -A TAXIDS_ALL=(
  ["Lepus_timidus"]=62621
  ["Lagopus_muta"]=64668
  ["Lyrurus_tetrix"]=1233216
  ["Tetrao_urogallus"]=100830
  ["Mustela_erminea"]=36723
  ["Rupicapra_rupicapra"]=64668
  ["Ovis_aries"]=9940
)
prepare_minimap2_index "all_references" ALL_SPECIES[@] TAXIDS_ALL



echo "🧬 Preparing Kraken2 database (PlusPFP-16GB)..."

KRAKEN_DB_DIR="${REFERENCE_DIR}/kraken2_db/k2_pluspf_16gb"
KRAKEN_DB_TAR="k2_pluspf_16gb_20250402.tar.gz"
KRAKEN_DB_URL="https://genome-idx.s3.amazonaws.com/kraken/${KRAKEN_DB_TAR}"

mkdir -p "$KRAKEN_DB_DIR"

if [[ ! -f "${KRAKEN_DB_DIR}/${KRAKEN_DB_TAR}" ]]; then
  echo "📦 Downloading Kraken2 DB (~15 GB)..."
  wget -c -O "${KRAKEN_DB_DIR}/${KRAKEN_DB_TAR}" "$KRAKEN_DB_URL"
else
  echo "✅ Kraken2 DB archive already exists."
fi

# Entpacken (nur wenn noch nicht entpackt)
if [[ ! -f "${KRAKEN_DB_DIR}/hash.k2d" ]]; then
  echo "📦 Extracting Kraken2 DB..."
  tar -xf "${KRAKEN_DB_DIR}/${KRAKEN_DB_TAR}" -C "$KRAKEN_DB_DIR"
else
  echo "✅ Kraken2 DB already extracted."
fi

# Taxdump
TAXONOMY_DIR="${REFERENCE_DIR}/kraken2_db/taxonomy"
mkdir -p "$TAXONOMY_DIR"
TAXDUMP_TAR="${TAXONOMY_DIR}/new_taxdump.tar.gz"

if [[ ! -f "$TAXDUMP_TAR" ]]; then
  echo "📥 Downloading NCBI taxonomy archive..."
  wget -O "$TAXDUMP_TAR" https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
else
  echo "✅ NCBI taxonomy archive already exists."
fi

echo "🧪 Testing kraken2..."
conda run -n "$ENV_NAME" kraken2 --version

echo "✅ Kraken2 setup complete."

echo "🎯 Setup complete!"