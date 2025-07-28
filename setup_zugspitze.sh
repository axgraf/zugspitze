#!/bin/bash
set -e

ENV_NAME="zugspitze_metagenome"
REQUIRED_PACKAGES=("nextflow" "seqkit" "kraken2" "pysam" "pandas" "plotly" "mosdepth" "matplotlib")
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
      TAXID=$(eval "echo \${$3[$SPECIES]}")
      zgrep "^>" "${REFERENCE_DIR}/${SPECIES}"/*.fna.gz | cut -d' ' -f1 | sed 's/^>//' | awk -v taxid="$TAXID" '{print $0 "\t" taxid}' >> "$REF2TAXID"
    done
  fi

  echo "🗺️ Generating mapping file for $name..."
  mkdir -p "$MAPPING_DIR"
  if [[ ! -f "$COMBINED_MAPPING" ]]; then
    > "$COMBINED_MAPPING"  # einmalig leeren

    for SPECIES in "${species_list[@]}"; do
      echo "📍 Processing $SPECIES..."
      zcat "${REFERENCE_DIR}/${SPECIES}"/*.fna.gz \
      | grep "^>" \
      | sed 's/^>//' \
      | awk -v sname="${SPECIES//_/ }" '
        function clean(x) { gsub(/,/, "", x); gsub(/[^[:alnum:]_:\- ]/, "", x); return x }
        {
          region = "unknown"
          for (i = 1; i < NF; i++) {
            if ($i == "chromosome" || $i == "chromosome:") { region = "chromosome " $(i+1); break }
            if ($i ~ /^scaffold_/) { region = $i; break }
            if ($i == "contig:" && $(i+1) ~ /^[A-Za-z0-9_]+$/) { region = $(i+1); break }
            if ($i ~ /^HAP1_SCAFFOLD_/) { region = $i; break }
            if ($i ~ /^SUPER_/) { region = $i; break }
            if ($i == "organelle:" && $(i+1) == "mitochondrion") { region = "mitochondrion"; break }
            if ($i ~ /mitochondrion/) { region = "mitochondrion"; break }
            if ($i ~ /HRSCAF_[0-9]+/) { region = $i; break }
          }
          print $1 "\t" sname "\t" clean(region)
        }' >> "$COMBINED_MAPPING"
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

# -----------------------------------------------
# 📥 Download reference genomes (if missing)
# -----------------------------------------------

echo "🌐 Downloading reference genomes..."

declare -A GENOMES=(
  ["Lepus_timidus"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/893/245/GCA_040893245.2_mLepTim1.1_pri/GCA_040893245.2_mLepTim1.1_pri_genomic.fna.gz"
  ["Lagopus_muta"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/343/835/GCF_023343835.1_bLagMut1_primary/GCF_023343835.1_bLagMut1_primary_genomic.fna.gz"
  ["Lyrurus_tetrix"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/043/882/375/GCA_043882375.1_ASM4388237v1/GCA_043882375.1_ASM4388237v1_genomic.fna.gz"
  ["Tetrao_urogallus"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/951/394/365/GCA_951394365.1_bTetUro1.1/GCA_951394365.1_bTetUro1.1_genomic.fna.gz"
  ["Mustela_erminea"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/829/155/GCF_009829155.1_mMusErm1.Pri/GCF_009829155.1_mMusErm1.Pri_genomic.fna.gz"
  ["Rupicapra_rupicapra"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/963/981/305/GCA_963981305.1_mRupRup1.1/GCA_963981305.1_mRupRup1.1_genomic.fna.gz"
  ["Ovis_aries"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/772/045/GCF_016772045.2_ARS-UI_Ramb_v3.0/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna.gz"
)

for SPECIES in "${!GENOMES[@]}"; do
  URL="${GENOMES[$SPECIES]}"
  TARGET_DIR="${REFERENCE_DIR}/${SPECIES}"
  FILENAME="$(basename "$URL")"
  FILEPATH="${TARGET_DIR}/${FILENAME}"

  mkdir -p "$TARGET_DIR"
  if [[ ! -f "$FILEPATH" ]]; then
    echo "🔽 Downloading $SPECIES..."
    wget -q -O "$FILEPATH" "$URL"
  else
    echo "✅ $SPECIES already downloaded."
  fi
done

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