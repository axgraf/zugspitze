#!/bin/bash
set -e

OUT_DIR="./mappings"
mkdir -p "$OUT_DIR"

### Lagopus muta
echo "Processing Lagopus muta..."
zcat Lagopus_muta/*.fna.gz \
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
  }' > "$OUT_DIR/lagopus_mapping.tsv"

### Lepus timidus
echo "Processing Lepus timidus..."
zcat Lepus_timidus/*.fna.gz \
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
  }' > "$OUT_DIR/lepus_mapping.tsv"

### Lyrurus tetrix
echo "Processing Lyrurus tetrix..."
zcat Lyrurus_tetrix/*.fna.gz \
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
  }' > "$OUT_DIR/lyrurus_mapping.tsv"


### Combine all into one file
echo "Combining mappings..."
cat "$OUT_DIR/lagopus_mapping.tsv" "$OUT_DIR/lepus_mapping.tsv" "$OUT_DIR/lyrurus_mapping.tsv" > "$OUT_DIR/combined_mapping.tsv"

echo "Final mapping file created: $OUT_DIR/combined_mapping.tsv"

