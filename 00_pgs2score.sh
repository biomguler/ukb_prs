#!/bin/bash
set -euo pipefail

#############################################################################
# Title:    Download and convert PGS Catalog scoring file to DKFZ beta format
# Author:   Murat Guler
# Usage:
#   bash 00_pgs2score.sh <PGS_ID> [BUILD] [OUT_DIR]
#
# Example:
#   bash 00_pgs2score.sh PGS000083 GRCh37 scores
#
# Output:
#   scores/PGS000083/PGS000083_beta.csv
#
# Expected output columns:
#   rsid,chr_name,chr_position,effect_allele,noneffect_allele,Beta,eaf,chr_pos
#
# Notes:
#   - BUILD defaults to GRCh37.
#   - The script downloads from the Harmonized PGS Catalog FTP directory.
#   - For harmonized files, hm_chr/hm_pos/hm_rsID are used when available.
#############################################################################

if [ "$#" -lt 1 ]; then
  echo "Usage: bash $0 <PGS_ID> [BUILD] [OUT_DIR]"
  echo ""
  echo "Example:"
  echo "  bash $0 PGS000083 GRCh37 scores"
  exit 1
fi

PGS_ID="$1"
BUILD="${2:-GRCh37}"
BASE_OUT_DIR="${3:-.}"

PGS_OUT_DIR="${BASE_OUT_DIR}/${PGS_ID}"
mkdir -p "$PGS_OUT_DIR"

SCORING_FILE="${PGS_OUT_DIR}/${PGS_ID}_hmPOS_${BUILD}.txt.gz"
OUT_BETA="${PGS_OUT_DIR}/${PGS_ID}_beta.csv"
LOG_FILE="${PGS_OUT_DIR}/${PGS_ID}_pgs2score.log"

URL="https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/${PGS_ID}/ScoringFiles/Harmonized/${PGS_ID}_hmPOS_${BUILD}.txt.gz"

#############################################################################
# Check dependencies
#############################################################################

for cmd in awk gzip; do
  if ! command -v "$cmd" &> /dev/null; then
    echo "ERROR: Required command not found: $cmd"
    exit 1
  fi
done

if command -v wget &> /dev/null; then
  DOWNLOAD_CMD="wget -q -O"
elif command -v curl &> /dev/null; then
  DOWNLOAD_CMD="curl -fsSL -o"
else
  echo "ERROR: Neither wget nor curl found. Please install one of them."
  exit 1
fi

#############################################################################
# Download
#############################################################################

echo "PGS ID:      $PGS_ID" | tee "$LOG_FILE"
echo "Build:       $BUILD" | tee -a "$LOG_FILE"
echo "URL:         $URL" | tee -a "$LOG_FILE"
echo "Output dir:  $PGS_OUT_DIR" | tee -a "$LOG_FILE"

if [ ! -s "$SCORING_FILE" ]; then
  echo "Downloading scoring file..." | tee -a "$LOG_FILE"

  if command -v wget &> /dev/null; then
    wget -q -O "$SCORING_FILE" "$URL"
  else
    curl -fsSL -o "$SCORING_FILE" "$URL"
  fi
else
  echo "Scoring file already exists, skipping download:" | tee -a "$LOG_FILE"
  echo "  $SCORING_FILE" | tee -a "$LOG_FILE"
fi

if [ ! -s "$SCORING_FILE" ]; then
  echo "ERROR: Downloaded scoring file is empty or missing: $SCORING_FILE"
  exit 1
fi

if ! gzip -t "$SCORING_FILE"; then
  echo "ERROR: Downloaded file is not a valid gzip file: $SCORING_FILE"
  echo "Check whether this PGS/build exists in the Harmonized directory."
  exit 1
fi

#############################################################################
# Convert
#############################################################################

echo "Converting scoring file to beta format..." | tee -a "$LOG_FILE"

gzip -cd "$SCORING_FILE" | awk -F'\t' -v OFS=',' '
BEGIN {
  print "rsid","chr_name","chr_position","effect_allele","noneffect_allele","Beta","eaf","chr_pos"
}

# Skip metadata lines
/^#/ {
  next
}

# Header line
NR > 1 && $0 !~ /^#/ && header_seen == 0 {
  header_seen = 1

  for (i = 1; i <= NF; i++) {
    col[$i] = i
  }

  required[1] = "rsID"
  required[2] = "chr_name"
  required[3] = "chr_position"
  required[4] = "effect_allele"
  required[5] = "other_allele"
  required[6] = "effect_weight"
  required[7] = "allelefrequency_effect"

  for (j = 1; j <= 7; j++) {
    if (!(required[j] in col)) {
      print "ERROR: Required column missing from PGS file: " required[j] > "/dev/stderr"
      exit 1
    }
  }

  has_hm_rsID = ("hm_rsID" in col)
  has_hm_chr  = ("hm_chr" in col)
  has_hm_pos  = ("hm_pos" in col)

  next
}

# Data lines
header_seen == 1 {
  rsid = $col["rsID"]
  chr  = $col["chr_name"]
  pos  = $col["chr_position"]

  if (has_hm_rsID && $col["hm_rsID"] != "") {
    rsid = $col["hm_rsID"]
  }

  if (has_hm_chr && $col["hm_chr"] != "") {
    chr = $col["hm_chr"]
  }

  if (has_hm_pos && $col["hm_pos"] != "") {
    pos = $col["hm_pos"]
  }

  effect_allele = toupper($col["effect_allele"])
  noneffect_allele = toupper($col["other_allele"])
  beta = $col["effect_weight"]
  eaf = $col["allelefrequency_effect"]

  gsub(/^chr/, "", chr)
  gsub(/^CHR/, "", chr)

  # Remove accidental CR characters
  gsub(/\r/, "", rsid)
  gsub(/\r/, "", chr)
  gsub(/\r/, "", pos)
  gsub(/\r/, "", effect_allele)
  gsub(/\r/, "", noneffect_allele)
  gsub(/\r/, "", beta)
  gsub(/\r/, "", eaf)

  # Drop incomplete rows that cannot be used by the PRS pipeline
  if (rsid == "" || chr == "" || pos == "" || effect_allele == "" || noneffect_allele == "" || beta == "") {
    skipped++
    next
  }

  print rsid, chr, pos, effect_allele, noneffect_allele, beta, eaf, chr ":" pos
  written++
}

END {
  if (header_seen != 1) {
    print "ERROR: No header line found in PGS scoring file." > "/dev/stderr"
    exit 1
  }

  print "Rows written: " written > "/dev/stderr"
  print "Rows skipped due to missing required values: " skipped + 0 > "/dev/stderr"
}
' > "$OUT_BETA" 2>> "$LOG_FILE"

#############################################################################
# Validate output
#############################################################################

N_ROWS=$(awk 'NR > 1 {n++} END {print n + 0}' "$OUT_BETA")

if [ "$N_ROWS" -eq 0 ]; then
  echo "ERROR: Converted beta file has no variant rows: $OUT_BETA" | tee -a "$LOG_FILE"
  exit 1
fi

echo "Conversion complete." | tee -a "$LOG_FILE"
echo "Rows in beta file: $N_ROWS" | tee -a "$LOG_FILE"
echo "Output beta file:" | tee -a "$LOG_FILE"
echo "  $OUT_BETA" | tee -a "$LOG_FILE"

echo ""
echo "Preview:"
head "$OUT_BETA"