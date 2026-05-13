#!/bin/bash

# List all beta files you want to process
for BETAFILE in *_beta.csv; do
  TRAIT_NAME=$(basename "$BETAFILE" _beta.csv)
  
  echo "Processing trait: $TRAIT_NAME"

  # Call the main PRS script, passing variables dynamically
  BETAS="$BETAFILE"
  OUT_PREFIX="$TRAIT_NAME"
  OUT_DIR="results/${OUT_PREFIX}"

  # You can source the main logic instead of duplicating
  export BETAS OUT_PREFIX OUT_DIR

  bash 01_dkfz_ukbb_prs_core.sh
done
