#!/bin/bash
set -euo pipefail
shopt -s nullglob

CORE_SCRIPT="01_dkfz_ukbb_prs_core.sh"

if [ ! -f "$CORE_SCRIPT" ]; then
  echo "ERROR: Core PRS script not found: $CORE_SCRIPT"
  exit 1
fi

beta_files=(*_beta.csv)

if [ "${#beta_files[@]}" -eq 0 ]; then
  echo "ERROR: No *_beta.csv files found in current directory."
  exit 1
fi

mkdir -p results logs

for BETAFILE in "${beta_files[@]}"; do
  TRAIT_NAME=$(basename "$BETAFILE" _beta.csv)

  echo "Processing trait: $TRAIT_NAME"

  BETAS="$(realpath "$BETAFILE")"
  OUT_PREFIX="$TRAIT_NAME"
  OUT_DIR="$(realpath results)/${OUT_PREFIX}"

  mkdir -p "$OUT_DIR"

  echo "  BETAS=$BETAS"
  echo "  OUT_PREFIX=$OUT_PREFIX"
  echo "  OUT_DIR=$OUT_DIR"

  BETAS="$BETAS" \
  OUT_PREFIX="$OUT_PREFIX" \
  OUT_DIR="$OUT_DIR" \
  bash "$CORE_SCRIPT" \
    > "logs/${TRAIT_NAME}.prs.log" \
    2> "logs/${TRAIT_NAME}.prs.err"

  echo "Finished trait: $TRAIT_NAME"
done

echo "All PRS jobs completed."