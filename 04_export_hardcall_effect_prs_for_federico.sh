#!/bin/bash
set -euo pipefail

#############################################################################
# Title:    Download and convert PGS Catalog scoring file to DKFZ beta format
# Author:   Murat Guler
# Usage:
#   bash 00_pgs2score.sh <PGS_ID> [BUILD] [OUT_DIR]
#
# Examples:
#
# Single PGS:
#   bash 00_pgs2score.sh PGS000083 GRCh37 scores
#
# Multiple PGS IDs:
#   for PGS_ID in PGS000083 PGS000018 PGS000021; do
#     bash 00_pgs2score.sh "$PGS_ID" GRCh37 scores
#   done
#
# Multiple PGS IDs from a text file, one ID per line:
#   while read -r PGS_ID; do
#     [ -z "$PGS_ID" ] && continue
#     bash 00_pgs2score.sh "$PGS_ID" GRCh37 scores
#   done < pgs_ids.txt
#
# Output:
#   scores/<PGS_ID>/<PGS_ID>_beta.csv
#
# Expected output columns:
#   rsid,chr_name,chr_position,effect_allele,noneffect_allele,Beta,eaf,chr_pos
#############################################################################
module load PLINK/2.00a6_amd_avx2
if [ "$#" -lt 3 ]; then
  echo "Usage: bash $0 <SNPQC_PFILE_PREFIX> <SCORE_FILE> <OUT_PREFIX> [KEEP_FILE]"
  exit 1
fi

PFILE_PREFIX="$1"
SCORE_FILE="$2"
OUT_PREFIX="$3"
KEEP_FILE="${4:-}"

OUT_DIR=$(dirname "$OUT_PREFIX")
mkdir -p "$OUT_DIR"

for cmd in plink2 awk sort wc; do
  if ! command -v "$cmd" &> /dev/null; then
    echo "ERROR: Required command not found: $cmd"
    exit 1
  fi
done

for ext in pgen pvar psam; do
  if [ ! -f "${PFILE_PREFIX}.${ext}" ]; then
    echo "ERROR: Missing PGEN file: ${PFILE_PREFIX}.${ext}"
    exit 1
  fi
done

if [ ! -f "$SCORE_FILE" ]; then
  echo "ERROR: SCORE_FILE not found: $SCORE_FILE"
  exit 1
fi

if [ -n "$KEEP_FILE" ] && [ ! -f "$KEEP_FILE" ]; then
  echo "ERROR: KEEP_FILE not found: $KEEP_FILE"
  exit 1
fi

HARD_PREFIX="${OUT_PREFIX}_tmp_hardcalls"
RAW_PREFIX="${OUT_PREFIX}_effect_allele_matrix"

SCORE_FILTERED="${OUT_PREFIX}_score_filtered_to_pgen.txt"
SCORE_IDS="${OUT_PREFIX}_score_variant_ids.txt"
EFFECT_ALLELE_FILE="${OUT_PREFIX}_effect_alleles.txt"
BETA_FILE="${OUT_PREFIX}_betas.txt"

OUT_TSV="${OUT_PREFIX}.tsv"

echo "Input PGEN prefix: $PFILE_PREFIX"
echo "Input score file:  $SCORE_FILE"
echo "Output prefix:     $OUT_PREFIX"

#############################################################################
# 1. Filter score file to variants present in the PGEN
#############################################################################

echo "Filtering score file to variants present in ${PFILE_PREFIX}.pvar..."

awk '
  BEGIN { OFS = "\t" }

  NR == FNR {
    if ($0 !~ /^#/ && NF >= 3) {
      keep[$3] = 1
    }
    next
  }

  FNR == 1 {
    print $0
    next
  }

  FNR > 1 && ($1 in keep) {
    print $0
  }
' "${PFILE_PREFIX}.pvar" "$SCORE_FILE" > "$SCORE_FILTERED"

N_SCORE_TOTAL=$(awk 'NR > 1 {n++} END {print n + 0}' "$SCORE_FILE")
N_SCORE_FILTERED=$(awk 'NR > 1 {n++} END {print n + 0}' "$SCORE_FILTERED")

echo "Variants in original score file: $N_SCORE_TOTAL"
echo "Variants overlapping PGEN:       $N_SCORE_FILTERED"

if [ "$N_SCORE_FILTERED" -eq 0 ]; then
  echo "ERROR: No score-file variants overlap the PGEN variant IDs."
  exit 1
fi

awk 'NR > 1 {print $1}' "$SCORE_FILTERED" > "$SCORE_IDS"

#############################################################################
# 2. Create effect-allele and beta helper files
#############################################################################

# For PLINK --export-allele:
# column 1 = variant ID
# column 2 = allele to count
awk 'NR > 1 {print $1, $2}' "$SCORE_FILTERED" > "$EFFECT_ALLELE_FILE"

# For AWK calculation:
# column 1 = variant ID
# column 2 = beta
awk 'NR > 1 {print $1, $3}' "$SCORE_FILTERED" > "$BETA_FILE"

#############################################################################
# 3. Make a dosage-free hard-call PGEN
#############################################################################

echo "Creating hard-call-only PGEN with erase-dosage..."

PLINK_KEEP_ARGS=()
if [ -n "$KEEP_FILE" ]; then
  PLINK_KEEP_ARGS=(--keep "$KEEP_FILE")
  echo "Applying sample keep file: $KEEP_FILE"
fi

plink2 --pfile "$PFILE_PREFIX" \
  "${PLINK_KEEP_ARGS[@]}" \
  --extract "$SCORE_IDS" \
  --make-pgen fill-missing-from-dosage erase-dosage \
  --out "$HARD_PREFIX"

#############################################################################
# 4. Export hard-called additive effect-allele counts
#############################################################################

echo "Exporting hard-called effect-allele dosage matrix..."

plink2 --pfile "$HARD_PREFIX" \
  --export A \
  --export-allele "$EFFECT_ALLELE_FILE" \
  --out "$RAW_PREFIX"

if [ ! -f "${RAW_PREFIX}.raw" ]; then
  echo "ERROR: Expected PLINK raw export not found: ${RAW_PREFIX}.raw"
  exit 1
fi

#############################################################################
# 5. Calculate per-individual unweighted and weighted hard-call PRS
#############################################################################

echo "Calculating hardcall_effect_allele_sum and hardcall_effect_allele_sumxbeta..."

awk '
  BEGIN {
    OFS = "\t"
  }

  NR == FNR {
    beta[$1] = $2 + 0
    next
  }

  FNR == 1 {
    iid_col = 0

    for (i = 1; i <= NF; i++) {
      if ($i == "IID") {
        iid_col = i
      }
    }

    if (iid_col == 0) {
      print "ERROR: IID column not found in PLINK .raw header." > "/dev/stderr"
      exit 1
    }

    n_score_cols = 0

    for (i = 7; i <= NF; i++) {
      raw_id = $i
      variant_id = raw_id

      # PLINK .raw columns may be either:
      #   variant_id
      # or
      #   variant_id_countedAllele
      #
      # Since our variant IDs already contain underscores, remove only the final
      # underscore-delimited allele if exact matching fails.
      if (!(variant_id in beta)) {
        stripped_id = raw_id
        sub(/_[^_]+$/, "", stripped_id)

        if (stripped_id in beta) {
          variant_id = stripped_id
        }
      }

      if (variant_id in beta) {
        score_col[i] = variant_id
        n_score_cols++
      } else {
        print "WARNING: No beta found for raw column: " raw_id > "/dev/stderr"
      }
    }

    if (n_score_cols == 0) {
      print "ERROR: No genotype columns in .raw could be matched to beta file." > "/dev/stderr"
      exit 1
    }

    print "IID", "hardcall_effect_allele_sum", "hardcall_effect_allele_sumxbeta"
    next
  }

  FNR > 1 {
    iid = $iid_col
    unweighted_sum = 0
    weighted_sum = 0

    for (i = 7; i <= NF; i++) {
      if (i in score_col) {
        g = $i

        # Missing hard calls are treated as 0 contribution,
        # matching no-mean-imputation behavior.
        if (g == "NA" || g == "." || g == "") {
          g = 0
        }

        # Sanity check: hard calls should be 0/1/2 after erase-dosage.
        if (g != 0 && g != 1 && g != 2) {
          non_integer_count++
        }

        variant_id = score_col[i]
        unweighted_sum += g
        weighted_sum += g * beta[variant_id]
      }
    }

    print iid, unweighted_sum, weighted_sum
  }

  END {
    if (non_integer_count > 0) {
      print "WARNING: Found " non_integer_count " non-0/1/2 genotype values in raw export." > "/dev/stderr"
    }
  }
' "$BETA_FILE" "${RAW_PREFIX}.raw" > "$OUT_TSV"

#############################################################################
# 6. Final report
#############################################################################

N_OUT=$(awk 'NR > 1 {n++} END {print n + 0}' "$OUT_TSV")

echo "Done."
echo "Output file:"
echo "  $OUT_TSV"
echo "Individuals written:"
echo "  $N_OUT"
echo ""
echo "Columns:"
echo "  IID"
echo "  hardcall_effect_allele_sum"
echo "  hardcall_effect_allele_sumxbeta"
echo ""
echo "Intermediate hard-call matrix:"
echo "  ${RAW_PREFIX}.raw"