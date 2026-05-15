#!/bin/bash
set -euo pipefail

############################
# Title:    PRS Calculation in the UKBB for DKFZ
# Author:   Murat Guler
# Contact:  murat.guler@dkfz.de
# Usage:    BETAS=<betas.csv> OUT_PREFIX=<name> OUT_DIR=<outdir> bash dkfz_ukbb_prs.sh
############################

############################
#     SOFTWARE MODULES     #
############################

module load SQLite/3.46.0-GCCcore-14.1.0 PLINK/2.00a6_amd_avx2 Micromamba/2.0.2-0

# Just one-time setup of micromamba environment for bgenix.
# If bgenix is already available, this environment setup is harmless.
eval "$(micromamba shell hook --shell bash)"

if ! micromamba env list | grep -q bgen_env; then
  echo "bgen_env environment not found, setting up micromamba environment..."
  micromamba create --quiet --yes --name bgen_env
  micromamba activate bgen_env
  micromamba install --yes conda-forge::bgenix
else
  echo "bgen_env environment already exists. Activating it..."
  micromamba activate bgen_env
fi

# Check if all required software is available
for cmd in sqlite3 plink2 bgenix cat-bgen awk sort comm wc head; do
  if ! command -v "$cmd" &> /dev/null; then
    echo "ERROR: $cmd is not installed or not in PATH. Please install it before running the script."
    exit 1
  fi
done

############################
#     USER CONFIGURATION   #
############################

: "${BETAS:?BETAS is not set}"
: "${OUT_PREFIX:?OUT_PREFIX is not set}"
: "${OUT_DIR:?OUT_DIR is not set}"

############################
# UKBB DATA CONFIGURATION  #
############################

UKBB_BASE="/omics/odcf/analysis/OE0540_projects/ukkb_joint_oe0136/669373/Genomics/Imputation/Imputation_from_genotype"
BGEN_PREFIX="${UKBB_BASE}/ukb22828_c"
MFI_PREFIX="${UKBB_BASE}/ukb_mfi_chr"
SAMPLE_FILE="${UKBB_BASE}/ukb22828_c1_b0_v3_s487188.sample"

# UKB-provided high-quality individuals from PCA QC
KEEPFAM="${UKBB_BASE}/usedinpca.txt"

############################
#   OUTPUT CONFIGURATION   #
############################

mkdir -p "$OUT_DIR"

RAW_OUT="${OUT_DIR}/${OUT_PREFIX}_raw"
SNPQC_OUT="${OUT_DIR}/${OUT_PREFIX}_snpQC"
SAMPLEQC_OUT="${OUT_DIR}/${OUT_PREFIX}_sampleQC"
PRS_OUT="${OUT_DIR}/${OUT_PREFIX}_PRS"

SCORE_FILE="${OUT_DIR}/${OUT_PREFIX}_score.txt"
EFFECT_ALLELE_FILE="${OUT_DIR}/${OUT_PREFIX}_effect_alleles_for_export.txt"

BGENTMP="${OUT_DIR}/${OUT_PREFIX}_initial_chr"
FINAL_BGEN="${OUT_DIR}/${OUT_PREFIX}_single_allelic"

MATCHED_VARIANTS="${OUT_DIR}/${OUT_PREFIX}_matched_variants.tsv"
UNMATCHED_BETAS="${OUT_DIR}/${OUT_PREFIX}_unmatched_betas.tsv"
ORIENTATION_COUNTS="${OUT_DIR}/${OUT_PREFIX}_orientation_counts.tsv"

SCORE_IDS="${OUT_DIR}/${OUT_PREFIX}_score_ids.txt"
RAW_PVAR_IDS="${OUT_DIR}/${OUT_PREFIX}_raw_pvar_ids.txt"
SNPQC_IDS="${OUT_DIR}/${OUT_PREFIX}_snpqc_ids.txt"

############################
#       SANITY CHECK       #
############################

for var in BETAS BGEN_PREFIX MFI_PREFIX SAMPLE_FILE KEEPFAM OUT_PREFIX OUT_DIR; do
  if [ -z "${!var}" ]; then
    echo "ERROR: Variable $var is not set."
    exit 1
  fi
done

if [ ! -f "$BETAS" ]; then
  echo "ERROR: BETAS file not found: $BETAS"
  exit 1
fi

if [ ! -f "$SAMPLE_FILE" ]; then
  echo "ERROR: SAMPLE_FILE not found: $SAMPLE_FILE"
  exit 1
fi

if [ ! -f "$KEEPFAM" ]; then
  echo "ERROR: KEEPFAM file not found: $KEEPFAM"
  exit 1
fi

############################
#        MAIN LOGIC        #
############################

echo "Step 1: Preparing SNP lists from beta file..."

# Expected beta file columns:
# rsid,chr_name,chr_position,effect_allele,noneffect_allele,Beta,eaf,chr_pos

awk -F, '
  NR > 1 {
    gsub(/\r/, "", $1)
    if ($1 != "") print $1
  }
' "$BETAS" > "${OUT_DIR}/rsidlist.txt"

awk -F, '
  NR > 1 {
    gsub(/\r/, "", $2)
    gsub(/\r/, "", $3)
    if ($2 != "" && $3 != "") {
      print sprintf("%02d", $2) ":" $3 "-" $3
    }
  }
' "$BETAS" > "${OUT_DIR}/chrposlist.txt"

echo "Step 2: Extracting candidate variants from UKB BGEN files..."

# Extract by position. The SQL join below enforces allele matching.
# This avoids losing variants due to rsID aliases/differences.
bgen_parts=()

for chr in {1..22}; do
  bgen_file="${BGEN_PREFIX}${chr}_b0_v3.bgen"
  out_bgen="${OUT_DIR}/chr_${chr}.bgen"

  if [ ! -f "$bgen_file" ]; then
    echo "ERROR: BGEN file not found: $bgen_file"
    exit 1
  fi

  bgenix -g "$bgen_file" \
         -incl-range "${OUT_DIR}/chrposlist.txt" \
         > "$out_bgen"

  bgen_parts+=("$out_bgen")
done

echo "Step 3: Combining extracted BGEN chunks..."

cat-bgen -g "${bgen_parts[@]}" -og "${BGENTMP}.bgen" -clobber
bgenix -g "${BGENTMP}.bgen" -index -clobber

rm -f "${OUT_DIR}/chr_"*.bgen

echo "Step 4: Importing beta file into BGEN index database..."

sqlite3 "${BGENTMP}.bgen.bgi" <<SQL
DROP TABLE IF EXISTS Betas;

CREATE TABLE Betas (
  rsid TEXT,
  chr_name INTEGER,
  chr_position INTEGER,
  effect_allele TEXT,
  noneffect_allele TEXT,
  Beta REAL,
  eaf REAL,
  chr_pos TEXT
);

.mode csv
.import --skip 1 "$BETAS" Betas

UPDATE Betas
SET
  rsid = trim(rsid),
  effect_allele = upper(trim(effect_allele)),
  noneffect_allele = upper(trim(noneffect_allele));
SQL

echo "Step 5: Matching BGEN variants to beta file with correct effect-allele orientation..."

sqlite3 "${BGENTMP}.bgen.bgi" <<SQL
DROP TABLE IF EXISTS Joined;

CREATE TABLE Joined AS
SELECT
    Variant.*,
    Betas.rsid AS beta_rsid,
    Betas.chr_name AS beta_chr_name,
    Betas.chr_position AS beta_chr_position,
    Betas.effect_allele AS effect_allele,
    Betas.noneffect_allele AS noneffect_allele,
    CAST(Betas.Beta AS REAL) AS Beta,
    CAST(Betas.eaf AS REAL) AS beta_eaf,
    CASE
      WHEN Variant.allele1 = Betas.noneffect_allele
       AND Variant.allele2 = Betas.effect_allele
      THEN 'effect_is_allele2'
      WHEN Variant.allele1 = Betas.effect_allele
       AND Variant.allele2 = Betas.noneffect_allele
      THEN 'effect_is_allele1'
    END AS orientation
FROM Variant
INNER JOIN Betas
  ON (
       Variant.chromosome = CAST(Betas.chr_name AS TEXT)
       OR Variant.chromosome = printf('%02d', Betas.chr_name)
     )
 AND Variant.position = Betas.chr_position
 AND (
      (
        Variant.allele1 = Betas.noneffect_allele
        AND Variant.allele2 = Betas.effect_allele
      )
      OR
      (
        Variant.allele1 = Betas.effect_allele
        AND Variant.allele2 = Betas.noneffect_allele
      )
 );
SQL

N_JOINED=$(sqlite3 "${BGENTMP}.bgen.bgi" "SELECT COUNT(*) FROM Joined;")
N_BETAS=$(sqlite3 "${BGENTMP}.bgen.bgi" "SELECT COUNT(*) FROM Betas;")

echo "Matched variants: ${N_JOINED} / ${N_BETAS}"

if [ "$N_JOINED" -eq 0 ]; then
  echo "ERROR: No variants matched between BGEN and beta file."
  exit 1
fi

echo "Step 6: Writing audit tables..."

sqlite3 -header -separator $'\t' "${BGENTMP}.bgen.bgi" "
SELECT
  orientation,
  COUNT(*) AS n_variants
FROM Joined
GROUP BY orientation
ORDER BY orientation;
" > "$ORIENTATION_COUNTS"

sqlite3 -header -separator $'\t' "${BGENTMP}.bgen.bgi" "
SELECT
  beta_rsid AS rsid,
  beta_chr_name AS chr_name,
  beta_chr_position AS chr_position,
  allele1 AS bgen_allele1,
  allele2 AS bgen_allele2,
  effect_allele,
  noneffect_allele,
  Beta,
  beta_eaf,
  orientation,
  CAST(beta_chr_name AS TEXT) || ':' || position || '_' || allele1 || '_' || allele2 AS expected_plink_variant_id
FROM Joined
ORDER BY beta_chr_name, beta_chr_position;
" > "$MATCHED_VARIANTS"

sqlite3 -header -separator $'\t' "${BGENTMP}.bgen.bgi" "
SELECT
  Betas.rsid,
  Betas.chr_name,
  Betas.chr_position,
  Betas.effect_allele,
  Betas.noneffect_allele,
  Betas.Beta,
  Betas.eaf,
  Betas.chr_pos
FROM Betas
WHERE NOT EXISTS (
  SELECT 1
  FROM Joined
  WHERE Joined.beta_rsid = Betas.rsid
    AND Joined.beta_chr_name = Betas.chr_name
    AND Joined.beta_chr_position = Betas.chr_position
);
" > "$UNMATCHED_BETAS"

echo "Orientation counts written to: $ORIENTATION_COUNTS"
echo "Matched variants written to:   $MATCHED_VARIANTS"
echo "Unmatched betas written to:    $UNMATCHED_BETAS"

echo "Step 7: Filtering BGEN to matched alleles..."

bgenix -g "${BGENTMP}.bgen" -table Joined > "${FINAL_BGEN}.bgen"
bgenix -g "${FINAL_BGEN}.bgen" -index -clobber

echo "Step 8a: Converting BGEN to hard-call-only PLINK2 PGEN..."

plink2 --bgen "${FINAL_BGEN}.bgen" ref-first \
  --hard-call-threshold 0.1 \
  --sample "$SAMPLE_FILE" \
  --memory 15000 \
  --set-all-var-ids @:#_\$r_\$a \
  --make-pgen fill-missing-from-dosage erase-dosage \
  --out "${RAW_OUT}"

echo "Step 8b: Calculating allele frequencies from hard-call-only PGEN..."

plink2 --pfile "${RAW_OUT}" \
  --memory 15000 \
  --freq \
  --out "${RAW_OUT}"

echo "Step 9: Identifying strand-ambiguous SNPs near allele frequency 0.5..."

# .afreq columns are expected to include:
# #CHROM ID REF ALT ALT_FREQS OBS_CT
# Print only variant IDs for --exclude.


awk 'NF >= 5 && $1 !~ /^#/ {ref=toupper($3); alt=toupper($4); af=$5+0; if (af > 0.4 && af < 0.6 && ((ref=="A" && alt=="T") || (ref=="T" && alt=="A") || (ref=="C" && alt=="G") || (ref=="G" && alt=="C"))) print $2}' \
  "${RAW_OUT}.afreq" > "${OUT_DIR}/exclrsIDs_ambiguous.txt"

echo "Ambiguous SNPs excluded: $(wc -l < "${OUT_DIR}/exclrsIDs_ambiguous.txt")"

echo "Step 10: Preparing unified UKB MFI INFO table..."

# UKB MFI v3 is expected to have 8 columns.
# After prepending chr, INFO is column 9 and generated variant ID is column 10.
# PLINK --extract-col-cond syntax:
# --extract-col-cond <file> <value-column> <ID-column>
for chr in {1..22}; do
  mfi_file="${MFI_PREFIX}${chr}_v3.txt"

  if [ ! -f "$mfi_file" ]; then
    echo "ERROR: MFI file not found: $mfi_file"
    exit 1
  fi

  awk -v chr=$chr 'BEGIN {FS="\t"; OFS="\t"} { print chr, $0, chr":"$3"_"$4"_"$5 }' \
    "$mfi_file"
done > "${OUT_DIR}/ukb_mfi_all_v3.tsv"

echo "Step 11: SNP QC..."

plink2 --pfile "${RAW_OUT}" \
  --memory 15000 \
  --exclude "${OUT_DIR}/exclrsIDs_ambiguous.txt" \
  --extract-col-cond "${OUT_DIR}/ukb_mfi_all_v3.tsv" 9 10 \
  --extract-col-cond-min 0.4 \
  --maf 0.005 \
  --write-snplist \
  --make-pgen \
  --out "${SNPQC_OUT}"

echo "Step 12: Sample QC..."

plink2 --pfile "${RAW_OUT}" \
  --memory 15000 \
  --extract "${SNPQC_OUT}.snplist" \
  --keep-fam "$KEEPFAM" \
  --write-samples \
  --out "${SAMPLEQC_OUT}"

echo "Step 13: Creating corrected PLINK score file..."

# CRITICAL FIX:
# Always score effect_allele with positive Beta.
# Do NOT score allele2 with flipped Beta.
{
  echo "ID ALLELE WEIGHTED UNWEIGHTED"
  sqlite3 -separator " " -noheader "${BGENTMP}.bgen.bgi" "
  SELECT
      CAST(beta_chr_name AS TEXT) || ':' || position || '_' || allele1 || '_' || allele2,
      effect_allele,
      printf('%.12g', Beta),
      1
  FROM Joined
  ORDER BY beta_chr_name, position;
  "
} > "$SCORE_FILE"

# Useful for later PLINK dosage export / Excel validation.
sqlite3 -separator " " -noheader "${BGENTMP}.bgen.bgi" "
SELECT
    CAST(beta_chr_name AS TEXT) || ':' || position || '_' || allele1 || '_' || allele2,
    effect_allele
FROM Joined
ORDER BY beta_chr_name, position;
" > "$EFFECT_ALLELE_FILE"

if [ "$(awk 'NR > 1 {n++} END {print n+0}' "$SCORE_FILE")" -eq 0 ]; then
  echo "ERROR: Score file has no variant rows."
  exit 1
fi

echo "Score file written to:          $SCORE_FILE"
echo "Effect allele export file:      $EFFECT_ALLELE_FILE"

echo "Step 14: Checking score-file variant IDs against RAW PLINK .pvar IDs..."

awk 'NR > 1 {print $1}' "$SCORE_FILE" | sort -u > "$SCORE_IDS"
awk 'NR > 1 && $0 !~ /^#/ {print $3}' "${RAW_OUT}.pvar" | sort -u > "$RAW_PVAR_IDS"

N_SCORE_IDS=$(wc -l < "$SCORE_IDS")
N_MATCHING_RAW_IDS=$(comm -12 "$SCORE_IDS" "$RAW_PVAR_IDS" | wc -l)

echo "Score variants: $N_SCORE_IDS"
echo "Score variants found in RAW pvar: $N_MATCHING_RAW_IDS"

if [ "$N_MATCHING_RAW_IDS" -eq 0 ]; then
  echo "ERROR: None of the score-file variant IDs match ${RAW_OUT}.pvar."
  echo "Check chromosome formatting and allele order in --set-all-var-ids versus SCORE_FILE."
  exit 1
fi

echo "Step 15: Checking score-file variant IDs against SNP-QC variant list..."

sort -u "${SNPQC_OUT}.snplist" > "$SNPQC_IDS"

N_SCORE_IDS_AFTER_QC=$(comm -12 "$SCORE_IDS" "$SNPQC_IDS" | wc -l)

echo "Score variants: $N_SCORE_IDS"
echo "Score variants surviving SNP QC: $N_SCORE_IDS_AFTER_QC"

if [ "$N_SCORE_IDS_AFTER_QC" -eq 0 ]; then
  echo "ERROR: None of the score-file variants survived SNP QC."
  echo "Check variant ID format, MFI matching, ambiguous SNP exclusion, and MAF filter."
  exit 1
fi

echo "Step 16: Calculating weighted and unweighted PRS..."

plink2 --pfile "${RAW_OUT}" \
  --memory 15000 \
  --extract "${SNPQC_OUT}.snplist" \
  --keep "${SAMPLEQC_OUT}.id" \
  --score "${SCORE_FILE}" 1 2 header-read no-mean-imputation ignore-dup-ids list-variants cols=+scoresums \
  --score-col-nums 3-4 \
  --out "${PRS_OUT}"

echo "PRS calculation complete."
echo "Main output: ${PRS_OUT}.sscore"

echo "PRS .sscore columns:"
head -n 1 "${PRS_OUT}.sscore"

echo ""
echo "Use these columns downstream:"
echo "  PRS_weighted   = WEIGHTED_SUM"
echo "  PRS_unweighted = UNWEIGHTED_SUM"
echo ""
echo "Audit files:"
echo "  $MATCHED_VARIANTS"
echo "  $UNMATCHED_BETAS"
echo "  $ORIENTATION_COUNTS"