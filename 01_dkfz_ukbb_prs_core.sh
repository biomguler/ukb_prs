#!/bin/bash
############################
# Title:    PRS Calculation in the UKBB for DKFZ
# Author:   Murat Guler
# Contact:  murat.guler@dkfz.de
# Usage:    bash dkfz_ukbb_prs.sh

############################
#     SOFTWARE MODULES     #
############################
# Available on CentOS clusters. Adjust for Debian if needed.
module load bgen/1.1.7 sqlite/3.38.5 plink/2.0_alpha5.12

############################
#     USER CONFIGURATION   #
############################

: "${BETAS:?BETAS is not set}"
: "${OUT_PREFIX:?OUT_PREFIX is not set}"
: "${OUT_DIR:?OUT_DIR is not set}"


############################
# UKBB DATA CONFIGURATION  #
############################

# Do not change input config if you are using ukbb data!

UKBB_BASE="/omics/odcf/analysis/OE0540_projects/ukkb_joint_oe0136/669373/Genomics/Imputation/Imputation_from_genotype"
BGEN_PREFIX="${UKBB_BASE}/ukb22828_c"       # e.g. ukb22828_c1_b0_v3.bgen
MFI_PREFIX="${UKBB_BASE}/ukb_mfi_chr"       # e.g. ukb_mfi_chr1_v3.txt
SAMPLE_FILE="${UKBB_BASE}/ukb22828_c1_b0_v3_s487188.sample"

# UKB-provided high-quality individuals (from PCA QC)
KEEPFAM="${UKBB_BASE}/usedinpca.txt"

############################
#   OUTPUT CONFIGURATION   #
############################

# Output file paths
RAW_OUT="${OUT_DIR}/${OUT_PREFIX}_raw"
SNPQC_OUT="${OUT_DIR}/${OUT_PREFIX}_snpQC"
SAMPLEQC_OUT="${OUT_DIR}/${OUT_PREFIX}_sampleQC"
PRS_OUT="${OUT_DIR}/${OUT_PREFIX}_PRS"
SCORE_FILE="${OUT_DIR}/${OUT_PREFIX}_score.txt"
BGENTMP="${OUT_DIR}/${OUT_PREFIX}_initial_chr"
FINAL_BGEN="${OUT_DIR}/${OUT_PREFIX}_single_allelic"

############################
#       SANITY CHECK       #
############################

for var in BETAS BGEN_PREFIX MFI_PREFIX SAMPLE_FILE KEEPFAM OUT_PREFIX OUT_DIR; do
  if [ -z "${!var}" ]; then
    echo "ERROR: Variable $var is not set. Please edit the script and set it at the top."
    exit 1
  fi
done

mkdir -p "$OUT_DIR"

############################
#        MAIN LOGIC        #
############################

# Extract SNPs of interest
awk -F, 'NR>1 { print $1 }' "$BETAS" > "${OUT_DIR}/rsidlist.txt"
awk -F, '{ if (NR>1) { print sprintf("%02d", $2)":"$3"-"$3 }}' "$BETAS" > "${OUT_DIR}/chrposlist.txt"

# Extract relevant SNPs from each chromosome's BGEN
cmd=""
for chr in {1..22}; do
  bgen_file="${BGEN_PREFIX}${chr}_b0_v3.bgen"
  out_bgen="${OUT_DIR}/chr_${chr}.bgen"
  bgenix -g "$bgen_file" \
         -incl-rsids "${OUT_DIR}/rsidlist.txt" \
         -incl-range "${OUT_DIR}/chrposlist.txt" > "$out_bgen"
  cmd="$cmd $out_bgen"
done

# Combine selected variants across chromosomes
cat-bgen -g $cmd -og "${BGENTMP}.bgen" -clobber
bgenix -g "${BGENTMP}.bgen" -index -clobber

# Clean up intermediate per-chromosome BGENs
rm "${OUT_DIR}/chr_"*.bgen

# Import betas and match with variants
sqlite3 "${BGENTMP}.bgen.bgi" "DROP TABLE IF EXISTS Betas;"
sqlite3 -separator "," "${BGENTMP}.bgen.bgi" ".import $BETAS Betas"

sqlite3 "${BGENTMP}.bgen.bgi" "DROP TABLE IF EXISTS Joined;"
sqlite3 -header -csv "${BGENTMP}.bgen.bgi" "
CREATE TABLE Joined AS 
 SELECT Variant.*, Betas.chr_name, Betas.Beta FROM Variant INNER JOIN Betas 
  ON Variant.chromosome = printf('%02d', Betas.chr_name)  
  AND Variant.position = Betas.chr_position 
  AND Variant.allele1 = Betas.noneffect_allele 
  AND Variant.allele2 = Betas.effect_allele 
 UNION 
 SELECT Variant.*, Betas.chr_name, -Betas.Beta FROM Variant INNER JOIN Betas 
  ON Variant.chromosome = printf('%02d', Betas.chr_name)  
  AND Variant.position = Betas.chr_position 
  AND Variant.allele1 = Betas.effect_allele 
  AND Variant.allele2 = Betas.noneffect_allele;"

# Filter BGEN to alleles of interest
bgenix -g "${BGENTMP}.bgen" -table Joined > "${FINAL_BGEN}.bgen"
bgenix -g "${FINAL_BGEN}.bgen" -index

# Convert BGEN to PLINK2 PGEN format
plink2 --bgen "${FINAL_BGEN}.bgen" ref-first \
  --hard-call-threshold 0.1 \
  --sample "$SAMPLE_FILE" \
  --memory 15000 \
  --set-all-var-ids @:#_\$r_\$a \
  --freq \
  --make-pgen \
  --out "${RAW_OUT}"


# Filter strand-ambiguous SNPs
awk '/^[^#]/ { if( $5>0.4 && $5<0.6 && ( ($3=="A" && $4=="T") || ($4=="T" && $3=="A") || ($3=="C" && $4=="G") || ($4=="G" && $3=="C") ) ) { print $0 }}' \
"${RAW_OUT}.afreq" > "${OUT_DIR}/exclrsIDs_ambiguous.txt"

# Prepare unified MFI table
for chr in {1..22}; do
  awk -v chr=$chr 'BEGIN {FS="\t"; OFS="\t"} { print chr,$0,chr":"$3"_"$4"_"$5 }' \
    "${MFI_PREFIX}${chr}_v3.txt"
done > "${OUT_DIR}/ukb_mfi_all_v3.tsv"

# SNP QC
plink2 --pfile "${RAW_OUT}" \
  --memory 15000 \
  --exclude "${OUT_DIR}/exclrsIDs_ambiguous.txt" \
  --extract-col-cond "${OUT_DIR}/ukb_mfi_all_v3.tsv" 9 10 --extract-col-cond-min 0.4 \
  --maf 0.005 \
  --write-snplist \
  --make-pgen \
  --out "${SNPQC_OUT}"

# Sample QC
plink2 --pfile "${RAW_OUT}" \
  --memory 15000 \
  --extract "${SNPQC_OUT}.snplist" \
  --keep-fam "$KEEPFAM" \
  --write-samples \
  --out "${SAMPLEQC_OUT}"

# Create score file from matched betas
sqlite3 -separator " " -list "${BGENTMP}.bgen.bgi" "
SELECT chr_name || ':' || position || '_' || allele1 || '_' || allele2, allele2, Beta FROM Joined;" \
> "${SCORE_FILE}"

# Calculate PRS
plink2 --pfile "${RAW_OUT}" \
  --memory 15000 \
  --extract "${SNPQC_OUT}.snplist" \
  --keep "${SAMPLEQC_OUT}.id" \
  --score "${SCORE_FILE}" list-variants ignore-dup-ids no-mean-imputation cols=+scoresums \
  --out "${PRS_OUT}"
