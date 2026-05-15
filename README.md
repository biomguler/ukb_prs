# DKFZ UKBB-PRS: Polygenic Risk Score Pipeline

**Author:** Murat Güler  
**Contact:** murat.guler@dkfz.de  
**Version:** 2.0 (Development with Enhanced QC)

---

## Overview

This is a comprehensive pipeline for calculating and analyzing polygenic risk scores (PRS) in the UK Biobank (UKBB) cohort. The pipeline performs:

1. **PGS Catalog Integration** — Download and format PGS Catalog scoring files
2. **Core PRS Calculation** — Compute PRS from BGEN genotype data  
3. **Batch Processing** — Automate trait-level PRS calculations
4. **Downstream Analysis** — Perform statistical evaluation and visualization
5. **QC Export** — Extract hardcall effect-allele matrices for validation

The pipeline was designed for the DKFZ ODCF cluster environment and follows best practices for PRS computation with extensive quality control.

---

## Pipeline Architecture

```
00_pgs2score.sh
    ↓
    Downloads PGS from Harmonized Catalog
    Outputs: scores/PGS<ID>/<PGS_ID>_beta.csv
    ↓
02_dkfz_ukbb_prs_runner.sh
    ↓
    Iterates over all *_beta.csv files
    ↓
01_dkfz_ukbb_prs_core.sh (main calculation)
    ├─ BGEN extraction & matching
    ├─ SNP/Sample QC
    ├─ Allele orientation checking
    ├─ PRS scoring
    └─ Outputs: results/<TRAIT>/*
    ↓
03_dkfz_ukbb_downstreamPRS.R
    ├─ Distribution analysis
    ├─ Quantile stratification
    ├─ Association testing
    ├─ ROC/AUC evaluation
    └─ Outputs: figures + statistics
    ↓
04_export_hardcall_effect_prs_for_federico.sh (QC export)
    └─ Outputs: hardcall effect-allele matrix for validation
```

---

## Script Details

### **00_pgs2score.sh** — Download PGS Catalog Scoring File

**Purpose:** Download and harmonize a polygenic risk score from the EBI Harmonized PGS Catalog.

**Usage:**
```bash
bash 00_pgs2score.sh <PGS_ID> [BUILD] [OUT_DIR]
```

**Arguments:**
- `<PGS_ID>` (required): PGS Catalog identifier (e.g., `PGS000083`)
- `[BUILD]` (optional): Genome build (`GRCh37` or `GRCh38`); defaults to `GRCh37`
- `[OUT_DIR]` (optional): Output directory path; defaults to current directory

**Examples:**
```bash
# Single PGS download
bash 00_pgs2score.sh PGS000083 GRCh37 scores

# Multiple PGS downloads
for PGS_ID in PGS000083 PGS000018 PGS000021; do
  bash 00_pgs2score.sh "$PGS_ID" GRCh37 scores
done

# From file
while read -r PGS_ID; do
  [ -z "$PGS_ID" ] && continue
  bash 00_pgs2score.sh "$PGS_ID" GRCh37 scores
done < pgs_ids.txt
```

**Requirements:**
- `wget` or `curl` — for downloading files
- `gzip` — for decompressing
- `awk` — for text processing

**Input:**
- Harmonized PGS file from EBI FTP (automatically downloaded)
- Format: tab-separated with columns including `rsID`, `hm_chr`, `hm_pos`, `effect_allele`, `other_allele`, `effect_weight`, `allelefrequency_effect`

**Output:**
```
scores/PGS000083/PGS000083_beta.csv
```

**Output Format:**
```csv
rsid,chr_name,chr_position,effect_allele,noneffect_allele,Beta,eaf,chr_pos
rs1234567,1,12345678,A,G,0.123,0.45,1:12345678
rs7654321,2,87654321,T,C,-0.089,0.32,2:87654321
```

**Quality Control:**
- Validates gzip format of downloaded file
- Skips variants with missing required fields (rsID, position, alleles, beta)
- Removes accidental carriage returns (CR) from Windows-formatted files
- Uses harmonized positions (hm_chr, hm_pos, hm_rsID) when available
- Reports number of variants written vs. skipped

**Notes:**
- Requires EBI Harmonized PGS Catalog access
- Supports both GRCh37 (hg19) and GRCh38 (hg38) builds
- Pre-harmonized variants avoid allele mismatch issues

---

### **01_dkfz_ukbb_prs_core.sh** — Core PRS Calculation

**Purpose:** Core PRS calculation pipeline. Processes UKBB BGEN genotypes, matches to beta file, performs QC, and computes weighted and unweighted PRS scores.

**Usage:**
```bash
BETAS=<betas.csv> OUT_PREFIX=<name> OUT_DIR=<outdir> bash 01_dkfz_ukbb_prs_core.sh
```

**Environment Variables:**
- `BETAS` (required): Path to beta file (`*_beta.csv`)
- `OUT_PREFIX` (required): Output prefix for result files
- `OUT_DIR` (required): Output directory

**Requirements:**
- `sqlite3` — for variant matching database
- `plink2` — for genotype processing
- `bgenix` + `cat-bgen` — for BGEN file handling
- `awk`, `sort`, `wc` — standard utilities
- Micromamba/conda — for bgenix environment setup

**Input Files:**
```
trait_beta.csv                    # Beta file
ukb22828_c[1-22]_b0_v3.bgen      # BGEN genotypes (chr 1-22)
ukb_mfi_chr[1-22]_v3.txt          # MFI quality metrics
ukb22828_c1_b0_v3_s487188.sample # Sample file
usedinpca.txt                     # QC-filtered sample list
```

**Beta File Format:**
```csv
rsid,chr_name,chr_position,effect_allele,noneffect_allele,Beta,eaf,chr_pos
rs6511720,19,11202306,T,G,-0.211427,0.1075,19:11202306
rs4420638,19,45422946,G,A,0.16801,0.1797,19:45422946
```

**Output Files:**

| File | Description |
|------|-------------|
| `<PREFIX>_raw_pvar_ids.txt` | Variant IDs before QC |
| `<PREFIX>_snpqc_ids.txt` | Variant IDs after SNP QC |
| `<PREFIX>_matched_variants.tsv` | Matched variants with allele info |
| `<PREFIX>_unmatched_betas.tsv` | Variants in beta file not in BGEN |
| `<PREFIX>_orientation_counts.tsv` | Allele orientation summary |
| `<PREFIX>_score.txt` | PLINK score file for PRS calculation |
| `<PREFIX>_effect_alleles_for_export.txt` | Effect allele reference |
| `<PREFIX>_PRS.sscore` | **Final PRS scores** |

**Processing Steps:**

1. **Prepare SNP Lists** — Extract rsID and chr:pos from beta file
2. **Extract BGEN Candidates** — Use `bgenix` with position ranges
3. **Combine BGEN Chunks** — Merge chr 1-22 into single BGEN file
4. **Build Matching Database** — Load beta file into SQLite BGEN index
5. **Match & Orient** — Join BGEN variants to betas, check allele orientation
6. **Filter Ambiguous SNPs** — Exclude A-T and C-G SNPs near MAF 0.5
7. **SNP QC** — Apply MFI INFO ≥ 0.4 filter and MAF ≥ 0.005
8. **Sample QC** — Keep only samples from PCA-filtered list
9. **Prepare Score File** — Create PLINK2 score format with correct effect allele
10. **Calculate PRS** — Compute weighted and unweighted PRS per individual

**Quality Control Features (Development Version):**
- **Allele Orientation Check:** Ensures effect allele correctly identified
- **Ambiguous SNP Filtering:** Removes strand-ambiguous variants near MAF 0.5
- **INFO Score Filter:** MFI INFO ≥ 0.4 (UKBB imputation quality)
- **MAF Filter:** Minor allele frequency ≥ 0.5% to reduce noise
- **Sample Filtering:** Restricts to PCA-QC individuals
- **Variant Matching Audit:** Reports orientation distribution and unmatched variants
- **Missing Data Handling:** No mean-imputation; missing genotypes treated as 0

---

### **02_dkfz_ukbb_prs_runner.sh** — Batch Processing Runner

**Purpose:** Automate PRS calculation across multiple traits by iterating over beta files and calling the core script.

**Usage:**
```bash
bash 02_dkfz_ukbb_prs_runner.sh
```

**How It Works:**
1. Finds all `*_beta.csv` files in the current directory
2. For each file, extracts trait name (e.g., `trait1_beta.csv` → `trait1`)
3. Calls `01_dkfz_ukbb_prs_core.sh` with appropriate environment variables
4. Captures stdout/stderr in separate log files

**Output Directory Structure:**
```
results/
├── trait1/
│   ├── trait1_raw_pvar_ids.txt
│   ├── trait1_snpqc_ids.txt
│   ├── trait1_PRS.sscore
│   └── [other intermediates]
├── trait2/
│   └── [similar files]
└── logs/
    ├── trait1.prs.log
    ├── trait1.prs.err
    ├── trait2.prs.log
    └── trait2.prs.err
```

**Requirements:**
- `01_dkfz_ukbb_prs_core.sh` must be in the same directory or in PATH
- Beta files must be in current directory (or modify script)

**Configuration:**
If beta files are in a different directory or runner script is elsewhere, modify:
```bash
# Edit this line to specify exact path
beta_files=(../path/to/*_beta.csv)
```

---

### **03_dkfz_ukbb_downstreamPRS.R** — Downstream Analysis & Visualization

**Purpose:** Comprehensive downstream analysis of PRS including distribution analysis, stratification into quantiles, association testing, and ROC/AUC evaluation.

**Usage:**
```bash
Rscript 03_dkfz_ukbb_downstreamPRS.R \
  <PRS.sscore> <phenotype.txt> <phenotype_colname> <out_dir> \
  <n_quantile_groups> <quantile_basis>
```

**Arguments:**
- `<PRS.sscore>` — PLINK2 .sscore file from core script
- `<phenotype.txt>` — Tab/space-separated phenotype file with IID and phenotype columns
- `<phenotype_colname>` — Column name for phenotype in the phenotype file
- `<out_dir>` — Output directory for figures and statistics
- `<n_quantile_groups>` — Number of quantile groups (e.g., 5 for quintiles, 10 for deciles)
- `<quantile_basis>` — Basis for quantile calculation: `controls`, `all`, `combined`, `cases+controls`

**Examples:**
```bash
# 5 quantiles based on controls
Rscript 03_dkfz_ukbb_downstreamPRS.R \
  results/trait1/trait1_PRS.sscore \
  phenotypes/trait1.tsv \
  PDAC \
  results/downstream/trait1 \
  5 \
  controls

# 10 deciles based on all samples
Rscript 03_dkfz_ukbb_downstreamPRS.R \
  results/trait1/trait1_PRS.sscore \
  phenotypes/trait1.tsv \
  PDAC \
  results/downstream/trait1 \
  10 \
  all
```

**Phenotype File Format:**
```
IID           PDAC  Age  Sex  PC1    PC2    PC3    ...PC10
0000001       0     50   1    0.001  0.002  0.001  ...
0000002       1     55   2   -0.001  0.003  0.002  ...
0000003       0     48   1    0.002  0.001  0.003  ...
```

**Requirements:**
- R packages: `data.table`, `tidyverse`, `pROC`, `broom`
- Package installation is automatic if missing

**Input Requirements:**
- PRS file must contain columns: `IID`, `WEIGHTED_SUM`, `UNWEIGHTED_SUM`
- Phenotype file must contain: `IID`, phenotype column, `Age`, `Sex`, `PC1`–`PC10`
- Phenotype values must be binary (0/1)

**Output Files:**

| File | Description |
|------|-------------|
| `<phenocol>_analysis_summary.txt` | Analysis metadata and sample counts |
| `<phenocol>_weighted_prs_distribution.pdf` | Weighted PRS density by case/control |
| `<phenocol>_unweighted_prs_distribution.pdf` | Unweighted PRS density by case/control |
| `<phenocol>_weighted_prs_boxplot.pdf` | Weighted PRS boxplot |
| `<phenocol>_unweighted_prs_boxplot.pdf` | Unweighted PRS boxplot |
| `<phenocol>_weighted_prs_wilcoxon_test.txt` | Wilcoxon test results (weighted) |
| `<phenocol>_unweighted_prs_wilcoxon_test.txt` | Wilcoxon test results (unweighted) |
| `<phenocol>_weighted_auc_model_comparison.tsv` | AUC for PRS-only, Age+Sex, Age+Sex+PRS |
| `<phenocol>_weighted_auc_model_comparison_plot.pdf` | ROC curve comparison |
| `<phenocol>_<n>groups_<basis>_quantile_or_summary.tsv` | Odds ratios by PRS quantile |
| `<phenocol>_<n>groups_<basis>_quantile_counts.tsv` | Sample counts per quantile |
| `<phenocol>_<n>groups_<basis>_quantile_logit_model.txt` | Full quantile model output |

**Key Analyses:**

1. **Distributions** — Density and boxplots of weighted/unweighted PRS by phenotype
2. **Wilcoxon Test** — Non-parametric test for PRS difference between cases/controls
3. **Quantile Stratification** — Divide samples into N quantiles based on PRS
   - Reference basis can be controls only or all samples
   - Handles tied values intelligently
4. **Odds Ratios** — Logistic regression for each quantile vs. reference
5. **Model Comparison** — ROC curves and AUC for:
   - PRS only
   - Age + Sex
   - Age + Sex + PRS

**Quantile Basis Options:**
- `controls` — Quantiles calculated from control PRS distribution
- `all` / `combined` / `cases+controls` — Quantiles from all samples

---

### **04_export_hardcall_effect_prs_for_federico.sh** — QC Hardcall Export

**Purpose:** Extract hardcall-only effect-allele genotypes and calculate unweighted/weighted hardcall PRS for validation and QC export (development feature).

**Usage:**
```bash
bash 04_export_hardcall_effect_prs_for_federico.sh \
  <SNPQC_PFILE_PREFIX> <SCORE_FILE> <OUT_PREFIX> [KEEP_FILE]
```

**Arguments:**
- `<SNPQC_PFILE_PREFIX>` — Path to QC-filtered PGEN files (without extension)
- `<SCORE_FILE>` — PLINK score file (whitespace-separated)
- `<OUT_PREFIX>` — Output file prefix
- `[KEEP_FILE]` (optional) — Sample IDs to keep (if subset desired)

**Example:**
```bash
bash 04_export_hardcall_effect_prs_for_federico.sh \
  results/trait1/trait1_snpQC \
  results/trait1/trait1_score.txt \
  results/exports/trait1_hardcall \
  pca_samples.txt
```

**Input Files:**
- `.pgen`, `.pvar`, `.psam` — PLINK2 binary genotype files
- Score file with columns: `ID`, `ALLELE`, `WEIGHTED`, `UNWEIGHTED`

**Output Files:**

| File | Description |
|------|-------------|
| `<PREFIX>_tmp_hardcalls.pgen` | Temporary: hard-call genotypes only |
| `<PREFIX>_effect_allele_matrix.raw` | Effect-allele count matrix (0/1/2) |
| `<PREFIX>.tsv` | **Final output: PRS scores** |

**Output TSV Format:**
```
IID    hardcall_effect_allele_sum    hardcall_effect_allele_sumxbeta
0000001    45    3.456
0000002    48    4.123
0000003    42    2.987
```

**Processing Steps:**

1. **Filter Score to PGEN Variants** — Keep only variants present in PGEN
2. **Create Effect-Allele Reference** — Map variant ID to effect allele
3. **Erase Dosage** — Convert to hard calls only (0/1/2)
4. **Export Dosage Matrix** — Extract effect-allele counts per individual
5. **Calculate PRS** — Compute both unweighted (allele count) and weighted (allele count × beta) sums

**Quality Control:**
- Hard-call conversion only (no dosages)
- Missing genotypes treated as 0 (no mean-imputation)
- Sanity checks for non-0/1/2 values after erase-dosage
- Per-individual reporting of variant counts

---

## Workflow: Complete Example

### 1. Obtain Beta File from PGS Catalog

```bash
bash 00_pgs2score.sh PGS000083 GRCh37 scores
# Output: scores/PGS000083/PGS000083_beta.csv
```

### 2. Run Core PRS Pipeline

Copy the beta file to the working directory:
```bash
cp scores/PGS000083/PGS000083_beta.csv ./

# Run via runner script
bash 02_dkfz_ukbb_prs_runner.sh

# Or run manually
BETAS=PGS000083_beta.csv OUT_PREFIX=PGS000083 OUT_DIR=results bash 01_dkfz_ukbb_prs_core.sh
```

### 3. Downstream Analysis

Prepare phenotype file, then run:
```bash
Rscript 03_dkfz_ukbb_downstreamPRS.R \
  results/PGS000083/PGS000083_PRS.sscore \
  phenotypes.tsv \
  PDAC \
  results/downstream \
  5 \
  controls
```

### 4. (Optional) Export for Validation

```bash
bash 04_export_hardcall_effect_prs_for_federico.sh \
  results/PGS000083/PGS000083_snpQC \
  results/PGS000083/PGS000083_score.txt \
  results/exports/PGS000083_hardcall
```

---

## Software Dependencies

### Cluster Modules (DKFZ ODCF)

```bash
module load SQLite/3.46.0-GCCcore-14.1.0
module load PLINK/2.00a6_amd_avx2
module load Micromamba/2.0.2-0
```

### System Tools

| Tool | Purpose | Install |
|------|---------|---------|
| `wget` / `curl` | Download files | `apt install wget curl` |
| `gzip` | Decompress | `apt install gzip` |
| `awk` | Text processing | `apt install gawk` |
| `sqlite3` | Database queries | Loaded via module |
| `bgenix` | BGEN filtering | Micromamba env |
| `cat-bgen` | BGEN merging | Micromamba env |
| `plink2` | Genotype processing | Loaded via module |

### R Packages

Auto-installed if missing:
- `data.table` — Fast data I/O
- `tidyverse` — Data manipulation
- `pROC` — ROC curve analysis
- `broom` — Model summaries

---

## Quality Control Summary (Development Features)

This pipeline includes comprehensive QC at multiple stages:

### **Variant-Level QC**
- ✅ Allele orientation checking (BGEN vs. beta file)
- ✅ Strand-ambiguous SNP filtering (A-T, C-G near MAF 0.5)
- ✅ Imputation quality (MFI INFO ≥ 0.4)
- ✅ MAF filtering (≥ 0.5%)
- ✅ Variant matching audit (reports unmatched variants)

### **Sample-Level QC**
- ✅ PCA-filtered sample restriction
- ✅ Complete-case handling in downstream analysis
- ✅ Per-individual variant count validation

### **Data Integrity**
- ✅ Gzip validation (script 00)
- ✅ Missing field detection (script 00)
- ✅ Windows CR character removal (script 00)
- ✅ SQLite join validation (script 01)
- ✅ Hard-call sanity checks (script 04)

### **Audit Output**
- `*_orientation_counts.tsv` — Allele orientation breakdown
- `*_matched_variants.tsv` — Full matched variant detail
- `*_unmatched_betas.tsv` — Variants in beta not in BGEN
- `*_effect_alleles_for_export.txt` — Effect allele reference

---

## Troubleshooting

### No variants matched between BGEN and beta file
- Check chromosome naming (1 vs chr1)
- Verify genome build matches (GRCh37 vs GRCh38)
- Confirm beta file format (CSV with required columns)

### Module load errors
- Verify DKFZ cluster environment
- Check available modules: `module avail`

### BGEN/MFI file not found
- Verify `UKBB_BASE` path in script
- Check UKB data accessibility

### Phenotype file issues
- Ensure IID column present
- Verify phenotype column name matches input
- Check for binary (0/1) encoding of phenotype

---

## Citation

If using this pipeline, please cite:

Güler, M., et al. DKFZ UKBB-PRS: A Pipeline for Polygenic Risk Score Calculation in UK Biobank. Version 2.0. 2025.

Reference pipeline: [Collister et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35251129/)

---

## License

See [LICENSE](LICENSE) file for details.

---

## Contact

For questions or issues, contact:
- **Author:** Murat Güler  
- **Email:** murat.guler@dkfz.de
