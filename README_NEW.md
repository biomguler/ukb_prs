# DKFZ UKBB-PRS: A Pipeline for Beginners

**Author:** Murat Güler  
**Date:** 2025-07-30

## Overview and Motivation

Using UK Biobank cohort genetic data can be challenging for various groups of researchers. One of the most common forms of usage is calculating polygenic risk scores (PRS). This pipeline calculates PRS for a given set of SNPs with minimal effort for the whole UKBB cohort. The pipeline was created for the DKFZ ODCF environment and is adapted from [Collister et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35251129/).

The pipeline consists of two steps:

1. **Step 1 (Unix based):** Calculate PRS for the whole UKBB cohort
2. **Step 2 (R based):** Downstream analysis and performance evaluation of PRS

---

## Step 1: Calculate PRS for Whole UKBB Cohort

This section explains how to use the `01_dkfz_ukbb_prs_core.sh` pipeline to compute polygenic risk scores (PRS) from UK Biobank BGEN data using a list of summary statistics (betas). The script was designed for reproducible use across different phenotypes/traits or beta files, with flexible input/output configuration.

### Requirements

The pipeline uses CentOS, but can be used on Debian by changing modules.

**Full list of tools needed:**
- `bash`
- `awk`
- `sqlite3`
- `bgenix`
- `cat-bgen`
- `plink2`

*Note: Do not worry! You do not need to install any of them—all are available in the cluster.*

### Input Files

| File | Description |
|------|-------------|
| `trait_beta.csv` | Summary statistics with rsID, chromosome, position, alleles, and beta |
| `ukb_imp_chr[1-22]_v3.bgen` | UKB imputed genotype data |
| `ukb_mfi_chr[1-22]_v3.txt` | MFI files with per-SNP quality |
| `ukbA_imp_chrN_v3_sP.sample` | Sample file (usually from UKB) |
| `usedinpca.txt` | List of sample IDs to keep (e.g., from PCA filtering) |

#### Structure of `trait_beta.csv`

**Important point:** The file must be named as `trait_beta.csv`. For example, if you have trait1, trait2, trait3, the file names should be `trait1_beta.csv`, `trait2_beta.csv`, `trait3_beta.csv`, and outputs will be written in `results/trait1`, `results/trait2`, `results/trait3`.

The CSV file should contain summary statistics with the following columns:

| rsid | chr_name | chr_position | effect_allele | noneffect_allele | Beta | eaf | chr_pos |
|------|----------|--------------|---------------|-----------------|------|-----|---------|
| rs6511720 | 19 | 11202306 | T | G | -0.211427 | 0.1075 | 19:11202306 |
| rs4420638 | 19 | 45422946 | G | A | 0.16801 | 0.1797 | 19:45422946 |
| rs629301 | 1 | 109818306 | T | G | 0.157518 | 0.7742 | 1:109818306 |
| rs2328223 | 20 | 17845921 | C | A | 0.14 | 0.2244 | 20:17845921 |

**Column Descriptions:**
- `chr_name`: Chromosome number (1–22)
- `chr_position`: Base-pair position
- `effect_allele` / `noneffect_allele`: Alleles corresponding to the direction of effect (Beta)
- `Beta`: Effect size for the effect allele (used in scoring)
- `eaf`: Effect allele frequency (optional, not used in this pipeline)
- `chr_pos`: Combined column useful for QC (optional)

### Run 01_dkfz_ukbb_prs_core.sh

A runner script `02_dkfz_ukbb_prs_runner.sh` has been created to run the PRS script.

If your beta files are named `trait*_beta.csv` and are in the same path as `01_dkfz_ukbb_prs_core.sh`, you do not need to modify anything.

If the beta files and script are in different paths, modify the runner script to provide exact paths. Note that all beta files must be in the same folder.

```bash
# List all beta files you want to process
for BETAFILE in your/path/*_beta.csv; do
  TRAIT_NAME=$(basename "$BETAFILE" _beta.csv)
  
  echo "Processing trait: $TRAIT_NAME"

  # Call the main PRS script, passing variables dynamically
  BETAS="$BETAFILE"
  OUT_PREFIX="$TRAIT_NAME"
  OUT_DIR="results/${OUT_PREFIX}"

  # You can source the main logic instead of duplicating
  export BETAS OUT_PREFIX OUT_DIR

  bash your/script/path/01_dkfz_ukbb_prs_core.sh
done
```

### How to Run

After editing the runner script `02_dkfz_ukbb_prs_runner.sh`, run it in your terminal:

**In the terminal:**
```bash
bash 02_dkfz_ukbb_prs_runner.sh
```

**As a job:**
```bash
#!/bin/bash
#BSUB -q long
#BSUB -J prs_example
#BSUB -n 2                  
#BSUB -R "rusage[mem=10GB]"
#BSUB -o /path/to/your/logs/prs.out
#BSUB -e /path/to/your/logs/prs.err

cd your/working/dir 

bash 02_dkfz_ukbb_prs_runner.sh
```

*Important note: Please do not run it on submission nodes. Run it as a job or on worker nodes.*

All intermediate and final outputs will be stored in the directory you defined in `OUT_DIR`.

### Output Files

Below is a summary of key output files, assuming your `OUT_PREFIX` is set to `"trait1"`:

| File | Description |
|------|-------------|
| `trait1_raw.*` | PLINK2 files before QC |
| `trait1_snpQC.*` | After SNP quality filtering |
| `trait1_sampleQC.*` | After sample filtering |
| `trait1_PRS.sscore` | Final PRS results (output from `--score`) |
| `trait1_score.txt` | Formatted score file used by PLINK2 |
| `trait1_initial_chr.bgen` | Combined BGEN across all chromosomes |
| `trait1_single_allelic.bgen` | Filtered BGEN (single-allelic, matched to betas) |

All results will be saved under: `results/trait1/`

### Notes

- The script **automatically filters out strand-ambiguous SNPs** using allele frequencies.
- SNPs are further filtered based on:
  - **Imputation INFO score** (from MFI files)
  - **Minor Allele Frequency** (MAF ≥ 0.005)
- Only SNPs that **match the alleles in your `trait_beta.csv`**, including allele order and orientation, are retained.

---

## Step 2: Downstream Analysis and Performance of PRS

After calculating PRS, you need to run downstream analysis. For this pipeline, several analyses are performed:

### Output Files from Downstream Analysis

#### PRS Distribution and Statistical Tests

| Output File | Description |
|-------------|-------------|
| `<phenocol>_weighted_prs_distribution.pdf` | Distribution plot of weighted PRS |
| `<phenocol>_unweighted_prs_distribution.pdf` | Distribution plot of unweighted PRS |
| `<phenocol>_weighted_prs_boxplot.pdf` | Boxplot of weighted PRS by phenotype |
| `<phenocol>_unweighted_prs_boxplot.pdf` | Boxplot of unweighted PRS by phenotype |
| `<phenocol>_weighted_prs_wilcoxon_test.txt` | Wilcoxon test results (weighted) |
| `<phenocol>_unweighted_prs_wilcoxon_test.txt` | Wilcoxon test results (unweighted) |

#### AUC Analysis

| Output File | Description |
|-------------|-------------|
| `<phenocol>_weighted_auc_model_comparison.tsv` | AUC comparison table (weighted) |
| `<phenocol>_unweighted_auc_model_comparison.tsv` | AUC comparison table (unweighted) |
| `<phenocol>_weighted_auc_model_comparison_plot.pdf` | AUC comparison plot (weighted) |
| `<phenocol>_unweighted_auc_model_comparison_plot.pdf` | AUC comparison plot (unweighted) |
| `<phenocol>_weighted_auc_delong_tests.txt` | DeLong test results (weighted) |
| `<phenocol>_unweighted_auc_delong_tests.txt` | DeLong test results (unweighted) |

#### Logistic Regression and Odds Ratios

| Output File | Description |
|-------------|-------------|
| `<phenocol>_weighted_logistic_regression.txt` | Logistic regression results (weighted) |
| `<phenocol>_unweighted_logistic_regression.txt` | Logistic regression results (unweighted) |
| `<phenocol>_weighted_odds_ratios_by_quantile.tsv` | Odds ratios by PRS quantile (weighted) |
| `<phenocol>_unweighted_odds_ratios_by_quantile.tsv` | Odds ratios by PRS quantile (unweighted) |
| `<phenocol>_weighted_quantile_or_forestplot.pdf` | Forest plot of ORs (weighted PRS) |
| `<phenocol>_unweighted_quantile_or_forestplot.pdf` | Forest plot of ORs (unweighted PRS) |

#### Summary of Output Files

| Analysis Type | Weighted Output | Unweighted Output |
|---------------|-----------------|-------------------|
| PRS Distribution | `_weighted_prs_distribution.pdf` | `_unweighted_prs_distribution.pdf` |
| Boxplot | `_weighted_prs_boxplot.pdf` | `_unweighted_prs_boxplot.pdf` |
| Wilcoxon Test | `_weighted_prs_wilcoxon_test.txt` | `_unweighted_prs_wilcoxon_test.txt` |
| AUC Table | `_weighted_auc_model_comparison.tsv` | `_unweighted_auc_model_comparison.tsv` |
| AUC Plot | `_weighted_auc_model_comparison_plot.pdf` | `_unweighted_auc_model_comparison_plot.pdf` |
| DeLong Test | `_weighted_auc_delong_tests.txt` | `_unweighted_auc_delong_tests.txt` |
| Logistic Regression | `_weighted_logistic_regression.txt` | `_unweighted_logistic_regression.txt` |
| Odds Ratio Table | `_weighted_odds_ratios_by_quantile.tsv` | `_unweighted_odds_ratios_by_quantile.tsv` |
| Odds Ratio Plot | `_weighted_quantile_or_forestplot.pdf` | `_unweighted_quantile_or_forestplot.pdf` |

### How to Run for Multiple Traits

You can run the script for multiple traits (e.g., trait1, trait2, trait3) in a loop:

```bash
for dir in ./results/*/; do
  trait=$(basename "$dir")
  prs_file="./results/${trait}/${trait}_PRS.sscore"
  pheno_file="./phenotypes/${trait}_pheno.txt"
  out_dir="./output_${trait}/"

  # Create output directory
  mkdir -p "$out_dir"

  # Run the script
  Rscript 03_dkfz_ukbb_downstreamPRS.R "$prs_file" "$pheno_file" "$trait" "$out_dir"
done
```

This assumes you followed the naming convention: `trait1_beta.csv`, `trait1_pheno.txt`, with phenotype column named `trait1` in the `trait1_pheno.txt` files.

**As a job:**
```bash
#!/bin/bash
#BSUB -q long
#BSUB -J prs_example
#BSUB -n 2                  
#BSUB -R "rusage[mem=10GB]"
#BSUB -o /path/to/your/logs/prs.out
#BSUB -e /path/to/your/logs/prs.err

module load R/4.3.0

cd your/working/dir 

for dir in ./results/*/; do
  trait=$(basename "$dir")
  prs_file="./results/${trait}/${trait}_PRS.sscore"
  pheno_file="./phenotypes/${trait}_pheno.txt"
  out_dir="./output_${trait}/"

  # Create output directory
  mkdir -p "$out_dir"

  # Run the script
  Rscript 03_dkfz_ukbb_downstreamPRS.R "$prs_file" "$pheno_file" "$trait" "$out_dir"
done
```

---

## Real Data Example (for Real Beginners)

*Important note: If you do not have experience with coding, please exactly follow these steps.*

To run the real data example in UKBB, follow the steps below:

### Prepare Files

Create a new folder (optional but highly recommended) and create beta files in the same format as described above.

In this example, we will use 3 phenotypes/traits: CLL, MM, PDAC.

If you run the code below in your terminal, it will create a new folder in your current directory and copy all scripts:

```bash
# Create your folder
mkdir -p prs_example && cd prs_example

# Copy scripts
cp /omics/odcf/analysis/OE0540_projects/ukkb_joint_oe0136/669373/Genomics/Imputation/Imputation_from_genotype/prs_dkfz_test/*_* .

# Add +x to scripts
chmod +x *.sh *.R
```

### Run First Script to Calculate PRS

To calculate PRS for the given beta files, run the one-line code below in your terminal:

```bash
bash 02_dkfz_ukbb_prs_runner.sh > step1.logs 2>&1
```

This code will create a `results` folder under the `prs_example` folder. For each trait, the script will create a subfolder under `results`.

Here's the folder tree of the results:

```
+---prs_example
    |  +---results
    |      |   +---CLL
    |      |   +---MM
    |      |   +---PDAC_TA_34
    |      |   +---PDAC_EUR_34
```

### Final Step: PRS Downstream Analysis

For this step, we need phenotype files for each trait. The code below will copy phenotype files to the `phenotypes` folder in your current working directory:

```bash
# Copy example phenotypes
cp -R /omics/odcf/analysis/OE0540_projects/ukkb_joint_oe0136/669373/Genomics/Imputation/Imputation_from_genotype/prs_dkfz_test/phenotypes .

# For PDAC we need two copies with different names (for each beta file) and need to fix phenotype columns
for trait in PDAC_TA_34 PDAC_EUR_34; do
  cp ./phenotypes/PDAC_pheno.txt ./phenotypes/${trait}_pheno.txt
  awk -v newcol="$trait" 'NR==1{$2=newcol} 1' OFS='\t' ./phenotypes/${trait}_pheno.txt > ./phenotypes/tmp && mv ./phenotypes/tmp ./phenotypes/${trait}_pheno.txt
done
```

After getting phenotype files, run the code below to get all downstream analysis results:

```bash
module load R/4.3.0

for dir in ./results/*/; do
  trait=$(basename "$dir")
  prs_file="./results/${trait}/${trait}_PRS.sscore"
  pheno_file="./phenotypes/${trait}_pheno.txt"
  out_dir="./output_${trait}/"

  # Create output directory
  mkdir -p "$out_dir"

  # Run the script and log output
  echo "===== Running: $trait =====" >> step2.logs
  Rscript 03_dkfz_ukbb_downstreamPRS.R "$prs_file" "$pheno_file" "$trait" "$out_dir" >> step2.logs 2>&1
  echo "===== Done: $trait =====" >> step2.logs
  echo "" >> step2.logs

done
```

After running the pipeline, you will have output folders named `output_<trait>` for each phenotype under the `./prs_example/` directory:

```
prs_example/
├── output_MM/
├── output_CLL/
├── output_PDAC_TA_34/
└── output_PDAC_EUR_34/
```

Each of these contains various result files for weighted and unweighted PRS analysis, including AUC comparison tables, plots, odds ratios by quantile, etc.

### Checking Results

To check all weighted PRS AUC for all models:

```bash
output_summary="all_weighted_auc_results.txt"

# Clear previous file if exists
> "$output_summary"

for file in ./output_*/**_weighted_auc_model_comparison.tsv; do
  echo "===== $(basename "$file" | cut -d'_' -f1) =====" >> "$output_summary"
  cat "$file" >> "$output_summary"
  echo "" >> "$output_summary"
done
```

This code will create an `all_weighted_auc_results.txt` file with content like:

| Phenotype | Model | AUC |
|-----------|-------|-----|
| MM | PRS_only | 0.5150 |
| | Age + Sex | 0.8156 |
| | Age + Sex + PRS | 0.8156 |
| CLL | PRS_only | 0.7090 |
| | Age + Sex | 0.7890 |
| | Age + Sex + PRS | 0.8365 |
| PDAC_EUR_34 | PRS_only | 0.6251 |
| | Age + Sex | 0.8676 |
| | Age + Sex + PRS | 0.8755 |
| PDAC_TA_34 | PRS_only | 0.6251 |
| | Age + Sex | 0.8676 |
| | Age + Sex + PRS | 0.8755 |

The AUC values are extracted from `output_<trait>/<trait>_weighted_auc_model_comparison.tsv`. You can use similar scripts to extract values from other output files.

---

## Additional Notes

- This pipeline is optimized for the DKFZ ODCF environment using CentOS
- All required tools are available on the cluster
- For questions or issues, refer to the original pipeline documentation
- Ensure proper naming conventions for input files to avoid errors
