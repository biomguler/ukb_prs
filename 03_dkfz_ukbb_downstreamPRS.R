#############################################################################
# Title:    PRS Evaluation and Visualization (Script Version)
# Author:   Murat Guler
# Contact:  murat.guler@dkfz.de
# Usage:    Rscript 03_dkfz_ukbb_downstreamPRS.R <PRS.sscore> <phenotype.txt> <phenotype_colname> <out_dir>
#############################################################################
# Arguments
#############################################################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript 03_dkfz_ukbb_downstreamPRS.R <PRS.sscore> <phenotype.txt> <phenotype_colname> <out_dir>")
}

prs_file <- args[1]
pheno_file <- args[2]
phenocol_name <- args[3]
out_dir <- args[4]
#############################################################################
# Required packages (will automatically install if missing)
#############################################################################
required_packages <- c("data.table", "tidyverse", "pROC", "broom")

# Function to check, install, and load packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Package '%s' not found. Attempting to install...", pkg))
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
        message(sprintf("Successfully installed '%s'.", pkg))
      }, error = function(e) {
        message(sprintf("Laura and Andronika! Failed to install package '%s': %s", pkg, e$message))
        next
      })
    }
    
    # Try to load the package
    success <- tryCatch({
      library(pkg, character.only = TRUE)
      TRUE
    }, error = function(e) {
      message(sprintf("Laura and Andronika! Failed to load package '%s': %s", pkg, e$message))
      FALSE
    })
    
    if (!success) {
      message(sprintf("Please check the package '%s' manually.", pkg))
    }
  }
}


# Execute the function
install_if_missing(required_packages)

#############################################################################
# Load Input Data
#############################################################################
# Load PRS
prs <- fread(prs_file) %>%
  rename_with(~ gsub("^#", "", .x)) %>%  # remove '#' from header
  select(IID, NAMED_ALLELE_DOSAGE_SUM, SCORE1_SUM) %>%
  rename(IID = IID, PRS_unweighted = NAMED_ALLELE_DOSAGE_SUM, PRS_weighted = SCORE1_SUM)

# Load phenotype
pheno <- fread(pheno_file)

# Check and rename the selected phenotype column
if (!(phenocol_name %in% colnames(pheno))) {
  stop(paste("Phenotype column", phenocol_name, "not found in phenotype file."))
}

pheno <- pheno %>%
  rename(phenotype_col = !!sym(phenocol_name))

# Merge
df <- inner_join(prs, pheno, by = "IID")

#############################################################################
# Plot PRS Distribution: Case vs Control
#############################################################################
df <- df %>%
  mutate(phenotype = factor(phenotype_col, levels = c(0, 1), labels = c("Control", "Case")))

pdf(file.path(out_dir,paste0(phenocol_name,"_","weighted_prs_distribution.pdf")))
ggplot(df, aes(x = PRS_weighted, fill = phenotype)) +
  geom_density(alpha = 0.8) +
  labs(title = "PRS-weighted Distribution: Cases vs Controls", x = "PRS", fill = "Phenotype") +
  theme_minimal()
dev.off()

pdf(file.path(out_dir,paste0(phenocol_name,"_","unweighted_prs_distribution.pdf")))
ggplot(df, aes(x = PRS_unweighted, fill = phenotype)) +
  geom_density(alpha = 0.8) +
  labs(title = "PRS-unweighted Distribution: Cases vs Controls", x = "PRS", fill = "Phenotype") +
  theme_minimal()
dev.off()

#############################################################################
# Boxplot
#############################################################################
pdf(file.path(out_dir,paste0(phenocol_name,"_","weighted_prs_boxplot.pdf")))
ggplot(df, aes(x = phenotype, y = PRS_weighted, fill = phenotype)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "PRS_weighted by Case-Control Status", x = "Phenotype", y = "PRS") +
  theme_minimal()
dev.off()

pdf(file.path(out_dir,paste0(phenocol_name,"_","unweighted_prs_boxplot.pdf")))
ggplot(df, aes(x = phenotype, y = PRS_unweighted, fill = phenotype)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "PRS_unweighted by Case-Control Status", x = "Phenotype", y = "PRS") +
  theme_minimal()
dev.off()

#############################################################################
# Wilcoxon test
#############################################################################
sink(file.path(out_dir,paste0(phenocol_name,"_","weighted_prs_wilcoxon_test.txt")))
print(wilcox.test(PRS_weighted ~ phenotype, data = df))
sink()

sink(file.path(out_dir,paste0(phenocol_name,"_","unweighted_prs_wilcoxon_test.txt")))
print(wilcox.test(PRS_unweighted ~ phenotype, data = df))
sink()


#############################################################################
# AUC Comparison: Continuous PRS vs Covariates
# --------------------------------------------
# Evaluate and compare model performance with AUC:
# - PRS only
# - Age + Sex only
# - Age + Sex + PRS
#############################################################################

# Weighted PRS
model_prs_only_weighted <- glm(phenotype ~ PRS_weighted, data = df, family = binomial())
model_base_weighted <- glm(phenotype ~ Age + Sex, data = df, family = binomial())
model_full_weighted <- glm(phenotype ~ PRS_weighted + Age + Sex, data = df, family = binomial())

roc_prs_only_weighted <- roc(df$phenotype, fitted(model_prs_only_weighted))
roc_base_weighted <- roc(df$phenotype, fitted(model_base_weighted))
roc_full_weighted <- roc(df$phenotype, fitted(model_full_weighted))

auc_comparison_weighted <- data.frame(
  Model = c("PRS_only", "Age+Sex", "Age+Sex+PRS"),
  AUC = c(auc(roc_prs_only_weighted), auc(roc_base_weighted), auc(roc_full_weighted))
)

fwrite(auc_comparison_weighted, file.path(out_dir, paste0(phenocol_name, "_weighted_auc_model_comparison.tsv")), sep = "\t")

pdf(file.path(out_dir, paste0(phenocol_name, "_weighted_auc_model_comparison_plot.pdf")))
plot(roc_prs_only_weighted, col = "blue", main = "AUC Comparison: Weighted PRS", lwd = 2)
plot(roc_base_weighted, col = "green", add = TRUE, lwd = 2)
plot(roc_full_weighted, col = "red", add = TRUE, lwd = 2)
legend("bottomright", legend = c("PRS only", "Age + Sex", "Age + Sex + PRS"),
       col = c("blue", "green", "red"), lwd = 2)
dev.off()

# Unweighted PRS
model_prs_only_unweighted <- glm(phenotype ~ PRS_unweighted, data = df, family = binomial())
model_base_unweighted <- glm(phenotype ~ Age + Sex, data = df, family = binomial())
model_full_unweighted <- glm(phenotype ~ PRS_unweighted + Age + Sex, data = df, family = binomial())

roc_prs_only_unweighted <- roc(df$phenotype, fitted(model_prs_only_unweighted))
roc_base_unweighted <- roc(df$phenotype, fitted(model_base_unweighted))
roc_full_unweighted <- roc(df$phenotype, fitted(model_full_unweighted))

auc_comparison_unweighted <- data.frame(
  Model = c("PRS_only", "Age+Sex", "Age+Sex+PRS"),
  AUC = c(auc(roc_prs_only_unweighted), auc(roc_base_unweighted), auc(roc_full_unweighted))
)

fwrite(auc_comparison_unweighted, file.path(out_dir, paste0(phenocol_name, "_unweighted_auc_model_comparison.tsv")), sep = "\t")

pdf(file.path(out_dir, paste0(phenocol_name, "_unweighted_auc_model_comparison_plot.pdf")))
plot(roc_prs_only_unweighted, col = "blue", main = "AUC Comparison: Unweighted PRS", lwd = 2)
plot(roc_base_unweighted, col = "green", add = TRUE, lwd = 2)
plot(roc_full_unweighted, col = "red", add = TRUE, lwd = 2)
legend("bottomright", legend = c("PRS only", "Age + Sex", "Age + Sex + PRS"),
       col = c("blue", "green", "red"), lwd = 2)
dev.off()

#############################################################################
# AUC Comparison (DeLong Test) - Weighted PRS
#############################################################################
delong_weighted_prs_vs_full <- roc.test(roc_prs_only_weighted, roc_full_weighted, method = "delong")
delong_weighted_base_vs_full <- roc.test(roc_base_weighted, roc_full_weighted, method = "delong")

sink(file.path(out_dir, paste0(phenocol_name, "_weighted_auc_delong_tests.txt")))
cat("DeLong Test: PRS only vs Age+Sex+PRS\n")
print(delong_weighted_prs_vs_full)
cat("\nDeLong Test: Age+Sex vs Age+Sex+PRS\n")
print(delong_weighted_base_vs_full)
sink()

#############################################################################
# AUC Comparison (DeLong Test) - Unweighted PRS
#############################################################################
delong_unweighted_prs_vs_full <- roc.test(roc_prs_only_unweighted, roc_full_unweighted, method = "delong")
delong_unweighted_base_vs_full <- roc.test(roc_base_unweighted, roc_full_unweighted, method = "delong")

sink(file.path(out_dir, paste0(phenocol_name, "_unweighted_auc_delong_tests.txt")))
cat("DeLong Test: PRS only vs Age+Sex+PRS\n")
print(delong_unweighted_prs_vs_full)
cat("\nDeLong Test: Age+Sex vs Age+Sex+PRS\n")
print(delong_unweighted_base_vs_full)
sink()

#############################################################################
# Logistic Regression with Covariates
#############################################################################
sink(file.path(out_dir,paste0(phenocol_name,"_","weighted_logistic_regression.txt")))
# Create formula dynamically
pc_terms <- paste0("PC", 1:10)
formula_str <- as.formula(paste("phenotype ~ PRS_weighted + Age + Sex +", paste(pc_terms, collapse = " + ")))

# Fit the model
glm_bin <- glm(formula_str, data = df, family = binomial())

print(summary(glm_bin))
sink()

sink(file.path(out_dir,paste0(phenocol_name,"_","unweighted_logistic_regression.txt")))
pc_terms <- paste0("PC", 1:10)
formula_str <- as.formula(paste("phenotype ~ PRS_unweighted + Age + Sex +", paste(pc_terms, collapse = " + ")))

# Fit the model
glm_bin <- glm(formula_str, data = df, family = binomial())

print(summary(glm_bin))
sink()
#############################################################################
# Quantile-based ORs
#############################################################################
# Weighted
df <- df %>%
  mutate(quantile_weighted = factor(ntile(PRS_weighted, 5), levels = 1:5))

df$quantile_weighted <- relevel(df$quantile_weighted, ref = "1")

logit_quantile_weighted <- glm(
  phenotype ~ quantile_weighted + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = df,
  family = binomial()
)

or_table_weighted <- tidy(logit_quantile_weighted, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(str_detect(term, "quantile_weighted")) %>%
  mutate(Quantile = str_replace(term, "quantile_weighted", "Quantile ")) %>%
  rename(
    OR = estimate,
    OR_lower_CI = conf.low,
    OR_upper_CI = conf.high
  ) %>%
  select(Quantile, OR, OR_lower_CI, OR_upper_CI)

case_control_counts_weighted <- df %>%
  group_by(quantile_weighted) %>%
  summarise(
    Cases = sum(phenotype == "Case"),
    Controls = sum(phenotype == "Control"),
    .groups = "drop"
  ) %>%
  mutate(Quantile = paste0("Quantile ", quantile_weighted)) %>%
  select(Quantile, Cases, Controls)

or_table_weighted_final <- full_join(or_table_weighted, case_control_counts_weighted, by = "Quantile", )


fwrite(or_table_weighted_final, file.path(out_dir,paste0(phenocol_name,"_","weighted_odds_ratios_by_quantile.tsv")), sep = "\t")

pdf(file.path(out_dir,paste0(phenocol_name,"_","weighted_quantile_or_forestplot.pdf")), width = 7, height = 5)
ggplot(or_table_weighted_final, aes(x = Quantile, y = OR, ymin = OR_lower_CI, ymax = OR_upper_CI)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = "Odds Ratios by PRS Quantile (vs 1st)", x = "PRS Quantile", y = "OR (95% CI)") +
  theme_minimal() 
dev.off()

# Unweighted
df <- df %>%
  mutate(quantile_unweighted = factor(ntile(PRS_unweighted, 5), levels = 1:5))

df$quantile_unweighted <- relevel(df$quantile_unweighted, ref = "1")

logit_quantile_unweighted <- glm(
  phenotype ~ quantile_unweighted + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
  data = df,
  family = binomial()
)

or_table_unweighted <- tidy(logit_quantile_unweighted, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(str_detect(term, "quantile_unweighted")) %>%
  mutate(Quantile = str_replace(term, "quantile_unweighted", "Quantile ")) %>%
  rename(
    OR = estimate,
    OR_lower_CI = conf.low,
    OR_upper_CI = conf.high
  ) %>%
  select(Quantile, OR, OR_lower_CI, OR_upper_CI)

case_control_counts_unweighted <- df %>%
  group_by(quantile_unweighted) %>%
  summarise(
    Cases = sum(phenotype == "Case"),
    Controls = sum(phenotype == "Control"),
    .groups = "drop"
  ) %>%
  mutate(Quantile = paste0("Quantile ", quantile_unweighted)) %>%
  select(Quantile, Cases, Controls)

or_table_unweighted_final <- full_join(or_table_unweighted, case_control_counts_unweighted, by = "Quantile")

fwrite(or_table_unweighted_final, file.path(out_dir, paste0(phenocol_name, "_", "unweighted_odds_ratios_by_quantile.tsv")), sep = "\t")

pdf(file.path(out_dir, paste0(phenocol_name, "_", "unweighted_quantile_or_forestplot.pdf")), width = 7, height = 5)
ggplot(or_table_unweighted_final, aes(x = Quantile, y = OR, ymin = OR_lower_CI, ymax = OR_upper_CI)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = "Odds Ratios by PRS Quantile (Unweighted vs 1st)", x = "PRS Quantile", y = "OR (95% CI)") +
  theme_minimal()
dev.off()

#End