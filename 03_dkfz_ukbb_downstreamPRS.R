#############################################################################
# Title:    PRS Evaluation and Visualization (Script Version)
# Author:   Murat Guler
# Contact:  murat.guler@dkfz.de
# Usage:
#   Rscript 03_dkfz_ukbb_downstreamPRS.R \
#     <PRS.sscore> <phenotype.txt> <phenotype_colname> <out_dir> \
#     <n_quantile_groups> <quantile_basis>
#
# Example:
#   Rscript 03_dkfz_ukbb_downstreamPRS.R \
#     PGS000083_PRS.sscore phenotype_PDAC.tsv PDAC results/downstream 5 controls
#
# quantile_basis options:
#   controls
#   control
#   all
#   combined
#   cases+controls
#   case_control
#############################################################################

#############################################################################
# Arguments
#############################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(
    paste(
      "Usage:",
      "Rscript 03_dkfz_ukbb_downstreamPRS.R",
      "<PRS.sscore> <phenotype.txt> <phenotype_colname> <out_dir>",
      "<n_quantile_groups> <quantile_basis>",
      "\nExample:",
      "Rscript 03_dkfz_ukbb_downstreamPRS.R PGS000083_PRS.sscore phenotype_PDAC.tsv PDAC results/downstream 5 controls"
    )
  )
}

prs_file <- args[1]
pheno_file <- args[2]
phenocol_name <- args[3]
out_dir <- args[4]
n_quantile_groups <- suppressWarnings(as.integer(args[5]))
quantile_basis_input <- args[6]

if (is.na(n_quantile_groups) || n_quantile_groups < 2) {
  stop("Argument <n_quantile_groups> must be an integer >= 2, e.g. 5 or 10.")
}

parse_quantile_basis <- function(x) {
  x_clean <- tolower(x)
  x_clean <- gsub("[^a-z0-9]", "", x_clean)

  if (x_clean %in% c("control", "controls", "ctrl")) {
    return("controls")
  }

  if (x_clean %in% c("all", "combined", "casescontrols", "casecontrol", "casesandcontrols", "both")) {
    return("all")
  }

  stop(
    paste0(
      "Invalid quantile_basis: ", x,
      "\nAllowed values: controls, control, all, combined, cases+controls, case_control"
    )
  )
}

quantile_basis <- parse_quantile_basis(quantile_basis_input)
quantile_basis_label <- ifelse(quantile_basis == "controls", "controls", "cases+controls")

quantile_prefix <- ifelse(
  n_quantile_groups == 10,
  "Decile",
  ifelse(n_quantile_groups == 5, "Quantile", "Quantile")
)

quantile_out_tag <- paste0(n_quantile_groups, "groups_", quantile_basis)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#############################################################################
# Required packages
#############################################################################

required_packages <- c("data.table", "tidyverse", "pROC", "broom")

install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Package '%s' not found. Attempting to install...", pkg))
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
        message(sprintf("Successfully installed '%s'.", pkg))
      }, error = function(e) {
        stop(sprintf("Failed to install package '%s': %s", pkg, e$message))
      })
    }

    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }
}

install_if_missing(required_packages)

#############################################################################
# Helper functions
#############################################################################

write_text_output <- function(x, file) {
  writeLines(capture.output(x), con = file)
}

stop_if_missing_cols <- function(data, required_cols, data_name) {
  missing_cols <- setdiff(required_cols, colnames(data))

  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required columns in ", data_name, ": ",
        paste(missing_cols, collapse = ", "),
        "\nAvailable columns: ",
        paste(colnames(data), collapse = ", ")
      )
    )
  }
}

make_prs_quantile <- function(score,
                              phenotype_binary,
                              n = 5,
                              basis = c("controls", "all"),
                              prefix = "Quantile") {
  basis <- match.arg(basis)

  if (basis == "controls") {
    reference_scores <- score[phenotype_binary == 0]
  } else {
    reference_scores <- score
  }

  reference_scores <- reference_scores[!is.na(reference_scores)]

  if (length(reference_scores) < 2) {
    stop("Not enough non-missing reference scores to create quantile groups.")
  }

  raw_breaks <- as.numeric(
    quantile(
      reference_scores,
      probs = seq(0, 1, length.out = n + 1),
      na.rm = TRUE,
      type = 7
    )
  )

  breaks <- unique(raw_breaks)

  if (length(breaks) < 2) {
    stop("Could not create quantile groups because all reference PRS values are identical.")
  }

  if ((length(breaks) - 1) < n) {
    warning(
      paste0(
        "Requested ", n, " groups, but only ", length(breaks) - 1,
        " groups could be created due to tied PRS values."
      )
    )
  }

  # Use -Inf/Inf so cases outside the control range are still assigned
  # when quantile_basis = controls.
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf

  labels <- paste0(prefix, " ", seq_len(length(breaks) - 1))

  cut(
    score,
    breaks = breaks,
    include.lowest = TRUE,
    labels = labels,
    ordered_result = TRUE
  )
}

plot_roc_comparison <- function(roc_prs, roc_base, roc_full, title, file) {
  pdf(file)
  plot(roc_prs, col = "blue", main = title, lwd = 2)
  plot(roc_base, col = "green", add = TRUE, lwd = 2)
  plot(roc_full, col = "red", add = TRUE, lwd = 2)
  legend(
    "bottomright",
    legend = c("PRS only", "Age + Sex", "Age + Sex + PRS"),
    col = c("blue", "green", "red"),
    lwd = 2
  )
  dev.off()
}

make_quantile_or_table <- function(df, quantile_col, model_formula) {
  df <- df %>%
    filter(!is.na(.data[[quantile_col]])) %>%
    mutate(
      !!quantile_col := droplevels(.data[[quantile_col]])
    )

  quantile_levels <- levels(df[[quantile_col]])

  if (length(quantile_levels) < 2) {
    stop(paste0("Only one level present in ", quantile_col, "; cannot fit quantile OR model."))
  }

  logit_quantile <- glm(
    model_formula,
    data = df,
    family = binomial()
  )

  or_table <- broom::tidy(
    logit_quantile,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(str_detect(term, quantile_col)) %>%
    mutate(
      Quantile = str_replace(term, paste0("^", quantile_col), ""),
      Quantile = str_replace_all(Quantile, "`", ""),
      Quantile = str_trim(Quantile)
    ) %>%
    rename(
      OR = estimate,
      OR_lower_CI = conf.low,
      OR_upper_CI = conf.high
    ) %>%
    select(Quantile, OR, OR_lower_CI, OR_upper_CI)

  count_table <- df %>%
    group_by(Quantile = .data[[quantile_col]]) %>%
    summarise(
      Cases = sum(phenotype_binary == 1),
      Controls = sum(phenotype_binary == 0),
      .groups = "drop"
    ) %>%
    mutate(Quantile = as.character(Quantile))

  final_table <- full_join(count_table, or_table, by = "Quantile") %>%
    mutate(
      Quantile = factor(Quantile, levels = quantile_levels)
    ) %>%
    arrange(Quantile) %>%
    mutate(
      Quantile = as.character(Quantile)
    )

  ref_quantile <- quantile_levels[1]

  final_table <- final_table %>%
    mutate(
      OR = ifelse(is.na(OR) & Quantile == ref_quantile, 1, OR),
      OR_lower_CI = ifelse(is.na(OR_lower_CI) & Quantile == ref_quantile, 1, OR_lower_CI),
      OR_upper_CI = ifelse(is.na(OR_upper_CI) & Quantile == ref_quantile, 1, OR_upper_CI)
    )

  list(
    model = logit_quantile,
    table = final_table
  )
}

#############################################################################
# Load input data
#############################################################################

message("Loading PRS file: ", prs_file)

prs <- fread(prs_file) %>%
  rename_with(~ gsub("^#", "", .x))

required_prs_cols <- c("IID", "WEIGHTED_SUM", "UNWEIGHTED_SUM")
stop_if_missing_cols(prs, required_prs_cols, "PRS .sscore file")

prs <- prs %>%
  mutate(IID = as.character(IID)) %>%
  select(IID, WEIGHTED_SUM, UNWEIGHTED_SUM) %>%
  rename(
    PRS_weighted = WEIGHTED_SUM,
    PRS_unweighted = UNWEIGHTED_SUM
  )

message("Loading phenotype file: ", pheno_file)

pheno <- fread(pheno_file) %>%
  rename_with(~ gsub("^#", "", .x))

if (!("IID" %in% colnames(pheno))) {
  stop(
    paste0(
      "Phenotype file must contain an IID column.\nAvailable columns: ",
      paste(colnames(pheno), collapse = ", ")
    )
  )
}

if (!(phenocol_name %in% colnames(pheno))) {
  stop(
    paste0(
      "Phenotype column '", phenocol_name, "' not found in phenotype file.\nAvailable columns: ",
      paste(colnames(pheno), collapse = ", ")
    )
  )
}

pheno <- pheno %>%
  mutate(IID = as.character(IID)) %>%
  rename(phenotype_col = all_of(phenocol_name))

#############################################################################
# Merge PRS and phenotype
#############################################################################

df <- inner_join(prs, pheno, by = "IID")

if (nrow(df) == 0) {
  stop("No overlapping IIDs between PRS file and phenotype file.")
}

message("Merged sample size: ", nrow(df))

df <- df %>%
  mutate(
    phenotype_col = as.numeric(phenotype_col)
  ) %>%
  filter(phenotype_col %in% c(0, 1)) %>%
  mutate(
    phenotype_binary = phenotype_col,
    phenotype = factor(
      phenotype_binary,
      levels = c(0, 1),
      labels = c("Control", "Case")
    )
  )

if (nrow(df) == 0) {
  stop("No samples with phenotype coded as 0/1 after filtering.")
}

message("Samples with valid 0/1 phenotype: ", nrow(df))
message("Cases: ", sum(df$phenotype_binary == 1))
message("Controls: ", sum(df$phenotype_binary == 0))

#############################################################################
# Check required covariates
#############################################################################

pc_terms <- paste0("PC", 1:10)

required_model_cols <- c(
  "IID",
  "PRS_weighted",
  "PRS_unweighted",
  "phenotype_binary",
  "phenotype",
  "Age",
  "Sex",
  pc_terms
)

stop_if_missing_cols(df, required_model_cols, "merged PRS/phenotype data")

df <- df %>%
  mutate(
    Age = as.numeric(Age),
    Sex = factor(Sex),
    across(all_of(pc_terms), as.numeric),
    PRS_weighted = as.numeric(PRS_weighted),
    PRS_unweighted = as.numeric(PRS_unweighted)
  )

#############################################################################
# Complete-case data for modelling
#############################################################################

model_cols <- c(
  "IID",
  "PRS_weighted",
  "PRS_unweighted",
  "phenotype_binary",
  "phenotype",
  "Age",
  "Sex",
  pc_terms
)

df_model <- df %>%
  select(all_of(model_cols)) %>%
  drop_na()

if (nrow(df_model) == 0) {
  stop("No complete cases available for modelling after removing missing PRS/covariate values.")
}

message("Complete-case sample size for modelling: ", nrow(df_model))
message("Complete-case cases: ", sum(df_model$phenotype_binary == 1))
message("Complete-case controls: ", sum(df_model$phenotype_binary == 0))
message("PRS grouping: ", n_quantile_groups, " groups")
message("PRS grouping basis: ", quantile_basis_label)

if (length(unique(df_model$phenotype_binary)) < 2) {
  stop("Complete-case model data does not contain both cases and controls.")
}

#############################################################################
# Basic summary output
#############################################################################

summary_file <- file.path(out_dir, paste0(phenocol_name, "_analysis_summary.txt"))

writeLines(
  c(
    paste0("PRS file: ", prs_file),
    paste0("Phenotype file: ", pheno_file),
    paste0("Phenotype column: ", phenocol_name),
    paste0("Merged sample size: ", nrow(df)),
    paste0("Valid phenotype sample size: ", nrow(df)),
    paste0("Complete-case model sample size: ", nrow(df_model)),
    paste0("Cases complete-case: ", sum(df_model$phenotype_binary == 1)),
    paste0("Controls complete-case: ", sum(df_model$phenotype_binary == 0)),
    paste0("PRS grouping n: ", n_quantile_groups),
    paste0("PRS grouping basis: ", quantile_basis_label),
    "",
    "PRS columns used:",
    "PRS_weighted   = WEIGHTED_SUM",
    "PRS_unweighted = UNWEIGHTED_SUM"
  ),
  con = summary_file
)

#############################################################################
# Plot PRS distributions
#############################################################################

pdf(file.path(out_dir, paste0(phenocol_name, "_weighted_prs_distribution.pdf")))
p <- ggplot(df, aes(x = PRS_weighted, fill = phenotype)) +
  geom_density(alpha = 0.6) +
  labs(
    title = "Weighted PRS Distribution: Cases vs Controls",
    x = "Weighted PRS",
    fill = "Phenotype"
  ) +
  theme_minimal()
print(p)
dev.off()

pdf(file.path(out_dir, paste0(phenocol_name, "_unweighted_prs_distribution.pdf")))
p <- ggplot(df, aes(x = PRS_unweighted, fill = phenotype)) +
  geom_density(alpha = 0.6, adjust = 4) +
  labs(
    title = "Unweighted PRS Distribution: Cases vs Controls",
    x = "Unweighted effect-allele count",
    fill = "Phenotype"
  ) +
  theme_minimal()

print(p)
dev.off()

#############################################################################
# Boxplots
#############################################################################

pdf(file.path(out_dir, paste0(phenocol_name, "_weighted_prs_boxplot.pdf")))
p <- ggplot(df, aes(x = phenotype, y = PRS_weighted, fill = phenotype)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Weighted PRS by Case-Control Status",
    x = "Phenotype",
    y = "Weighted PRS"
  ) +
  theme_minimal()
print(p)
dev.off()

pdf(file.path(out_dir, paste0(phenocol_name, "_unweighted_prs_boxplot.pdf")))
p <- ggplot(df, aes(x = phenotype, y = PRS_unweighted, fill = phenotype)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Unweighted PRS by Case-Control Status",
    x = "Phenotype",
    y = "Unweighted effect-allele count"
  ) +
  theme_minimal()
print(p)
dev.off()

#############################################################################
# Wilcoxon tests
#############################################################################

weighted_wilcox <- wilcox.test(PRS_weighted ~ phenotype, data = df)
unweighted_wilcox <- wilcox.test(PRS_unweighted ~ phenotype, data = df)

write_text_output(
  weighted_wilcox,
  file.path(out_dir, paste0(phenocol_name, "_weighted_prs_wilcoxon_test.txt"))
)

write_text_output(
  unweighted_wilcox,
  file.path(out_dir, paste0(phenocol_name, "_unweighted_prs_wilcoxon_test.txt"))
)

#############################################################################
# AUC comparison: PRS only, Age+Sex, Age+Sex+PRS
#############################################################################

model_prs_only_weighted <- glm(
  phenotype_binary ~ PRS_weighted,
  data = df_model,
  family = binomial()
)

model_base_weighted <- glm(
  phenotype_binary ~ Age + Sex,
  data = df_model,
  family = binomial()
)

model_full_weighted <- glm(
  phenotype_binary ~ PRS_weighted + Age + Sex,
  data = df_model,
  family = binomial()
)

roc_prs_only_weighted <- pROC::roc(
  response = df_model$phenotype_binary,
  predictor = fitted(model_prs_only_weighted),
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

roc_base_weighted <- pROC::roc(
  response = df_model$phenotype_binary,
  predictor = fitted(model_base_weighted),
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

roc_full_weighted <- pROC::roc(
  response = df_model$phenotype_binary,
  predictor = fitted(model_full_weighted),
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

auc_comparison_weighted <- data.frame(
  Model = c("PRS_only", "Age+Sex", "Age+Sex+PRS"),
  AUC = c(
    as.numeric(pROC::auc(roc_prs_only_weighted)),
    as.numeric(pROC::auc(roc_base_weighted)),
    as.numeric(pROC::auc(roc_full_weighted))
  )
)

fwrite(
  auc_comparison_weighted,
  file.path(out_dir, paste0(phenocol_name, "_weighted_auc_model_comparison.tsv")),
  sep = "\t"
)

plot_roc_comparison(
  roc_prs = roc_prs_only_weighted,
  roc_base = roc_base_weighted,
  roc_full = roc_full_weighted,
  title = "AUC Comparison: Weighted PRS",
  file = file.path(out_dir, paste0(phenocol_name, "_weighted_auc_model_comparison_plot.pdf"))
)

model_prs_only_unweighted <- glm(
  phenotype_binary ~ PRS_unweighted,
  data = df_model,
  family = binomial()
)

model_base_unweighted <- glm(
  phenotype_binary ~ Age + Sex,
  data = df_model,
  family = binomial()
)

model_full_unweighted <- glm(
  phenotype_binary ~ PRS_unweighted + Age + Sex,
  data = df_model,
  family = binomial()
)

roc_prs_only_unweighted <- pROC::roc(
  response = df_model$phenotype_binary,
  predictor = fitted(model_prs_only_unweighted),
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

roc_base_unweighted <- pROC::roc(
  response = df_model$phenotype_binary,
  predictor = fitted(model_base_unweighted),
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

roc_full_unweighted <- pROC::roc(
  response = df_model$phenotype_binary,
  predictor = fitted(model_full_unweighted),
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

auc_comparison_unweighted <- data.frame(
  Model = c("PRS_only", "Age+Sex", "Age+Sex+PRS"),
  AUC = c(
    as.numeric(pROC::auc(roc_prs_only_unweighted)),
    as.numeric(pROC::auc(roc_base_unweighted)),
    as.numeric(pROC::auc(roc_full_unweighted))
  )
)

fwrite(
  auc_comparison_unweighted,
  file.path(out_dir, paste0(phenocol_name, "_unweighted_auc_model_comparison.tsv")),
  sep = "\t"
)

plot_roc_comparison(
  roc_prs = roc_prs_only_unweighted,
  roc_base = roc_base_unweighted,
  roc_full = roc_full_unweighted,
  title = "AUC Comparison: Unweighted PRS",
  file = file.path(out_dir, paste0(phenocol_name, "_unweighted_auc_model_comparison_plot.pdf"))
)

#############################################################################
# AUC comparison: DeLong tests
#############################################################################

delong_weighted_prs_vs_full <- pROC::roc.test(
  roc_prs_only_weighted,
  roc_full_weighted,
  method = "delong"
)

delong_weighted_base_vs_full <- pROC::roc.test(
  roc_base_weighted,
  roc_full_weighted,
  method = "delong"
)

writeLines(
  c(
    "DeLong Test: PRS only vs Age+Sex+PRS",
    capture.output(delong_weighted_prs_vs_full),
    "",
    "DeLong Test: Age+Sex vs Age+Sex+PRS",
    capture.output(delong_weighted_base_vs_full)
  ),
  con = file.path(out_dir, paste0(phenocol_name, "_weighted_auc_delong_tests.txt"))
)

delong_unweighted_prs_vs_full <- pROC::roc.test(
  roc_prs_only_unweighted,
  roc_full_unweighted,
  method = "delong"
)

delong_unweighted_base_vs_full <- pROC::roc.test(
  roc_base_unweighted,
  roc_full_unweighted,
  method = "delong"
)

writeLines(
  c(
    "DeLong Test: PRS only vs Age+Sex+PRS",
    capture.output(delong_unweighted_prs_vs_full),
    "",
    "DeLong Test: Age+Sex vs Age+Sex+PRS",
    capture.output(delong_unweighted_base_vs_full)
  ),
  con = file.path(out_dir, paste0(phenocol_name, "_unweighted_auc_delong_tests.txt"))
)

#############################################################################
# Logistic regression with Age, Sex, and PCs
#############################################################################

weighted_formula <- as.formula(
  paste(
    "phenotype_binary ~ PRS_weighted + Age + Sex +",
    paste(pc_terms, collapse = " + ")
  )
)

unweighted_formula <- as.formula(
  paste(
    "phenotype_binary ~ PRS_unweighted + Age + Sex +",
    paste(pc_terms, collapse = " + ")
  )
)

glm_weighted <- glm(
  weighted_formula,
  data = df_model,
  family = binomial()
)

glm_unweighted <- glm(
  unweighted_formula,
  data = df_model,
  family = binomial()
)

write_text_output(
  summary(glm_weighted),
  file.path(out_dir, paste0(phenocol_name, "_weighted_logistic_regression.txt"))
)

write_text_output(
  summary(glm_unweighted),
  file.path(out_dir, paste0(phenocol_name, "_unweighted_logistic_regression.txt"))
)

#############################################################################
# Quantile-based ORs
#############################################################################

df_model <- df_model %>%
  mutate(
    quantile_weighted = make_prs_quantile(
      score = PRS_weighted,
      phenotype_binary = phenotype_binary,
      n = n_quantile_groups,
      basis = quantile_basis,
      prefix = quantile_prefix
    )
  )

df_model$quantile_weighted <- droplevels(df_model$quantile_weighted)

df_model$quantile_weighted <- relevel(
  df_model$quantile_weighted,
  ref = levels(df_model$quantile_weighted)[1]
)

weighted_quantile_formula <- as.formula(
  paste(
    "phenotype_binary ~ quantile_weighted + Age + Sex +",
    paste(pc_terms, collapse = " + ")
  )
)

weighted_quantile_result <- make_quantile_or_table(
  df = df_model,
  quantile_col = "quantile_weighted",
  model_formula = weighted_quantile_formula
)

or_table_weighted_final <- weighted_quantile_result$table

fwrite(
  or_table_weighted_final,
  file.path(
    out_dir,
    paste0(
      phenocol_name,
      "_weighted_odds_ratios_by_",
      quantile_out_tag,
      ".tsv"
    )
  ),
  sep = "\t"
)

pdf(
  file.path(
    out_dir,
    paste0(
      phenocol_name,
      "_weighted_",
      quantile_out_tag,
      "_or_forestplot.pdf"
    )
  ),
  width = 7,
  height = 5
)

p <- ggplot(
  or_table_weighted_final,
  aes(x = Quantile, y = OR, ymin = OR_lower_CI, ymax = OR_upper_CI)
) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = paste0(
      "Odds Ratios by Weighted PRS Groups\n",
      n_quantile_groups, " groups based on ", quantile_basis_label
    ),
    x = paste0("Weighted PRS group"),
    y = "OR (95% CI)"
  ) +
  theme_minimal()

print(p)
dev.off()

n_unique_unweighted <- n_distinct(df_model$PRS_unweighted)

if (n_unique_unweighted < n_quantile_groups) {
  warning(
    paste0(
      "Unweighted PRS has only ", n_unique_unweighted,
      " unique values, fewer than requested groups n = ", n_quantile_groups,
      ". Quantile groups may be coarse because tied values are not split."
    )
  )
}

df_model <- df_model %>%
  mutate(
    quantile_unweighted = make_prs_quantile(
      score = PRS_unweighted,
      phenotype_binary = phenotype_binary,
      n = n_quantile_groups,
      basis = quantile_basis,
      prefix = quantile_prefix
    )
  )

df_model$quantile_unweighted <- droplevels(df_model$quantile_unweighted)

df_model$quantile_unweighted <- relevel(
  df_model$quantile_unweighted,
  ref = levels(df_model$quantile_unweighted)[1]
)

unweighted_quantile_formula <- as.formula(
  paste(
    "phenotype_binary ~ quantile_unweighted + Age + Sex +",
    paste(pc_terms, collapse = " + ")
  )
)

unweighted_quantile_result <- make_quantile_or_table(
  df = df_model,
  quantile_col = "quantile_unweighted",
  model_formula = unweighted_quantile_formula
)

or_table_unweighted_final <- unweighted_quantile_result$table

fwrite(
  or_table_unweighted_final,
  file.path(
    out_dir,
    paste0(
      phenocol_name,
      "_unweighted_odds_ratios_by_",
      quantile_out_tag,
      ".tsv"
    )
  ),
  sep = "\t"
)

pdf(
  file.path(
    out_dir,
    paste0(
      phenocol_name,
      "_unweighted_",
      quantile_out_tag,
      "_or_forestplot.pdf"
    )
  ),
  width = 7,
  height = 5
)

p <- ggplot(
  or_table_unweighted_final,
  aes(x = Quantile, y = OR, ymin = OR_lower_CI, ymax = OR_upper_CI)
) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = paste0(
      "Odds Ratios by Unweighted PRS Groups\n",
      n_quantile_groups, " groups based on ", quantile_basis_label
    ),
    x = paste0("Unweighted PRS group"),
    y = "OR (95% CI)"
  ) +
  theme_minimal()

print(p)
dev.off()

#############################################################################
# End
#############################################################################

message("Downstream PRS analysis complete.")
message("Output directory: ", out_dir)
message("PRS grouping n: ", n_quantile_groups)
message("PRS grouping basis: ", quantile_basis_label)