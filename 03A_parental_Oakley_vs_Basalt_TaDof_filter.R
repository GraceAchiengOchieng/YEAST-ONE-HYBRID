# ============================================================
# STEP 3A: Identify parental/background TaDof candidates
# Oakley vs Basalt under Low N
#
# Starting population:
#   64 root-expressed TaDof genes
#
# Selection criteria:
#   1. Mean absolute log2(Oakley/Basalt) >= 1
#      approximately a >=2-fold difference
#
#   2. Direction is consistent across T1-T4:
#      Oakley > Basalt
#             OR
#      Basalt > Oakley
#
#      "Mixed" genes are excluded
#
# Inputs:
#   1. TaDof_root_expression_full_summary.csv
#   2. TaDof_exhaustive_candidate_analysis.csv
# ============================================================


# ------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------

library(readr)
library(dplyr)


# ------------------------------------------------------------
# 1. Ask user to select the root-expression summary
# ------------------------------------------------------------

cat("\nPlease select the STEP 2 root-expression summary file...\n")
cat("Example: TaDof_root_expression_full_summary.csv\n\n")

root_file <- file.choose()

cat("\nSelected root-expression file:\n")
print(root_file)


# ------------------------------------------------------------
# 2. Ask user to select the exhaustive TaDof analysis file
# ------------------------------------------------------------

cat("\nPlease select the exhaustive TaDof candidate analysis file...\n")
cat("Example: TaDof_exhaustive_candidate_analysis.csv\n\n")

candidate_file <- file.choose()

cat("\nSelected candidate-analysis file:\n")
print(candidate_file)


# ------------------------------------------------------------
# 3. Read the files
# ------------------------------------------------------------

root_summary <- read_csv(
  root_file,
  show_col_types = FALSE
)

candidate_table <- read_csv(
  candidate_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 4. Basic checks
# ------------------------------------------------------------

cat("\n============================================\n")
cat("INPUT CHECK\n")
cat("============================================\n")

cat(
  "Genes in root-expression summary:",
  nrow(root_summary),
  "\n"
)

cat(
  "Genes in exhaustive candidate table:",
  nrow(candidate_table),
  "\n"
)


# ------------------------------------------------------------
# 5. Confirm required columns exist
# ------------------------------------------------------------

required_root_columns <- c(
  "GeneID",
  "Root_Expressed"
)

required_candidate_columns <- c(
  "GeneID",
  "Oakley_vs_Basalt_mean_abs_log2ratio",
  "Oakley_vs_Basalt_consistency",
  "Oakley_vs_Basalt_LowN_T1_log2ratio",
  "Oakley_vs_Basalt_LowN_T2_log2ratio",
  "Oakley_vs_Basalt_LowN_T3_log2ratio",
  "Oakley_vs_Basalt_LowN_T4_log2ratio"
)


missing_root_columns <- setdiff(
  required_root_columns,
  colnames(root_summary)
)

missing_candidate_columns <- setdiff(
  required_candidate_columns,
  colnames(candidate_table)
)


if (length(missing_root_columns) > 0) {
  stop(
    paste(
      "Missing columns in root-expression file:",
      paste(missing_root_columns, collapse = ", ")
    )
  )
}


if (length(missing_candidate_columns) > 0) {
  stop(
    paste(
      "Missing columns in exhaustive candidate file:",
      paste(missing_candidate_columns, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 6. Keep ONLY the root-expressed Dofs
# ------------------------------------------------------------

root_expressed <- root_summary %>%
  filter(Root_Expressed == TRUE) %>%
  select(GeneID)


cat("\nRoot-expressed TaDofs entering Step 3A:",
    nrow(root_expressed), "\n")


# ------------------------------------------------------------
# 7. Join the 64 root-expressed genes to the
# exhaustive candidate table
# ------------------------------------------------------------

parental_input <- candidate_table %>%
  inner_join(
    root_expressed,
    by = "GeneID"
  )


cat(
  "Root-expressed Dofs successfully matched:",
  nrow(parental_input),
  "\n"
)


# ------------------------------------------------------------
# 8. Apply parental/background criteria
#
# Criterion 1:
# Mean absolute Oakley/Basalt log2 ratio >= 1
#
# Criterion 2:
# Consistent direction across time
# Oakley>Basalt OR Basalt>Oakley
#
# Mixed genes are excluded
# ------------------------------------------------------------

parental_candidates <- parental_input %>%
  
  filter(
    Oakley_vs_Basalt_mean_abs_log2ratio >= 1,
    Oakley_vs_Basalt_consistency != "Mixed"
  ) %>%
  
  arrange(
    desc(Oakley_vs_Basalt_mean_abs_log2ratio)
  )


# ------------------------------------------------------------
# 9. Create a clean results table
# ------------------------------------------------------------

parental_results <- parental_candidates %>%
  
  select(
    GeneID,
    Description,
    Mean_CPM_all,
    Basalt_LowN_mean_CPM,
    Oakley_LowN_mean_CPM,
    Oakley_vs_Basalt_mean_abs_log2ratio,
    Oakley_vs_Basalt_consistency,
    Oakley_vs_Basalt_LowN_T1_log2ratio,
    Oakley_vs_Basalt_LowN_T2_log2ratio,
    Oakley_vs_Basalt_LowN_T3_log2ratio,
    Oakley_vs_Basalt_LowN_T4_log2ratio
  )


# ------------------------------------------------------------
# 10. Print the parental candidates
# ------------------------------------------------------------

cat("\n\n============================================\n")
cat("STEP 3A RESULTS\n")
cat("============================================\n")

cat(
  "Root-expressed Dofs tested:",
  nrow(parental_input),
  "\n"
)

cat(
  "Parental/background candidates:",
  nrow(parental_results),
  "\n"
)

cat("============================================\n\n")


print(
  parental_results %>%
    select(
      GeneID,
      Oakley_vs_Basalt_mean_abs_log2ratio,
      Oakley_vs_Basalt_consistency
    )
)


# ------------------------------------------------------------
# 11. View full table in RStudio
# ------------------------------------------------------------

View(parental_results)


# ------------------------------------------------------------
# 12. Also identify genes that failed this route
#
# They are NOT discarded from the overall analysis.
# They may still qualify through hybrid, nitrogen-response,
# or high-expression routes.
# ------------------------------------------------------------

non_parental_candidates <- parental_input %>%
  
  filter(
    !(
      Oakley_vs_Basalt_mean_abs_log2ratio >= 1 &
      Oakley_vs_Basalt_consistency != "Mixed"
    )
  )


# ------------------------------------------------------------
# 13. Ask where to save results
# ------------------------------------------------------------

cat(
  "\nPlease choose a folder where STEP 3A results should be saved...\n"
)

output_folder <- choose.dir(
  caption =
    "Choose folder for Step 3A parental TaDof results"
)


# ------------------------------------------------------------
# 14. Save results
# ------------------------------------------------------------

write_csv(
  parental_results,
  file.path(
    output_folder,
    "TaDof_Step3A_parental_candidates.csv"
  )
)


write_csv(
  non_parental_candidates,
  file.path(
    output_folder,
    "TaDof_Step3A_not_parental_candidates.csv"
  )
)


# ------------------------------------------------------------
# 15. Final confirmation
# ------------------------------------------------------------

cat("\n\n============================================\n")
cat("STEP 3A COMPLETE\n")
cat("============================================\n")

cat(
  "Starting root-expressed Dofs:",
  nrow(parental_input),
  "\n"
)

cat(
  "Parental/background candidates:",
  nrow(parental_results),
  "\n"
)

cat(
  "Other root-expressed Dofs:",
  nrow(non_parental_candidates),
  "\n"
)

cat(
  "\nResults saved in:\n",
  output_folder,
  "\n"
)

cat("============================================\n")
