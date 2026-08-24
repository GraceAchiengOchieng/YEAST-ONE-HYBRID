# ============================================================
# STEP 3C: Identify highly expressed constitutive TaDof genes
#
# Goal:
#   Find root-expressed Dofs that are highly abundant across
#   the experiment, even if they are not strongly differential.
#
# Starting population:
#   Root-expressed TaDof genes from STEP 2
#
# Main selection:
#   Top 10% by Overall_Mean_CPM
#
# Additional characterization:
#   Expression stability across biological groups
# ============================================================


# ------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------

library(readr)
library(dplyr)


# ------------------------------------------------------------
# 1. Ask user to select STEP 2 root-expression summary file
# ------------------------------------------------------------

cat("\nPlease select the STEP 2 root-expression summary file...\n")
cat("Example: TaDof_root_expression_full_summary.csv\n\n")

root_file <- file.choose()

cat("\nSelected file:\n")
print(root_file)


# ------------------------------------------------------------
# 2. Read the file
# ------------------------------------------------------------

root_summary <- read_csv(
  root_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 3. Check required columns
# ------------------------------------------------------------

required_columns <- c(
  "GeneID",
  "Overall_Mean_CPM",
  "Root_Expressed"
)

missing_columns <- setdiff(
  required_columns,
  colnames(root_summary)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "Missing required columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 4. Keep only root-expressed Dofs
# ------------------------------------------------------------

root_expressed <- root_summary %>%
  
  filter(
    Root_Expressed == TRUE
  )


cat("\n============================================\n")
cat("STEP 3C INPUT\n")
cat("============================================\n")

cat(
  "Root-expressed TaDofs entering analysis:",
  nrow(root_expressed),
  "\n"
)


# ------------------------------------------------------------
# 5. Identify biological-group CPM columns
#
# We exclude summary/statistical columns and keep only
# columns containing group-level CPM values.
# ------------------------------------------------------------

non_group_columns <- c(
  "GeneID",
  "Overall_Mean_CPM",
  "Max_Group_Mean_CPM",
  "Highest_Expression_Group",
  "Highest_Group_Mean_CPM",
  "Pass_Overall_CPM",
  "Pass_Group_CPM",
  "Root_Expressed"
)

group_columns <- setdiff(
  colnames(root_expressed),
  non_group_columns
)


cat(
  "Biological-group CPM columns detected:",
  length(group_columns),
  "\n"
)


# ------------------------------------------------------------
# 6. Make sure detected group columns are numeric
# ------------------------------------------------------------

numeric_group_columns <- group_columns[
  sapply(
    root_expressed[group_columns],
    is.numeric
  )
]

cat(
  "Numeric biological-group columns used:",
  length(numeric_group_columns),
  "\n"
)


# ------------------------------------------------------------
# 7. Calculate expression stability metrics
#
# Mean_Group_CPM:
#   average expression across biological groups
#
# SD_Group_CPM:
#   standard deviation across groups
#
# CV:
#   coefficient of variation = SD / mean
#
# Lower CV = more stable / constitutive expression
# ------------------------------------------------------------

root_expression_stats <- root_expressed %>%
  
  rowwise() %>%
  
  mutate(
    
    Mean_Group_CPM =
      mean(
        c_across(
          all_of(numeric_group_columns)
        ),
        na.rm = TRUE
      ),
    
    SD_Group_CPM =
      sd(
        c_across(
          all_of(numeric_group_columns)
        ),
        na.rm = TRUE
      ),
    
    CV_Group_CPM =
      ifelse(
        Mean_Group_CPM > 0,
        SD_Group_CPM / Mean_Group_CPM,
        NA_real_
      ),
    
    Min_Group_CPM =
      min(
        c_across(
          all_of(numeric_group_columns)
        ),
        na.rm = TRUE
      ),
    
    Max_Group_CPM =
      max(
        c_across(
          all_of(numeric_group_columns)
        ),
        na.rm = TRUE
      )
    
  ) %>%
  
  ungroup()


# ------------------------------------------------------------
# 8. Rank genes by overall mean CPM
# ------------------------------------------------------------

root_expression_stats <- root_expression_stats %>%
  
  arrange(
    desc(Overall_Mean_CPM)
  ) %>%
  
  mutate(
    
    Expression_Rank =
      row_number()
    
  )


# ------------------------------------------------------------
# 9. Define top 10%
#
# For 64 genes this should retain about 6-7 genes.
# ceiling() ensures at least the upper 10% are retained.
# ------------------------------------------------------------

top_n <- ceiling(
  nrow(root_expression_stats) * 0.10
)

cat(
  "\nTop 10% corresponds to:",
  top_n,
  "genes\n"
)


# ------------------------------------------------------------
# 10. Select the high-expression safety-net candidates
# ------------------------------------------------------------

high_expression_candidates <- root_expression_stats %>%
  
  slice_head(
    n = top_n
  ) %>%
  
  mutate(
    High_Expression_Safety_Net = TRUE
  )


# ------------------------------------------------------------
# 11. Create a clean results table
# ------------------------------------------------------------

high_expression_results <- high_expression_candidates %>%
  
  select(
    GeneID,
    Overall_Mean_CPM,
    Mean_Group_CPM,
    SD_Group_CPM,
    CV_Group_CPM,
    Min_Group_CPM,
    Max_Group_CPM,
    Expression_Rank
  )


# ------------------------------------------------------------
# 12. Rank the top-expression genes by stability as well
#
# Lower CV = more stable expression across biological groups
# ------------------------------------------------------------

high_expression_by_stability <- high_expression_results %>%
  
  arrange(
    CV_Group_CPM
  )


# ------------------------------------------------------------
# 13. Print results
# ------------------------------------------------------------

cat("\n\n============================================\n")
cat("STEP 3C RESULTS\n")
cat("============================================\n")

cat(
  "Root-expressed Dofs tested:",
  nrow(root_expression_stats),
  "\n"
)

cat(
  "Top 10% high-expression Dofs:",
  nrow(high_expression_results),
  "\n"
)

cat("============================================\n\n")


print(
  high_expression_results %>%
    select(
      GeneID,
      Overall_Mean_CPM,
      CV_Group_CPM,
      Expression_Rank
    )
)


# ------------------------------------------------------------
# 14. View results in RStudio
# ------------------------------------------------------------

View(high_expression_results)

View(high_expression_by_stability)


# ------------------------------------------------------------
# 15. Save full ranked root-expressed list
# ------------------------------------------------------------

full_ranked_table <- root_expression_stats %>%
  
  select(
    GeneID,
    Overall_Mean_CPM,
    Mean_Group_CPM,
    SD_Group_CPM,
    CV_Group_CPM,
    Min_Group_CPM,
    Max_Group_CPM,
    Expression_Rank,
    everything()
  )


# ------------------------------------------------------------
# 16. Ask where results should be saved
# ------------------------------------------------------------

cat(
  "\nPlease choose a folder where STEP 3C results should be saved...\n"
)

output_folder <- choose.dir(
  caption =
    "Choose folder for Step 3C high-expression TaDof results"
)


# ------------------------------------------------------------
# 17. Save results
# ------------------------------------------------------------

write_csv(
  high_expression_results,
  file.path(
    output_folder,
    "TaDof_Step3C_high_expression_candidates.csv"
  )
)


write_csv(
  high_expression_by_stability,
  file.path(
    output_folder,
    "TaDof_Step3C_high_expression_by_stability.csv"
  )
)


write_csv(
  full_ranked_table,
  file.path(
    output_folder,
    "TaDof_Step3C_all_root_expressed_ranked.csv"
  )
)


# ------------------------------------------------------------
# 18. Final confirmation
# ------------------------------------------------------------

cat("\n\n============================================\n")
cat("STEP 3C COMPLETE\n")
cat("============================================\n")

cat(
  "Starting root-expressed Dofs:",
  nrow(root_expression_stats),
  "\n"
)

cat(
  "High-expression safety-net candidates:",
  nrow(high_expression_results),
  "\n"
)

cat(
  "\nResults saved in:\n",
  output_folder,
  "\n"
)

cat("============================================\n")
