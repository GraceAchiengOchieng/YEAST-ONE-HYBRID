# ============================================================
# STEP 2: Filter the 100 TaDof genes for detectable
# root expression using RNA-seq counts
#
# Retention rule:
#   Overall mean CPM >= 1
#          OR
#   Mean CPM >= 2 in at least one biological group
#
# Input files:
#   1. wheat_counts_for_analysis.csv
#   2. TaDof_100_genes_RefSeqv2.1.csv
#
# Output:
#   Root-expressed TaDof genes and expression statistics
# ============================================================


# ------------------------------------------------------------
# 0. Load packages
# ------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(stringr)


# ------------------------------------------------------------
# 1. Ask user to select the RNA-seq count file
# ------------------------------------------------------------

cat("\nPlease select the wheat RNA-seq COUNT file...\n")

count_file <- file.choose()

cat("\nSelected count file:\n")
print(count_file)


# ------------------------------------------------------------
# 2. Ask user to select the 100-TaDof gene list from Step 1
# ------------------------------------------------------------

cat("\nPlease select the TaDof 100-gene list from STEP 1...\n")

dof_file <- file.choose()

cat("\nSelected TaDof file:\n")
print(dof_file)


# ------------------------------------------------------------
# 3. Read both files
# ------------------------------------------------------------

counts <- read_csv(count_file, show_col_types = FALSE)

dof100 <- read_csv(dof_file, show_col_types = FALSE)


# ------------------------------------------------------------
# 4. Inspect the files
# ------------------------------------------------------------

cat("\n============================================\n")
cat("INPUT FILE INFORMATION\n")
cat("============================================\n")

cat("Genes in count matrix:", nrow(counts), "\n")
cat("Columns in count matrix:", ncol(counts), "\n")

cat("\nFirst columns in count matrix:\n")
print(head(colnames(counts), 10))

cat("\nNumber of genes in TaDof list:", nrow(dof100), "\n")


# ------------------------------------------------------------
# 5. Rename the first count-matrix column to GeneID
#
# In your uploaded file this column is currently called:
# "Unnamed: 0"
# ------------------------------------------------------------

colnames(counts)[1] <- "GeneID"


# Rename first column of TaDof file as GeneID too
colnames(dof100)[1] <- "GeneID"


# ------------------------------------------------------------
# 6. Check for duplicated gene IDs
# ------------------------------------------------------------

if (anyDuplicated(counts$GeneID) > 0) {
  
  warning(
    "Duplicated gene IDs were found in the count matrix."
  )
  
}

if (anyDuplicated(dof100$GeneID) > 0) {
  
  warning(
    "Duplicated gene IDs were found in the TaDof list."
  )
  
}


# ------------------------------------------------------------
# 7. Keep only the 100 TaDof genes
# ------------------------------------------------------------

dof_counts <- counts %>%
  filter(GeneID %in% dof100$GeneID)


cat("\n============================================\n")
cat("MATCHING TaDof GENES\n")
cat("============================================\n")

cat(
  "TaDof genes found in count matrix:",
  nrow(dof_counts),
  "\n"
)


# Check whether any of the 100 genes are missing

missing_dofs <- setdiff(
  dof100$GeneID,
  dof_counts$GeneID
)

cat(
  "TaDof genes missing from count matrix:",
  length(missing_dofs),
  "\n"
)

if (length(missing_dofs) > 0) {
  print(missing_dofs)
}


# ------------------------------------------------------------
# 8. Identify sample columns
# ------------------------------------------------------------

sample_columns <- setdiff(
  colnames(dof_counts),
  "GeneID"
)


# ------------------------------------------------------------
# 9. Calculate library size for every RNA-seq sample
#
# IMPORTANT:
# CPM must be calculated using ALL genes in the library,
# not only the 100 Dof genes.
# ------------------------------------------------------------

library_sizes <- colSums(
  counts[, sample_columns],
  na.rm = TRUE
)


cat("\n============================================\n")
cat("LIBRARY SIZES\n")
cat("============================================\n")

print(summary(library_sizes))


# ------------------------------------------------------------
# 10. Calculate CPM for the 100 TaDof genes
#
# CPM = raw count / total sample reads * 1,000,000
# ------------------------------------------------------------

dof_cpm <- dof_counts

for (sample in sample_columns) {
  
  dof_cpm[[sample]] <-
    dof_counts[[sample]] /
    library_sizes[[sample]] *
    1000000
}


# ------------------------------------------------------------
# 11. Calculate overall mean CPM for each Dof
#
# This asks:
# Is the gene reasonably expressed across the complete
# root experiment?
# ------------------------------------------------------------

dof_cpm$Overall_Mean_CPM <- rowMeans(
  dof_cpm[, sample_columns],
  na.rm = TRUE
)


# ------------------------------------------------------------
# 12. Convert expression table to long format
# ------------------------------------------------------------

cpm_long <- dof_cpm %>%
  
  select(
    GeneID,
    all_of(sample_columns)
  ) %>%
  
  pivot_longer(
    cols = -GeneID,
    names_to = "Sample",
    values_to = "CPM"
  )


# ------------------------------------------------------------
# 13. Derive biological group from the sample names
#
# Examples:
#
# B1_1_filtered    -> B1
# B1_2_filtered    -> B1
#
# O3_1_filtered    -> O3
#
# N4_3_filtered    -> N4
#
# BO4_1_filtered   -> BO4
# BO4_1a_filtered  -> BO4
#
# This allows replicate CPM values to be averaged.
# ------------------------------------------------------------

cpm_long <- cpm_long %>%
  
  mutate(
    
    Sample_clean =
      str_remove(
        Sample,
        "_filtered$"
      ),
    
    Biological_Group =
      case_when(
        
        # Hybrid / BO samples
        str_detect(
          Sample_clean,
          "^BO[1-6]"
        ) ~
          str_extract(
            Sample_clean,
            "^BO[1-6]"
          ),
        
        # Basalt
        str_detect(
          Sample_clean,
          "^B[1-6]"
        ) ~
          str_extract(
            Sample_clean,
            "^B[1-6]"
          ),
        
        # Oakley
        str_detect(
          Sample_clean,
          "^O[1-6]"
        ) ~
          str_extract(
            Sample_clean,
            "^O[1-6]"
          ),
        
        # npf212 mutant
        str_detect(
          Sample_clean,
          "^N[1-6]"
        ) ~
          str_extract(
            Sample_clean,
            "^N[1-6]"
          ),
        
        # Anything that does not fit the naming system
        TRUE ~ NA_character_
      )
  )


# ------------------------------------------------------------
# 14. Show samples for which the group could not be determined
# ------------------------------------------------------------

unassigned_samples <- cpm_long %>%
  
  filter(
    is.na(Biological_Group)
  ) %>%
  
  distinct(Sample)


cat("\n============================================\n")
cat("UNASSIGNED SAMPLE NAMES\n")
cat("============================================\n")

print(unassigned_samples)


# ------------------------------------------------------------
# 15. Calculate mean CPM for every gene in every
# biological group
#
# Example:
#
# TraesCS...  B1   3.4 CPM
# TraesCS...  B2   0.8 CPM
# TraesCS...  O1   4.2 CPM
# etc.
# ------------------------------------------------------------

group_mean_cpm <- cpm_long %>%
  
  filter(
    !is.na(Biological_Group)
  ) %>%
  
  group_by(
    GeneID,
    Biological_Group
  ) %>%
  
  summarise(
    Mean_CPM = mean(
      CPM,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 16. Find the highest biological-group mean CPM
# for every Dof gene
# ------------------------------------------------------------

max_group_cpm <- group_mean_cpm %>%
  
  group_by(GeneID) %>%
  
  summarise(
    
    Max_Group_Mean_CPM =
      max(
        Mean_CPM,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 17. Identify WHICH group produced the maximum CPM
# ------------------------------------------------------------

max_group_identity <- group_mean_cpm %>%
  
  group_by(GeneID) %>%
  
  slice_max(
    order_by = Mean_CPM,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup() %>%
  
  rename(
    Highest_Expression_Group =
      Biological_Group,
    
    Highest_Group_Mean_CPM =
      Mean_CPM
  )


# ------------------------------------------------------------
# 18. Combine expression statistics
# ------------------------------------------------------------

expression_summary <- dof_cpm %>%
  
  select(
    GeneID,
    Overall_Mean_CPM
  ) %>%
  
  left_join(
    max_group_cpm,
    by = "GeneID"
  ) %>%
  
  left_join(
    max_group_identity,
    by = "GeneID"
  )


# ------------------------------------------------------------
# 19. Apply the root-expression filter
#
# PASS if:
#
# Overall mean CPM >= 1
#
# OR
#
# At least one biological group mean CPM >= 2
# ------------------------------------------------------------

expression_summary <- expression_summary %>%
  
  mutate(
    
    Pass_Overall_CPM =
      Overall_Mean_CPM >= 1,
    
    Pass_Group_CPM =
      Max_Group_Mean_CPM >= 2,
    
    Root_Expressed =
      Pass_Overall_CPM |
      Pass_Group_CPM
  )


# ------------------------------------------------------------
# 20. Extract the root-expressed Dofs
# ------------------------------------------------------------

root_expressed_dofs <- expression_summary %>%
  
  filter(
    Root_Expressed == TRUE
  ) %>%
  
  arrange(
    desc(Overall_Mean_CPM)
  )


# ------------------------------------------------------------
# 21. Extract genes failing the filter
# ------------------------------------------------------------

not_root_expressed <- expression_summary %>%
  
  filter(
    Root_Expressed == FALSE
  ) %>%
  
  arrange(
    desc(Overall_Mean_CPM)
  )


# ------------------------------------------------------------
# 22. Print results
# ------------------------------------------------------------

cat("\n\n============================================\n")
cat("STEP 2 RESULTS\n")
cat("============================================\n")

cat(
  "Starting TaDof genes:",
  nrow(dof100),
  "\n"
)

cat(
  "TaDof genes found in count matrix:",
  nrow(dof_counts),
  "\n"
)

cat(
  "Root-expressed TaDofs:",
  nrow(root_expressed_dofs),
  "\n"
)

cat(
  "TaDofs failing expression filter:",
  nrow(not_root_expressed),
  "\n"
)

cat("============================================\n")


# ------------------------------------------------------------
# 23. Display root-expressed genes in RStudio
# ------------------------------------------------------------

View(root_expressed_dofs)


# ------------------------------------------------------------
# 24. Create wide biological-group CPM table
#
# Useful for later analysis
# ------------------------------------------------------------

group_cpm_wide <- group_mean_cpm %>%
  
  pivot_wider(
    names_from = Biological_Group,
    values_from = Mean_CPM
  )


# ------------------------------------------------------------
# 25. Combine group CPMs with the filter statistics
# ------------------------------------------------------------

root_expression_full <- expression_summary %>%
  
  left_join(
    group_cpm_wide,
    by = "GeneID"
  ) %>%
  
  arrange(
    desc(Overall_Mean_CPM)
  )


# ------------------------------------------------------------
# 26. Ask user where results should be saved
# ------------------------------------------------------------

cat(
  "\nPlease choose a folder where STEP 2 results should be saved...\n"
)

output_folder <- choose.dir(
  caption =
    "Choose folder for Step 2 TaDof expression results"
)


# ------------------------------------------------------------
# 27. Save results
# ------------------------------------------------------------

write_csv(
  root_expressed_dofs,
  file.path(
    output_folder,
    "TaDof_root_expressed_genes.csv"
  )
)


write_csv(
  not_root_expressed,
  file.path(
    output_folder,
    "TaDof_not_root_expressed_genes.csv"
  )
)


write_csv(
  root_expression_full,
  file.path(
    output_folder,
    "TaDof_root_expression_full_summary.csv"
  )
)


write_csv(
  group_mean_cpm,
  file.path(
    output_folder,
    "TaDof_group_mean_CPM_long.csv"
  )
)


write_csv(
  group_cpm_wide,
  file.path(
    output_folder,
    "TaDof_group_mean_CPM_wide.csv"
  )
)


write_csv(
  unassigned_samples,
  file.path(
    output_folder,
    "unassigned_RNAseq_samples.csv"
  )
)


# ------------------------------------------------------------
# 28. Final message
# ------------------------------------------------------------

cat("\n\n============================================\n")
cat("STEP 2 COMPLETE\n")
cat("============================================\n")

cat(
  "Starting Dof genes:",
  nrow(dof100),
  "\n"
)

cat(
  "Root-expressed Dof genes:",
  nrow(root_expressed_dofs),
  "\n"
)

cat(
  "Not root-expressed:",
  nrow(not_root_expressed),
  "\n"
)

cat(
  "\nResults saved in:\n",
  output_folder,
  "\n"
)

cat("============================================\n")
