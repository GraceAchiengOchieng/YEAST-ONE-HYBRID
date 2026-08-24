# ============================================================
# STEP 5A: Map IWGSC RefSeq v2.1 TaDof IDs to older RefSeq IDs
#
# Goal:
#   Add older RefSeq / Ensembl IDs to the 17 candidate genes
#   using the official IWGSC correspondence table.
#
# Inputs:
#   1. TaDof_Step4_master_17_candidates.csv
#   2. Official IWGSC RefSeq correspondence file
# ============================================================

library(readr)
library(dplyr)

# ------------------------------------------------------------
# 1. Choose the 17-gene master candidate table
# ------------------------------------------------------------

cat("\nPlease select the STEP 4 master candidate table...\n")

master_file <- file.choose()

cat("\nSelected master file:\n")
print(master_file)


# ------------------------------------------------------------
# 2. Choose the official IWGSC correspondence file
# ------------------------------------------------------------

cat("\nPlease select the official IWGSC RefSeq correspondence file...\n")

mapping_file <- file.choose()

cat("\nSelected correspondence file:\n")
print(mapping_file)


# ------------------------------------------------------------
# 3. Read both files
# ------------------------------------------------------------

master <- read_csv(
  master_file,
  show_col_types = FALSE
)

mapping <- read_csv(
  mapping_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 4. Inspect correspondence-file columns
# ------------------------------------------------------------

cat("\n============================================\n")
cat("CORRESPONDENCE FILE COLUMNS\n")
cat("============================================\n\n")

print(colnames(mapping))


# ------------------------------------------------------------
# 5. Interactive selection of the v2.1 column
# ------------------------------------------------------------

cat("\nSelect the column containing the CURRENT RefSeq v2.1 IDs\n")
cat("Example IDs look like: TraesCS2B03G1484200\n\n")

v2_col_number <- menu(
  choices = colnames(mapping),
  title = "Choose the RefSeq v2.1 gene-ID column"
)

v2_column <- colnames(mapping)[v2_col_number]

cat("\nSelected v2.1 column:", v2_column, "\n")


# ------------------------------------------------------------
# 6. Interactive selection of the older-ID column
# ------------------------------------------------------------

cat("\nSelect the column containing the OLDER RefSeq IDs\n")
cat("Example IDs look like: TraesCS2B02G592600\n\n")

old_col_number <- menu(
  choices = colnames(mapping),
  title = "Choose the older RefSeq gene-ID column"
)

old_column <- colnames(mapping)[old_col_number]

cat("\nSelected old-ID column:", old_column, "\n")


# ------------------------------------------------------------
# 7. Create a clean two-column correspondence table
# ------------------------------------------------------------

id_mapping <- mapping %>%
  
  select(
    all_of(v2_column),
    all_of(old_column)
  )

colnames(id_mapping) <- c(
  "GeneID",
  "RefSeq_v1_ID"
)


# ------------------------------------------------------------
# 8. Remove blank and duplicated mappings
# ------------------------------------------------------------

id_mapping <- id_mapping %>%
  
  filter(
    !is.na(GeneID),
    !is.na(RefSeq_v1_ID),
    GeneID != "",
    RefSeq_v1_ID != ""
  ) %>%
  
  distinct()


# ------------------------------------------------------------
# 9. Check whether one v2.1 ID maps to multiple old IDs
# ------------------------------------------------------------

multiple_mappings <- id_mapping %>%
  
  count(GeneID) %>%
  
  filter(n > 1)

cat("\n============================================\n")
cat("MULTIPLE MAPPINGS CHECK\n")
cat("============================================\n")

print(multiple_mappings)


# ------------------------------------------------------------
# 10. Join the old IDs to the 17-gene master table
# ------------------------------------------------------------

master_mapped <- master %>%
  
  select(
    -any_of("RefSeq_v1_ID")
  ) %>%
  
  left_join(
    id_mapping,
    by = "GeneID"
  )


# ------------------------------------------------------------
# 11. Check mapping success
# ------------------------------------------------------------

mapped_genes <- master_mapped %>%
  filter(!is.na(RefSeq_v1_ID))

unmapped_genes <- master_mapped %>%
  filter(is.na(RefSeq_v1_ID))


cat("\n============================================\n")
cat("STEP 5A MAPPING RESULTS\n")
cat("============================================\n")

cat(
  "Candidate genes:",
  nrow(master),
  "\n"
)

cat(
  "Successfully mapped:",
  nrow(mapped_genes),
  "\n"
)

cat(
  "Not mapped:",
  nrow(unmapped_genes),
  "\n"
)

cat("============================================\n")


# ------------------------------------------------------------
# 12. Display mapped IDs
# ------------------------------------------------------------

print(
  master_mapped %>%
    select(
      GeneID,
      RefSeq_v1_ID,
      Candidate_Route
    )
)

View(master_mapped)


# ------------------------------------------------------------
# 13. Show any unmapped genes separately
# ------------------------------------------------------------

if (nrow(unmapped_genes) > 0) {
  
  cat("\nUnmapped genes:\n")
  
  print(
    unmapped_genes$GeneID
  )
}


# ------------------------------------------------------------
# 14. Ask where to save result
# ------------------------------------------------------------

cat(
  "\nPlease choose a folder where STEP 5A results should be saved...\n"
)

output_folder <- choose.dir(
  caption =
    "Choose folder for Step 5A RefSeq ID mapping"
)


# ------------------------------------------------------------
# 15. Save mapped master table
# ------------------------------------------------------------

write_csv(
  master_mapped,
  file.path(
    output_folder,
    "TaDof_Step5A_master_candidates_with_old_IDs.csv"
  )
)

write_csv(
  unmapped_genes,
  file.path(
    output_folder,
    "TaDof_Step5A_unmapped_genes.csv"
  )
)


# ------------------------------------------------------------
# 16. Final message
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 5A COMPLETE\n")
cat("============================================\n")

cat(
  "Mapped:",
  nrow(mapped_genes),
  "of",
  nrow(master),
  "candidate genes\n"
)

cat(
  "\nResults saved in:\n",
  output_folder,
  "\n"
)

cat("============================================\n")
