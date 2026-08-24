# ============================================================
# STEP 1: Identify all TaDof genes in IWGSC RefSeq v2.1
# using PF02701 and IPR003851
# ============================================================

# Load packages
library(readr)
library(dplyr)

# ------------------------------------------------------------
# 1. Ask user to select/upload the annotation file
# ------------------------------------------------------------

cat("Please select the IWGSC RefSeq v2.1 functional annotation CSV file...\n")

annotation_file <- file.choose()

cat("Selected file:\n")
print(annotation_file)

# ------------------------------------------------------------
# 2. Read the annotation file
# ------------------------------------------------------------

ann <- read_csv(annotation_file)

# Check column names
cat("\nColumn names detected:\n")
print(colnames(ann))

# ------------------------------------------------------------
# 3. Identify Dof genes using Pfam PF02701
# ------------------------------------------------------------

dof_pfam <- ann %>%
  filter(f.name == "PF02701") %>%
  distinct(g2.identifier) %>%
  arrange(g2.identifier)

cat("\nNumber of unique genes identified by PF02701:",
    nrow(dof_pfam), "\n")

# ------------------------------------------------------------
# 4. Identify Dof genes using InterPro IPR003851
# ------------------------------------------------------------

dof_interpro <- ann %>%
  filter(f.name == "IPR003851") %>%
  distinct(g2.identifier) %>%
  arrange(g2.identifier)

cat("Number of unique genes identified by IPR003851:",
    nrow(dof_interpro), "\n")

# ------------------------------------------------------------
# 5. Compare Pfam and InterPro gene sets
# ------------------------------------------------------------

pfam_only <- setdiff(
  dof_pfam$g2.identifier,
  dof_interpro$g2.identifier
)

interpro_only <- setdiff(
  dof_interpro$g2.identifier,
  dof_pfam$g2.identifier
)

cat("\nGenes found by PF02701 but not IPR003851:\n")
print(pfam_only)

cat("\nGenes found by IPR003851 but not PF02701:\n")
print(interpro_only)

# ------------------------------------------------------------
# 6. Create final unique TaDof gene list
# ------------------------------------------------------------

TaDof_genes <- ann %>%
  filter(f.name %in% c("PF02701", "IPR003851")) %>%
  distinct(g2.identifier) %>%
  arrange(g2.identifier)

cat("\nFinal number of unique TaDof genes:",
    nrow(TaDof_genes), "\n")

# ------------------------------------------------------------
# 7. Create summary table
# ------------------------------------------------------------

summary_table <- data.frame(
  Method = c(
    "Pfam PF02701",
    "InterPro IPR003851",
    "Final unique Dof genes"
  ),
  Number_of_genes = c(
    nrow(dof_pfam),
    nrow(dof_interpro),
    nrow(TaDof_genes)
  )
)

print(summary_table)

# ------------------------------------------------------------
# 8. View final Dof list
# ------------------------------------------------------------

View(TaDof_genes)

# ------------------------------------------------------------
# 9. Ask user where to save results
# ------------------------------------------------------------

output_folder <- choose.dir(
  caption = "Choose a folder where results should be saved"
)

# ------------------------------------------------------------
# 10. Save results
# ------------------------------------------------------------

write_csv(
  TaDof_genes,
  file.path(output_folder, "TaDof_100_genes_RefSeqv2.1.csv")
)

write_csv(
  summary_table,
  file.path(output_folder, "TaDof_domain_summary_RefSeqv2.1.csv")
)

write_csv(
  dof_pfam,
  file.path(output_folder, "TaDof_PF02701_genes.csv")
)

write_csv(
  dof_interpro,
  file.path(output_folder, "TaDof_IPR003851_genes.csv")
)

# ------------------------------------------------------------
# 11. Final confirmation
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 1 COMPLETE\n")
cat("============================================\n")
cat("PF02701 genes:", nrow(dof_pfam), "\n")
cat("IPR003851 genes:", nrow(dof_interpro), "\n")
cat("Final unique TaDof genes:", nrow(TaDof_genes), "\n")
cat("Results saved to:\n", output_folder, "\n")
cat("============================================\n")
