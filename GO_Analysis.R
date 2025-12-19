## ===============================
## GO Functional Annotation
## ===============================

#Authors: Diogo Belbute nº66156, Luís Mota nº66176, Sara Cunha nº66039

## ---- Libraries ----
library(readr)
library(openxlsx)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

## --- Directory and OutputFolder ---
setwd("/Applications/Mestrado/Omics_Aproaches/") #Change if needed
output_folder <- "GO_Analysis" #Change if needed
dir.create(output_folder, showWarnings = FALSE)

##--Saving Files Function--
generate_unique_filename <- function(base_name, ext = ".tsv", folder = ".") {
  i <- 0
  repeat {
    fn <- if (i == 0) paste0(base_name, ext) else paste0(base_name, "_", i, ext)
    full <- file.path(folder, fn)
    if (!file.exists(full)) return(full)
    i <- i + 1
  }
}

## ---- Data Import ----
# Import gene list
df <- read_tsv("Normalization_2/Fallopian_EC.tsv") # Change accordingly
hit_genes <- df %>%
  filter(hit == "Hit") %>%
  pull(Gene.names)

universe_genes <- df$Gene.names

# Clean up (important)
universe_genes <- unique(na.omit(universe_genes))
hit_genes      <- unique(na.omit(hit_genes))

# Check length (Usually need 20 genes)
length(hit_genes)
length(universe_genes)

## ---- GO Enrichment Analysis ----
ego_hits <- enrichGO(
  gene          = hit_genes,
  universe      = universe_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 1,
  minGSSize     = 10,
  maxGSSize     = 300,
  readable      = TRUE
)

# Simplify GO terms to reduce redundancy
ego_hits_simplified <- simplify(
  ego_hits,
  cutoff = 0.8,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

## ---- Results Visualization ----
dotplot(
  ego_hits_simplified,
  x = "FoldEnrichment",
  color = "p.adjust",
  size = "Count",
  showCategory = 20) + 
  labs(title = "Fallopian_EC_GO")
ggsave(
  filename = generate_unique_filename("Fallopian_EC_GO_Dotplot", ".png", output_folder),
  width = 10, height = 8)

## ---- Save Results ----
GO_Results <- as.data.frame(ego_hits_simplified)
write_tsv(
  GO_Results,
  generate_unique_filename("Fallopian_HGSC_GO", ".tsv", output_folder))

write.xlsx(
  GO_Results,
  generate_unique_filename("Fallopian_HGSC_GO", ".xlsx", output_folder))

## --- Protein GO Annotation ---
go_res <- as.data.frame(ego_hits_simplified)
go_long <- go_res %>%
  dplyr::select(
    GO_ID   = ID,
    GO_Term = Description,
    p.adjust,
    FoldEnrichment,
    geneID
  ) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(Gene = geneID)

go_annotated <- go_long %>%
  dplyr::left_join(
    df %>% dplyr::select(Protein.IDs, Gene.names),
    by = c("Gene" = "Gene.names"))

GO_Protein_Table <- go_annotated %>%
  dplyr::select(
    Protein.IDs,
    Gene,
    GO_ID,
    GO_Term,
    FoldEnrichment,
    p.adjust
  ) %>%
  distinct()

## --- Save Annotation ---
write_tsv(GO_Protein_Table,
          generate_unique_filename("HGSC_EC_GO_Proteins", ".tsv", output_folder))
write.xlsx(GO_Protein_Table,
           generate_unique_filename("HGSC_EC_GO_Proteins", ".xlsx", output_folder))
