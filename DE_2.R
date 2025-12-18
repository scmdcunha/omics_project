## ===============================
##  LFQ NORMALIZATION PIPELINE (Robust PCA)
## ===============================

#Authors: Diogo Belbute nº*, Luís Mota nº66176, Sara Cunha nº*

## ---- Libraries ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(limma)
library(openxlsx)
library(readr)

## --- Directory and OutputFolder ---
setwd("/Applications/Mestrado/Omics_Aproaches/")
output_folder <- "Normalization_2"
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

## ===============================
## 1. Normalization
## ===============================
normalize_lfq <- function(df_long, min_fraction = 0.7) {
  
  df_norm <- df_long %>%
    mutate(log2_LFQ = ifelse(LFQ_Intensity > 0, log2(LFQ_Intensity), NA)) %>%
    group_by(Tissue, Protein.IDs) %>%
    filter(sum(!is.na(log2_LFQ)) / n() >= min_fraction) %>%
    ungroup() %>%
    group_by(LFQ_column) %>%
    mutate(log2_LFQ_norm = log2_LFQ - median(log2_LFQ, na.rm = TRUE)) %>%
    ungroup()
  
  df_norm
}

## ===============================
## 2. Prepare limma input
## ===============================
prepare_limma <- function(df_norm) {
  
  expr <- df_norm %>%
    select(Protein.IDs, LFQ_column, log2_LFQ_norm) %>%
    pivot_wider(names_from = LFQ_column, values_from = log2_LFQ_norm) %>%
    column_to_rownames("Protein.IDs")
  
  samples <- df_norm %>%
    select(LFQ_column, Tissue) %>%
    distinct() %>%
    arrange(LFQ_column)
  
  expr <- expr[, samples$LFQ_column, drop = FALSE]
  
  design <- model.matrix(~ 0 + Tissue, data = samples)
  colnames(design) <- gsub("Tissue", "", colnames(design))
  
  list(expr = expr, design = design, samples = samples)
}

## ===============================
## 3. Fit limma model
## ===============================
fit_limma <- function(expr, design) {
  
  contrast_matrix <- makeContrasts(
    Fallopian_vs_HGSC = Fallopian - HGSC,
    Fallopian_vs_EC   = Fallopian - EC,
    HGSC_vs_EC        = HGSC - EC,
    levels = design
  )
  
  fit <- lmFit(expr, design)
  fit <- contrasts.fit(fit, contrast_matrix)
  fit <- eBayes(fit)
  
  fit
}

## ===============================
## 4. Boxplot QC
## ===============================
plot_boxplots <- function(df_norm, output_folder = "Normalization") {
  
  dir.create(output_folder, showWarnings = FALSE)
  
  p_raw <- ggplot(df_norm, aes(Tissue, log2_LFQ, colour = Tissue)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Boxplot of log2 LFQ Intensities")
  
  ggsave(file.path(output_folder, "Boxplot_log2_LFQ.png"),
         p_raw,
         width = 10, height = 6)
  
  p_norm <- ggplot(df_norm, aes(Tissue, log2_LFQ_norm, colour = Tissue)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Boxplot of normalized log2 LFQ Intensities")
  
  ggsave(file.path(output_folder, "Boxplot_log2_LFQ_norm.png"), p_norm, width = 10, height = 6)
  
  list(raw = p_raw, normalized = p_norm)
}

## ===============================
## 5. PCA (safe)
## ===============================
plot_pca <- function(expr, samples) {
  
  mat <- as.matrix(expr)
  storage.mode(mat) <- "numeric"
  
  mat[!is.finite(mat)] <- NA
  mat <- t(apply(mat, 1, function(x) {
    if (all(is.na(x))) return(rep(NA, length(x)))
    x[is.na(x)] <- min(x, na.rm = TRUE) - 1
    x
  }))
  
  pca <- prcomp(t(mat), scale. = FALSE)
  var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  
  df <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  ) %>% left_join(samples, by = c("Sample" = "LFQ_column"))
  
  ggplot(df, aes(PC1, PC2, colour = Tissue)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      title = "PCA of normalized LFQ data",
      x = paste0("PC1 (", var[1], "%)"),
      y = paste0("PC2 (", var[2], "%)")
    )
}

## ===============================
## 6. RUN PIPELINE
## ===============================
OT_Processed <- read.delim("Proteomics/OT_Processed.tsv")

df_norm     <- normalize_lfq(OT_Processed)
limma_input <- prepare_limma(df_norm)

plot_pca(limma_input$expr, limma_input$samples)
plot_boxplots(df_norm, output_folder)

## ===============================
## 7. limma analysis
## ===============================
fit <- fit_limma(limma_input$expr, limma_input$design)

metadata <- OT_Processed %>%
  select(Protein.IDs, Gene.names, Protein.names) %>%
  distinct()

## ===============================
## Detection table
## ===============================
detection_table <- df_norm %>%
  group_by(Protein.IDs, Tissue) %>%
  summarise(
    n_detected = sum(!is.na(log2_LFQ_norm)),
    mean_log2_LFQ = mean(log2_LFQ_norm, na.rm = TRUE),
    .groups = "drop"
  )

## ===============================
## Tissue-specific proteins
## ===============================
tissue_specific_proteins <- detection_table %>%
  filter(n_detected > 0) %>%
  group_by(Protein.IDs) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  left_join(metadata, by = "Protein.IDs")

write_tsv(
  tissue_specific_proteins,
  generate_unique_filename("Tissue_Specific_Proteins", ".tsv", output_folder)
)

write.xlsx(
  tissue_specific_proteins,
  generate_unique_filename("Tissue_Specific_Proteins", ".xlsx", output_folder)
)

## ===============================
## Contrast-aware filtering
## ===============================
filter_contrast_detected <- function(results, detection_table,
                                     tissue_A, tissue_B, min_reps = 2) {
  
  keep <- detection_table %>%
    filter(Tissue %in% c(tissue_A, tissue_B), n_detected >= min_reps) %>%
    group_by(Protein.IDs) %>%
    filter(n() == 2) %>%
    pull(Protein.IDs)
  
  results %>%
    filter(Protein.IDs %in% keep) %>%
    filter(!is.na(logFC))
}

## ===============================
## Differential expression results
## ===============================
results_F_vs_HGSC <- topTable(fit, "Fallopian_vs_HGSC", Inf) %>%
  rownames_to_column("Protein.IDs") %>%
  left_join(metadata, by = "Protein.IDs") %>%
  filter_contrast_detected(detection_table, "Fallopian", "HGSC")%>%
  mutate(
    hit = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Hit", "No Hit")
  )

results_F_vs_EC <- topTable(fit, "Fallopian_vs_EC", Inf) %>%
  rownames_to_column("Protein.IDs") %>%
  left_join(metadata, by = "Protein.IDs") %>%
  filter_contrast_detected(detection_table, "Fallopian", "EC")%>%
  mutate(
    hit = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Hit", "No Hit")
  )

results_HGSC_vs_EC <- topTable(fit, "HGSC_vs_EC", Inf) %>%
  rownames_to_column("Protein.IDs") %>%
  left_join(metadata, by = "Protein.IDs") %>%
  filter_contrast_detected(detection_table, "HGSC", "EC")%>%
  mutate(
    hit = ifelse( P.Value < 0.05 & abs(logFC) > 1, "Hit", "No Hit")
  )

## ===============================
## Save DE results
## ===============================
write_tsv(results_F_vs_HGSC,
          generate_unique_filename("Fallopian_HGSC", ".tsv", output_folder))
write.xlsx(results_F_vs_HGSC,
           generate_unique_filename("Fallopian_HGSC", ".xlsx", output_folder))

write_tsv(results_F_vs_EC,
          generate_unique_filename("Fallopian_EC", ".tsv", output_folder))
write.xlsx(results_F_vs_EC,
           generate_unique_filename("Fallopian_EC", ".xlsx", output_folder))

write_tsv(results_HGSC_vs_EC,
          generate_unique_filename("HGSC_EC", ".tsv", output_folder))
write.xlsx(results_HGSC_vs_EC,
           generate_unique_filename("HGSC_EC", ".xlsx", output_folder))

## ===============================
## Volcano plots (filtered DE proteins)
## ===============================

## Fallopian vs HGSC
volcano_F_vs_HGSC <- ggplot(
  results_F_vs_HGSC,
  aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(
    aes(color = adj.P.Val < 0.05 & abs(logFC) > 1),
    alpha = 0.6) +
  scale_color_manual(
    values = c("FALSE" = "grey70", "TRUE" = "red")) +
  theme_minimal() +
  labs(
    title = "Fallopian vs HGSC",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value",
    color = "Significant")

ggsave(
  filename = generate_unique_filename("F_HGSC_Volcano", ".png", output_folder),
  plot = volcano_F_vs_HGSC,
  width = 10,
  height = 8)

## Fallopian vs EC
volcano_F_vs_EC <- ggplot(
  results_F_vs_EC,
  aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(
    aes(color = adj.P.Val < 0.05 & abs(logFC) > 1),
    alpha = 0.6) +
  scale_color_manual(
    values = c("FALSE" = "grey70", "TRUE" = "red")) +
  theme_minimal() +
  labs(
    title = "Fallopian vs EC",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value",
    color = "Significant")

ggsave(
  filename = generate_unique_filename("F_EC_Volcano", ".png", output_folder),
  plot = volcano_F_vs_EC,
  width = 10,
  height = 8)

## HGSC vs EC
volcano_HGSC_vs_EC <- ggplot(
  results_HGSC_vs_EC,
  aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(
    aes(color = P.Value < 0.05 & abs(logFC) > 1),
    alpha = 0.6) +
  scale_color_manual(
    values = c("FALSE" = "grey70", "TRUE" = "red")) +
  theme_minimal() +
  labs(
    title = "HGSC vs EC",
    x = "log2 Fold Change",
    y = "-log10 p-value",
    color = "Significant")

ggsave(
  filename = generate_unique_filename("HGSC_EC_Volcano", ".png", output_folder),
  plot = volcano_HGSC_vs_EC,
  width = 10,
  height = 8)

