#Omics Approaches

#Libraries
library(tidyverse)
library(dplyr)
library(reshape2)
library(openxlsx)

#Working Directory
setwd("/Applications/Mestrado/Omics_Aproaches/")

# Output Folder
output_folder <- "Proteomics"
dir.create(output_folder, showWarnings = FALSE)

# Utility Function
generate_unique_filename <- function(base_name, ext = ".tsv", folder = ".") {
  i <- 0
  repeat {
    fn <- if (i == 0) paste0(base_name, ext) else paste0(base_name, "_", i, ext)
    full <- file.path(folder, fn)
    if (!file.exists(full)) return(full)
    i <- i + 1
  }
}

#Data
OT_Data_raw <- read.delim("MaxQuant_proteinGroups_OrbitrapEliteData.txt", sep = "\t", dec = ".") %>%
  filter(!grepl("^REV_", Protein.IDs))
QEP_Data_raw <- read.delim("MaxQuant_proteinGroups_QExactivePlusData.txt", sep = "\t", dec = ".") %>%
  filter(!grepl("^REV_", Protein.IDs))

#Data Transformation Function
transform_data <- function(df) {
  
  library(dplyr)
  library(tidyr)
  
  # 1) Identify LFQ columns
  lfq_cols <- grep("^LFQ", colnames(df), value = TRUE)
  
  # 2) Remove all 'b' columns upfront
  lfq_cols <- lfq_cols[!grepl("b$", lfq_cols)]
  
  # 3) Extract T-number and priority (.2 > base)
  info <- tibble(
    lfq_col  = lfq_cols,
    # Extract numeric part of T-number
    T_number = as.numeric(stringr::str_extract(lfq_cols, "(?<=T|t)\\d+")),
    # Priority: .2 columns > base
    priority = ifelse(grepl("\\.2$", lfq_cols), 2, 1)
  )
  
  # 4) Keep only ONE column per T-number, preferring .2
  info <- info %>%
    group_by(T_number) %>%
    slice_max(priority, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(T_number)
  
  # 5) Assign tissue types based on T_number
  info <- info %>%
    mutate(Tissue = case_when(
      T_number %in% 1:4   ~ "Fallopian",
      T_number %in% 5:14  ~ "HGSC",
      T_number %in% 15:24 ~ "EC",
      TRUE ~ NA_character_
    ))
  
  # 6) Pivot data to long format
  df_long <- df %>%
    select(Protein.IDs, Protein.names, Gene.names, all_of(info$lfq_col)) %>%
    pivot_longer(
      cols = all_of(info$lfq_col),
      names_to = "LFQ_column",
      values_to = "LFQ_Intensity"
    ) %>%
    left_join(info, by = c("LFQ_column" = "lfq_col")) %>%
    arrange(T_number) %>%
    select(-priority)
  
  return(df_long)
}

#Function to calculate summary statistics
calculate_summary_stats <- function(df_long) {
  df_long %>%
    group_by(Protein.IDs, Protein.names, Gene.names, Tissue) %>%
    summarise(
      mean_LFQ = mean(LFQ_Intensity, na.rm = TRUE),
      median_LFQ = median(LFQ_Intensity, na.rm = TRUE),
      n_samples = sum(!is.na(LFQ_Intensity)),
      .groups = "drop"
    )
}

# Transform
OT_Data_long <- transform_data(OT_Data_raw)
QEP_Data_long <- transform_data(QEP_Data_raw)

# Calculate Summary Statistics
OT_Summary <- calculate_summary_stats(OT_Data_long)
QEP_Summary <- calculate_summary_stats(QEP_Data_long)

#Export Results
write_tsv(OT_Data_long, file = generate_unique_filename("OT_Processed", ".tsv", output_folder))
write_tsv(QEP_Data_long, file = generate_unique_filename("QEP_Processed", ".tsv", output_folder))

write.xlsx(OT_Data_long, file = generate_unique_filename("OT_Processed", ".xlsx", output_folder))
write.xlsx(QEP_Data_long, file = generate_unique_filename("QEP_Processed", ".xlsx", output_folder))

write_tsv(OT_Summary, file = generate_unique_filename("OT_Summary_Statistics", ".tsv", output_folder))
write_tsv(QEP_Summary, file = generate_unique_filename("QEP_Summary_Statistics", ".tsv", output_folder))
