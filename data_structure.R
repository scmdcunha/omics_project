### --- Setup and Data Import --- ###

# Set working directory
setwd("~/omics_project/")

# Check working dir and files
cat("Working directory:", getwd(), "\n")
files_needed <- c("MaxQuant_proteinGroups_QExactivePlusData.txt",
                  "MaxQuant_proteinGroups_OrbitrapEliteData.txt")
for (f in files_needed) {
  if (!file.exists(f)) {
    stop(paste("File not found:", f, "\nPlease check the filename and working directory."))
  } else {
    cat("Found file:", f, "\n")
  }
}

# Read MaxQuant output files (fast and safe)
mq_qexactive <- read.delim("MaxQuant_proteinGroups_QExactivePlusData.txt",
                           sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mq_orbitrap <- read.delim("MaxQuant_proteinGroups_OrbitrapEliteData.txt",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Print basic structure
cat("\n=== Dataset dimensions ===\n")
cat("Q Exactive Plus dataset:\n")
print(dim(mq_qexactive))
cat("Orbitrap Elite dataset:\n")
print(dim(mq_orbitrap))

# Quick look at columns and first rows
cat("\n=== First 8 columns (Q Exactive Plus) ===\n")
print(colnames(mq_qexactive)[1:8])
print(head(mq_qexactive[, 1:8], n = 3))

cat("\n=== First 8 columns (Orbitrap) ===\n")
print(colnames(mq_orbitrap)[1:8])
print(head(mq_orbitrap[, 1:8], n = 3))

# Verify presence of key columns and show example names if present
cat("\n=== Checking for key columns ===\n")
print(grep("Gene.names", colnames(mq_qexactive), value = TRUE))
print(grep("LFQ.intensity", colnames(mq_qexactive), value = TRUE))
print(grep("Gene.names", colnames(mq_orbitrap), value = TRUE))
print(grep("LFQ.intensity", colnames(mq_orbitrap), value = TRUE))

### --- Data Cleaning and Gene Extraction --- ###

clean_and_extract_genes <- function(df) {
  
  # If columns exist, filter properly
  filter_cols <- c("Reverse", "Potential.contaminant", "Only.identified.by.site")
  have_cols <- filter_cols %in% colnames(df)
  
  if (all(have_cols)) {
    df_clean <- df[
      df$Reverse != "+" &
        df$Potential.contaminant != "+" &
        df$Only.identified.by.site != "+",
    ]
  } else {
    warning("Filter columns missing — using unfiltered data.")
    df_clean <- df
  }
  
  # Extract gene names
  if (!"Gene.names" %in% colnames(df_clean)) {
    stop("Column 'Gene.names' not found in the data.")
  }
  
  # separate genes like "TP53;MDM2"
  genes <- unique(na.omit(unlist(strsplit(df_clean$Gene.names, ";"))))
  
  return(genes)
}


# Apply to both datasets
genes_qexactive <- clean_and_extract_genes(mq_qexactive)
genes_orbitrap  <- clean_and_extract_genes(mq_orbitrap)

# Combine gene lists from both instruments
genes_combined <- sort(unique(c(genes_qexactive, genes_orbitrap)))

# Export lists for GO analysis
write.table(genes_qexactive, "gene_list_QExactive.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(genes_orbitrap, "gene_list_Orbitrap.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(genes_combined, "gene_list_combined.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("\nExported gene lists for GO analysis:\n")
cat(" - gene_list_QExactive.txt (", length(genes_qexactive), "genes)\n")
cat(" - gene_list_Orbitrap.txt (", length(genes_orbitrap), "genes)\n")
cat(" - gene_list_combined.txt (", length(genes_combined), "genes)\n")

