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

# 1. Keep only essential columns
lfq_cols_qex <- grep("^LFQ.intensity", colnames(mq_qexactive), value = TRUE)
mq_qex_small <- mq_qexactive[, c("Protein.IDs", "Protein.names", "Gene.names", lfq_cols_qex)]

lfq_cols_orb <- grep("^LFQ.intensity", colnames(mq_orbitrap), value = TRUE)
mq_orbitrap_small <- mq_orbitrap[, c("Protein.IDs", "Protein.names", "Gene.names", lfq_cols_orb)]

# 2. Remove rows with multiple proteins (Protein.IDs contains ";"), but keep isoforms
remove_multi_proteins <- function(df) {
  df <- df[!grepl("^REV__", df$Protein.IDs), ]
  df <- df[!grepl("^CON__", df$Protein.IDs), ]
  df <- df[!grepl(";", df$Protein.IDs), ]  # eliminate multi-identifications
  return(df)
}

mq_qex_small <- remove_multi_proteins(mq_qex_small)
mq_orbitrap_small <- remove_multi_proteins(mq_orbitrap_small)

# 3. Remove rows with multiple gene names
mq_qex_small <- mq_qex_small[!grepl(";", mq_qex_small$Gene.names), ]
mq_orbitrap_small <- mq_orbitrap_small[!grepl(";", mq_orbitrap_small$Gene.names), ]

# 4. Remove rows with no LFQ signal
# mq_qex_small <- mq_qex_small[rowSums(mq_qex_small[, lfq_cols_qex] > 0, na.rm = TRUE) > 0, ]
# mq_orbitrap_small <- mq_orbitrap_small[rowSums(mq_orbitrap_small[, lfq_cols_orb] > 0, na.rm = TRUE) > 0, ]

# 5. Extract gene lists
genes_qexactive <- sort(unique(mq_qex_small$Gene.names))
genes_orbitrap  <- sort(unique(mq_orbitrap_small$Gene.names))
genes_combined  <- sort(unique(c(genes_qexactive, genes_orbitrap)))

# 6. Export
write.table(genes_qexactive, "gene_list_QExactive.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(genes_orbitrap, "gene_list_Orbitrap.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(genes_combined, "gene_list_combined.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("\nExported cleaned gene lists:\n")
cat(" - QExactive:", length(genes_qexactive), "genes\n")
cat(" - Orbitrap:", length(genes_orbitrap), "genes\n")
cat(" - Combined:", length(genes_combined), "genes\n")
