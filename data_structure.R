### --- Setup and Data Import --- ###

# Set working directory
setwd("~/omics_project/")

# Read MaxQuant output files
mq_qexactive <- read.delim("MaxQuant_proteinGroups_QExactivePlusData.txt", 
                           sep = "\t", header = TRUE)
mq_orbitrap <- read.delim("MaxQuant_proteinGroups_OrbitrapEliteData.txt", 
                          sep = "\t", header = TRUE)

# Check basic structure
cat("Q Exactive Plus dataset:\n")
dim(mq_qexactive)
cat("Orbitrap Elite dataset:\n")
dim(mq_orbitrap)

# Quick look at columns and first rows
head(mq_qexactive[, 1:8])
head(mq_orbitrap[, 1:8])

# Verify presence of key columns
grep("Gene.names", colnames(mq_qexactive), value = TRUE)
grep("LFQ.intensity", colnames(mq_qexactive), value = TRUE)
grep("Gene.names", colnames(mq_orbitrap), value = TRUE)
grep("LFQ.intensity", colnames(mq_orbitrap), value = TRUE)

### --- Data Cleaning and Gene Extraction --- ###

# Function to clean and extract valid gene names
clean_and_extract_genes <- function(df) {
  df_clean <- subset(df, 
                     is.na(Reverse) & 
                       is.na(Potential.contaminant) & 
                       is.na(Only.identified.by.site))
  genes <- unique(na.omit(df_clean$Gene.names))
  return(genes)
}

# Apply to both datasets
genes_qexactive <- clean_and_extract_genes(mq_qexactive)
genes_orbitrap  <- clean_and_extract_genes(mq_orbitrap)

# Combine gene lists from both instruments
genes_combined <- unique(c(genes_qexactive, genes_orbitrap))

# Export lists for GO analysis
write.table(genes_qexactive, "gene_list_QExactive.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(genes_orbitrap, "gene_list_Orbitrap.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(genes_combined, "gene_list_combined.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Exported gene lists for GO analysis.\n")
