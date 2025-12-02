library(ggplot2)
library(reshape2)

# Selecionar colunas LFQ
lfq_cols_qex <- grep("LFQ.intensity", colnames(mq_qexactive), value = TRUE)

# Extrair apenas LFQ
lfq_qex <- mq_qexactive[, lfq_cols_qex]

# Transformar log2 (adicionar 1 para evitar log(0))
lfq_log2 <- log2(lfq_qex + 1)

# Converter para formato longo para ggplot
lfq_melt <- melt(lfq_log2)

ggplot(lfq_melt, aes(value)) +
  geom_histogram(bins = 60, fill = "steelblue") +
  theme_minimal() +
  labs(title = "Distribution of LFQ Intensities (log2)",
       x = "log2(LFQ intensity)", y = "Frequency")
