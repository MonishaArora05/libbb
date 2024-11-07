if(!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("pasilla")
BiocManager::install("ggplot2")
BiocManager::install("ggrepel")

library("DESeq2")
library("pasilla")
library("ggplot2")
library("ggrepel")

count_csv <- system.file("extdata", 
                         "pasilla_gene_counts.tsv", 
                         package = "pasilla", mustWork = TRUE)

annotation_csv <- system.file("extdata", 
                              "pasilla_sample_annotation.csv", 
                              package = "pasilla", mustWork = TRUE)

count_data <- read.csv(count_csv, sep = "\t", row.names = "gene_id")
count_matrix <- as.matrix(count_data)

annotation_data <- read.csv(annotation_csv, row.names = 1)
annotation_data$condition <- factor(annotation_data$condition)

rownames(annotation_data) <- sub("fb", "", rownames(annotation_data))

all(rownames(annotation_data) %in% colnames(count_matrix))

all(rownames(annotation_data) == colnames(count_matrix))

count_matrix <- count_matrix[, rownames(annotation_data)]

all(rownames(annotation_data) == colnames(count_matrix))

deseq <- DESeqDataSetFromMatrix(countData = count_matrix, 
                                colData = annotation_data, 
                                design = ~ condition)

diff_exp_analysis <- DESeq(deseq)

diff_exp_result <- results(diff_exp_analysis)

write.csv("./DESeq_Analysis.csv", header = TRUE)

dataframe <- read.csv("./DESeq_Analysis.csv", header = TRUE)

dataframe$expressed <- "NO"

dataframe$expressed[dataframe$log2FoldChange > 0.1 & dataframe$pvalue < 0.05] <- "UP"
dataframe$expressed[dataframe$log2FoldChange < -0.1 & dataframe$pvalue < 0.05] <- "DOWN"

upregulated_genes <- rownames(dataframe[dataframe$expressed == "UP", ])
downregulated_genes <- rownames(dataframe[dataframe$expressed == "DOWN", ])

write(upregulated_genes, file = "upregualated_genes.txt")
write(downregulated_genes, file = "downregulated_genes.txt")

ggplot(data = dataframe, aes(x = log2FoldChange, y = -log10(pvalue), col = expressed)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_point(size = 2) +
  theme_minimal() + 
  scale_color_manual(values = c("turquoise", "grey", "pink"), 
                     labels = c("Downregulated", "Not Significant", "Upregulated")) + 
  labs(x = "log2 Fold Change", y = "-log10 P-value", color = "Expression")
