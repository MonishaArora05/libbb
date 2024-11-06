# Check if BiocManager is installed, and install if necessary
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install necessary Bioconductor packages
BiocManager::install("DESeq2")
BiocManager::install("pasilla")

# Load the libraries
library(DESeq2)
library(pasilla)

# Input for count data and metadata
pasCts <- system.file("extdata", "pasilla_gene_counts.tsv", package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata", "pasilla_sample_annotation.csv", package="pasilla", mustWork=TRUE)

# Read in the count matrix and sample metadata
count_csv <- read.csv(pasCts, sep="\t", row.names="gene_id")
cts <- as.matrix(count_csv)
coldata <- read.csv(pasAnno, row.names=1)

# Convert condition to a factor and adjust rownames
coldata$condition <- factor(coldata$condition)
rownames(coldata) <- sub("fb", "", rownames(coldata))

# Ensure data alignment
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)

# Perform differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# Save results to CSV
write.csv(res, file="DESeq_Analysis.csv")

# Read in results and label differentially expressed genes
df <- read.csv("./DESeq_Analysis.csv", header=TRUE)
df$Diffexpressed <- "NO"
df$Diffexpressed[df$log2FoldChange > 0.1 & df$pvalue < 0.05] <- "UP"
df$Diffexpressed[df$log2FoldChange < -0.1 & df$pvalue < 0.05] <- "DOWN"

# Load ggplot2 for visualization
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# Plotting results
ggplot(data = df, aes(x = log2FoldChange, y = -log10(pvalue), col = Diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = 'dashed') + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "green", linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.00003), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
  geom_point(size = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 10, colour = "black", face = "bold"),
        axis.title.x = element_text(size = 10, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 8, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 10, colour = "black", face="bold"),
        legend.text = element_text(size = 8, colour="black", face="bold")) +
  scale_color_manual(values = c("#00AFBB", "grey", "pink"),
                     labels = c("Downregulated", "Not Significant", "Upregulated"))

# Export upregulated and downregulated genes
upregulated_genes <- rownames(df[df$Diffexpressed == "UP", ])
downregulated_genes <- rownames(df[df$Diffexpressed == "DOWN", ])
write(upregulated_genes, file = "upregulated_genes.txt")
write(downregulated_genes, file = "downregulated_genes.txt")

# Reading additional biological process data (requires .txt files in the same directory)
biological_processes <- "./Biological_Processes.txt"
david_data <- read.delim(biological_processes, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
molecular_data <- read.delim("./Molecular_Functions.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
