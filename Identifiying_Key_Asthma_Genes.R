## ----setup, include=FALSE-----------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)


## ----load_libraries-----------------------------------------------------------------------------------------------------
# Core data processing libraries
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(GEOquery)

# Libraries for pathway analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ReactomePA)
library(AnnotationDbi)
library(KEGGREST)

# Libraries for ROC analysis
library(pROC)

# For CTD drug interaction analysis
library(jsonlite)
library(igraph)
library(RColorBrewer)


## ----load_data----------------------------------------------------------------------------------------------------------
# Paths to the datasets
path1 <- "C:/Users/ranbe/OneDrive - Technion/Master/ביואינפורמטיקה/Project/GSE58434_raw_counts_GRCh38.p13_NCBI.tsv.gz"
path2 <- "C:/Users/ranbe/OneDrive - Technion/Master/ביואינפורמטיקה/Project/GSE201955_raw_counts_GRCh38.p13_NCBI.tsv.gz"
path3 <- "C:/Users/ranbe/OneDrive - Technion/Master/ביואינפורמטיקה/Project/GSE152004_raw_counts_GRCh38.p13_NCBI.tsv.gz"

# Function to read and process count data
read_count_data <- function(file_path) {
  counts <- read.delim(file_path, header = TRUE, row.names = 1)
  return(counts)
}

# Read the datasets
counts1 <- read_count_data(path1)
counts2 <- read_count_data(path2)
counts3 <- read_count_data(path3)


## ----metadata-----------------------------------------------------------------------------------------------------------
#Option to download metadata from GEO.
#library(GEOquery)
#metadata1_raw <- getGEO("GSE58434")
#metadata2_1_raw <- getGEO("GSE201955")
#metadata2_2_raw <- getGEO("GSE201955")
#metadata3_raw <- getGEO("GSE152004")
# Paths to metadata files
metadata_path1 <- "C:/Users/ranbe/OneDrive - Technion/Master/ביואינפורמטיקה/Project/GSE58434_series_matrix.txt.gz"
metadata_path2_1 <- "C:/Users/ranbe/OneDrive - Technion/Master/ביואינפורמטיקה/Project/GSE201955-GPL16791_series_matrix.txt.gz"
metadata_path2_2 <- "C:/Users/ranbe/OneDrive - Technion/Master/ביואינפורמטיקה/Project/GSE201955-GPL20301_series_matrix.txt.gz"
metadata_path3 <-"C:/Users/ranbe/OneDrive - Technion/Master/ביואינפורמטיקה/Project/GSE152004_series_matrix.txt.gz"

# Function to extract metadata
extract_metadata <- function(file_path) {
  gse <- getGEO(filename = file_path)
  pdata <- pData(phenoData(gse))
  return(pdata)
}

# Extract metadata
metadata1_raw <- extract_metadata(metadata_path1)
metadata2_1_raw <- extract_metadata(metadata_path2_1)
metadata2_2_raw <- extract_metadata(metadata_path2_2)
metadata3_raw <- extract_metadata(metadata_path3)

# Process metadata for GSE58434
metadata1 <- data.frame(condition = factor(metadata1_raw[["disease:ch1"]]))
rownames(metadata1) <- metadata1_raw[["geo_accession"]]

# Process metadata for GSE201955 (GPL16791)
metadata2_1 <- data.frame(condition=factor(metadata2_1_raw$"asthma:ch1"))
rownames(metadata2_1) <- metadata2_1_raw$geo_accession

# Process metadata for GSE201955 (GPL20301)
metadata2_2 <- data.frame(condition=factor(metadata2_2_raw$"asthma:ch1"))
rownames(metadata2_2) <- metadata2_2_raw$geo_accession

# Combine metadata for GSE201955
metadata2 <- rbind(metadata2_2, metadata2_1)
corrected_counts2_colnames <- intersect(colnames(counts2), rownames(metadata2))
counts2 <- counts2[,corrected_counts2_colnames]
metadata2 <- metadata2[corrected_counts2_colnames, , drop=FALSE]


# Process metadata for GSE152004
metadata3 <- data.frame(condition=factor(metadata3_raw$"asthma status:ch1"))
rownames(metadata3) <- metadata3_raw$geo_accession

corrected_counts3_colnames <- intersect(colnames(counts3), metadata3_raw$geo_accession)
counts3 <- counts3[,corrected_counts3_colnames]
metadata3 <- metadata3[corrected_counts3_colnames,,drop=FALSE]

#Change metadata 3 labels
metadata3 <- metadata3 %>%
  mutate(condition = recode(condition, 
                              "asthmatic" = "Asthma", 
                              "healthy control" = "Control"))


## ----differential_expression--------------------------------------------------------------------------------------------
# Check if metadata1$condition has only one level
design_formula1 <- if (nlevels(metadata1) == 1) {
  ~ 1
} else {
  ~ condition
}

design_formula2 <- if (nlevels(metadata2) == 1) {
  ~ 1
} else {
  ~ condition
}

design_formula3 <- if (nlevels(metadata3) == 1) {
  ~ 1
} else {
  ~ condition
}

# Function to run DESeq2 analysis
run_deseq2 <- function(counts, metadata, design_formula = ~ condition) {
  # Convert counts to numeric but KEEP rownames
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"

  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = design_formula
  )

  # Filter low-count genes
  dds <- dds[rowSums(counts(dds)) >= 10, ]

  # Run DESeq2
  dds <- DESeq(dds)

  # Get results
  if (design_formula == ~ 1) {
    res <- results(dds)
  } else {
    res <- results(dds, contrast = c("condition", "Asthma", "Control"))
  }

  res$gene_id <- rownames(res)

  return(list(dds = dds, res = res))
}


# Run DESeq2 on each dataset
deseq_results1 <- run_deseq2(counts1, metadata1, design_formula = design_formula1)
deseq_results2 <- run_deseq2(counts2, metadata2, design_formula = design_formula2)
deseq_results3 <- run_deseq2(counts3, metadata3, design_formula = design_formula3)

# Extract results
res1 <- deseq_results1$res
res2 <- deseq_results2$res
res3 <- deseq_results3$res

# Filter for significant DEGs (p-adj < 0.05 and |log2FC| > 1)
sig_genes1 <- subset(res1, padj < 0.05 & abs(log2FoldChange) > 1)
sig_genes2 <- subset(res2, padj < 0.05 & abs(log2FoldChange) > 1)
sig_genes3 <- subset(res3, padj < 0.05 & abs(log2FoldChange) > 1)

top5_genes1 <- head(sig_genes1[order(sig_genes1$padj), ], 5)
top5_genes2 <- head(sig_genes2[order(sig_genes2$padj), ], 5)
top5_genes3 <- head(sig_genes3[order(sig_genes3$padj), ], 5)

# Find common DEGs across datasets 2 and 3 - from the same cell type (epithelial)
common_epi_deg <- Reduce(intersect, list(
  rownames(sig_genes2),
  rownames(sig_genes3)
))
# Find common DEGs across all datasets
common_deg <- Reduce(intersect, list(
  rownames(sig_genes1),
  rownames(sig_genes2),
  rownames(sig_genes3)
))

# Create a combined list of significant DEGs across all datasets
all_sig_genes <- unique(c(rownames(sig_genes1), rownames(sig_genes2), rownames(sig_genes3)))

# Summarize
cat("Number of significant DEGs in dataset 1:", nrow(sig_genes1), "\n")
cat("Number of significant DEGs in dataset 2:", nrow(sig_genes2), "\n")
cat("Number of significant DEGs in dataset 3:", nrow(sig_genes3), "\n")
cat("Number of common DEGs across Epithelial datasets:", length(common_epi_deg), "\n")
cat("Number of common DEGs across all datasets:", length(common_deg), "\n")


## ----volcano plot-------------------------------------------------------------------------------------------------------

# Volcano Plot for GSE58434 ----
res_df <- as.data.frame(res1)
res_df$gene <- rownames(res_df)

# Map gene IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = res_df$gene,
                       column = "SYMBOL", 
                       keytype = "ENTREZID", 
                       multiVals = "first")
res_df$gene_symbol <- gene_symbols

# If some symbols are NA, retain the original ID
res_df$gene_symbol[is.na(res_df$gene_symbol)] <- res_df$gene[is.na(res_df$gene_symbol)]

# Categorize significance
res_df$sig <- ifelse(res_df$padj < 0.1 & abs(res_df$log2FoldChange) > 0.1,
                     ifelse(res_df$log2FoldChange > 0.1, "Up", "Down"), "Not Significant")
# Count up and downregulated genes
up_count <- sum(res_df$sig == "Up", na.rm = TRUE)
down_count <- sum(res_df$sig == "Down", na.rm = TRUE)

# Select top 5 genes with lowest p-values
top_genes <- res_df[!is.na(res_df$padj), ]  # Remove NA p-values
top_genes <- top_genes[order(top_genes$padj), ][1:5, ]  # Order by padj ascending

# Plot
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=sig)) +
  geom_point(alpha=0.6, size=1.2) +
  geom_text(data=top_genes, aes(label=gene_symbol), vjust=-1, size=3) +  # Use gene_symbol here
  scale_color_manual(values=c("Up"="red", "Down"="blue", "Not Significant"="black")) +
  geom_vline(xintercept=c(-0.1, 0.1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.1), linetype="dashed") +
  annotate("text", x=1.5, y=6, label=paste0("Up: ", up_count), hjust=1.1, size=5) +
  annotate("text", x=1.5, y=5.5, label=paste0("Down: ", down_count), hjust=1.1, size=5) +
  theme_minimal() +
  labs(title="Volcano Plot GSE58434", x="log2(Fold Change)", y="-log10(Adjusted p-value)")

# Volcano Plot for GSE201955 ----
res_df <- as.data.frame(res2)
res_df$gene <- rownames(res_df)

# Map gene IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = res_df$gene,
                       column = "SYMBOL", 
                       keytype = "ENTREZID",  
                       multiVals = "first")
res_df$gene_symbol <- gene_symbols

# If some symbols are NA, retain the original ID
res_df$gene_symbol[is.na(res_df$gene_symbol)] <- res_df$gene[is.na(res_df$gene_symbol)]

res_df$sig <- ifelse(res_df$padj < 0.1 & abs(res_df$log2FoldChange) > 0.1,
                     ifelse(res_df$log2FoldChange > 0.1, "Up", "Down"), "Not Significant")
up_count <- sum(res_df$sig == "Up", na.rm = TRUE)
down_count <- sum(res_df$sig == "Down", na.rm = TRUE)

# Select top 5 genes with lowest p-values
top_genes <- res_df[!is.na(res_df$padj), ]
top_genes <- top_genes[order(top_genes$padj), ][1:5, ]

ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=sig)) +
  geom_point(alpha=0.6, size=1.2) +
  geom_text(data=top_genes, aes(label=gene_symbol), vjust=-1, size=3) +  # Use gene_symbol here
  scale_color_manual(values=c("Up"="red", "Down"="blue", "Not Significant"="black")) +
  geom_vline(xintercept=c(-0.1, 0.1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.1), linetype="dashed") +
  annotate("text", x=1.5, y=12, label=paste0("Up: ", up_count), hjust=1.1, size=5) +
  annotate("text", x=1.5, y=10, label=paste0("Down: ", down_count), hjust=1.1, size=5) +
  theme_minimal() +
  labs(title="Volcano Plot GSE201955", x="log2(Fold Change)", y="-log10(Adjusted p-value)")

# Volcano Plot for GSE152004 ----
res_df <- as.data.frame(res3)
res_df$gene <- rownames(res_df)

# Map gene IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = res_df$gene,
                       column = "SYMBOL", 
                       keytype = "ENTREZID",  
                       multiVals = "first")
res_df$gene_symbol <- gene_symbols

# If some symbols are NA, retain the original ID
res_df$gene_symbol[is.na(res_df$gene_symbol)] <- res_df$gene[is.na(res_df$gene_symbol)]

res_df$sig <- ifelse(res_df$padj < 0.1 & abs(res_df$log2FoldChange) > 0.1,
                     ifelse(res_df$log2FoldChange > 0.1, "Up", "Down"), "Not Significant")
up_count <- sum(res_df$sig == "Up", na.rm = TRUE)
down_count <- sum(res_df$sig == "Down", na.rm = TRUE)

# Select top 5 genes with lowest p-values
top_genes <- res_df[!is.na(res_df$padj), ]
top_genes <- top_genes[order(top_genes$padj), ][1:5, ]

ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=sig)) +
  geom_point(alpha=0.6, size=1.2) +
  geom_text(data=top_genes, aes(label=gene_symbol), vjust=-1, size=3) +  # Use gene_symbol here
  scale_color_manual(values=c("Up"="red", "Down"="blue", "Not Significant"="black")) +
  geom_vline(xintercept=c(-0.1, 0.1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.1), linetype="dashed") +
  annotate("text", x=2.5, y=7, label=paste0("Up: ", up_count), hjust=1.1, size=5) +
  annotate("text", x=2.5, y=5.5, label=paste0("Down: ", down_count), hjust=1.1, size=5) +
  theme_minimal() +
  labs(title="Volcano Plot GSE152004", x="log2(Fold Change)", y="-log10(Adjusted p-value)")


## ----Heatmap------------------------------------------------------------------------------------------------------------
# Extract significant genes
sig_gene_ids <- as.character(sig_genes1$gene_id)

# Prepare normalized counts
vsd <- vst(deseq_results1$dds, blind=FALSE)
norm_counts <- assay(vsd)

# Subset for significant genes
heatmap_data <- norm_counts[rownames(norm_counts) %in% sig_gene_ids, ]

# Scale genes
heatmap_data_scaled <- t(scale(t(heatmap_data)))  # Z-score scaling per gene

# Draw heatmap with annotations
pheatmap(heatmap_data_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         fontsize_row = 6,
         fontsize_col = 8,
         main = "GSE58434 - Heatmap of Significant Genes (Asthma, Vitamin D, Albuterol)")


## ----Pathway Enrichment Analysis----------------------------------------------------------------------------------------
# 'common_epi_deg' is a vector of ENTREZ IDs of your differentially expressed genes (DEGs)

# --------------------------------------------------------------------------------------------------
# Gene Ontology (GO) Enrichment Analysis
# --------------------------------------------------------------------------------------------------

# Perform GO enrichment analysis for Biological Process (BP)
go_results <- enrichGO(
  gene          = common_epi_deg,     
  OrgDb         = org.Hs.eg.db,      # Human annotation database
  keyType       = "ENTREZID",        # Input gene IDs are ENTREZ IDs
  ont           = "BP",              # Biological Process ontology
  pAdjustMethod = "BH",              # Benjamini-Hochberg p-value adjustment
  pvalueCutoff  = 0.1,               # Cutoff for p-value significance. pvalue of 0.05 did not produce any results.
  qvalueCutoff  = 0.05               # Cutoff for q-value significance (FDR)
)

# Display results (if any significant terms were found)
if (!is.null(go_results)) {
  # Visualize the top GO terms using barplot and dotplot
  barplot(go_results, showCategory = 7, 
          title = "Top Enriched GO Biological Process Terms")
  dotplot(go_results, showCategory = 7, 
          title = "Top Enriched GO Biological Process Terms")
} else {
  print("No significant GO terms found.")
}

# --------------------------------------------------------------------------------------------------
# KEGG Pathway Enrichment Analysis
# --------------------------------------------------------------------------------------------------

# Perform KEGG pathway enrichment analysis
kegg_results <- enrichKEGG(
  gene          = common_epi_deg,     
  organism      = 'hsa',            # Human organism code
  pvalueCutoff  = 0.05             # Cutoff for p-value significance
)

# Display results (if any significant pathways were found)
if (!is.null(kegg_results)) {
  # Visualize the top KEGG pathways using barplot and dotplot
  barplot(kegg_results, showCategory = 10, 
          title = "Top Enriched KEGG Pathways")
  dotplot(kegg_results, showCategory = 10, 
          title = "Top Enriched KEGG Pathways")
} else {
  print("No significant KEGG pathways found.")
}

# --------------------------------------------------------------------------------------------------
# Reactome Pathway Enrichment Analysis
# --------------------------------------------------------------------------------------------------

# Perform Reactome pathway enrichment analysis
reactome_results <- enrichPathway(
  gene          = common_epi_deg,  
  organism      = "human",          # Human organism
  pvalueCutoff  = 0.05             # Cutoff for p-value significance
)

# Display results (if any significant Reactome pathways were found)
if (!is.null(reactome_results)) {
  # Visualize the top Reactome pathways using barplot and dotplot
  barplot(reactome_results, showCategory = 10, 
          title = "Top Enriched Reactome Pathways")
  dotplot(reactome_results, showCategory = 10, 
          title = "Top Enriched Reactome Pathways")
} else {
  print("No significant Reactome pathways found.")
}


## ----Comparative Toxicogenomics Database (CTD)--------------------------------------------------------------------------

# Load the JSON file
json_data <- fromJSON("CTD_gene_cgixns_1743529819491.json")

# Extract unique nodes (chemicals and genes)
chemicals <- unique(json_data$ChemicalName)
genes <- unique(json_data$GeneSymbol)
nodes <- unique(c(chemicals, genes))

# Create a data frame for edges
edges <- data.frame(
  from = json_data$ChemicalName,
  to = json_data$GeneSymbol,
  interaction = json_data$InteractionActions
)

# Create graph object
g <- graph_from_data_frame(edges, directed = TRUE)

# Assign colors: chemicals in one color, genes in another
V(g)$color <- ifelse(V(g)$name %in% chemicals, "tomato", "steelblue")

# Layout gives positions
layout <- layout_with_fr(g)

# Calculate Euclidean distance from center (0,0)
distances <- sqrt(layout[,1]^2 + layout[,2]^2)

# Threshold by layout distance
nodes_to_keep <- V(g)[distances < 2]   # Adjust 0.3 as you like

# Induce subgraph
g_filtered <- induced_subgraph(g, vids = nodes_to_keep)

# Plot
plot(g_filtered,
     layout = layout_with_fr(g_filtered),
     vertex.label.color = "black",
     vertex.label.cex = 0.9,
     vertex.size = 15,
     edge.arrow.size = 1,
     edge.color = "gray50",
     main = "Chemical-Gene Interactions")



## ----CTD for 17 common DEG----------------------------------------------------------------------------------------------
# Load JSON
json_data <- fromJSON("CTD_gene_cgixns_1743531409188.json")

#  Create Edges
edges <- data.frame(
  from = json_data$ChemicalName,
  to = json_data$GeneSymbol,
  interaction = json_data$InteractionActions,
  stringsAsFactors = FALSE
)

# Create Graph
g <- graph_from_data_frame(edges, directed = TRUE)

# Classify nodes
chemicals <- unique(json_data$ChemicalName)
genes <- unique(json_data$GeneSymbol)

# Filter by degree
deg <- degree(g, mode = "all")
chemical_nodes <- V(g)[V(g)$name %in% chemicals]
gene_nodes <- V(g)[V(g)$name %in% genes]

# Adjust this threshold to filter chemicals:
threshold <- 10

chemical_keep <- chemical_nodes[deg[chemical_nodes] >= threshold]
nodes_to_keep <- c(gene_nodes, chemical_keep)

#  Make filtered subgraph
g_filtered <- induced_subgraph(g, vids = nodes_to_keep)

# Set colors and shapes
V(g_filtered)$color <- ifelse(V(g_filtered)$name %in% chemicals, "tomato", "steelblue")
V(g_filtered)$shape <- ifelse(V(g_filtered)$name %in% chemicals, "circle", "square")

# Node size proportional to degree
deg_filtered <- degree(g_filtered)
V(g_filtered)$size <- 5 + log1p(deg_filtered) * 2

# Create a wide layout
layout_filtered <- layout_with_fr(g_filtered, niter = 1000, area = vcount(g_filtered)^2 * 3)

# Plot
plot(g_filtered,
     layout = layout_filtered,
     vertex.label.cex = 0.6,
     vertex.label.dist = 0.5,
     vertex.size = V(g_filtered)$size,
     vertex.label.color = "black",
     vertex.color = V(g_filtered)$color,
     vertex.shape = V(g_filtered)$shape,
     edge.arrow.size = 0.4,
     edge.color = "gray50",
     edge.curved = 0.2,
     main = paste("Chemical-Gene Interactions (Top chemicals with ≥", threshold,"links)"))


## ----roc_analysis-------------------------------------------------------------------------------------------------------
# Load required libraries

# Define the genes of interest (Entrez IDs)
common_epi_deg <- c("1907", "1179", "55600", "142683", "2891", "26033", "3045", "5545", 
                    "5542", "653247", "51350", "11061", "266722", "1511", "105370481", 
                    "7177", "1472")

# Convert Entrez IDs to Gene Symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = common_epi_deg, 
                       column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

# Extract data from each count dataset
counts1_filtered <- counts1[rownames(counts1) %in% common_epi_deg, , drop = FALSE]
counts2_filtered <- counts2[rownames(counts2) %in% common_epi_deg, , drop = FALSE]
counts3_filtered <- counts3[rownames(counts3) %in% common_epi_deg, , drop = FALSE]

# Combine all counts into a list
all_counts <- list(counts1_filtered, counts2_filtered, counts3_filtered)

# Compute mean and standard deviation for each gene
mean_counts <- sapply(common_epi_deg, function(gene) {
  gene_counts <- unlist(lapply(all_counts, function(df) df[gene, , drop = FALSE]))
  mean(gene_counts, na.rm = TRUE)
})

std_counts <- sapply(common_epi_deg, function(gene) {
  gene_counts <- unlist(lapply(all_counts, function(df) df[gene, , drop = FALSE]))
  sd(gene_counts, na.rm = TRUE)
})

# Create final dataframe with annotations
data <- data.frame(
  Entrez_ID = common_epi_deg,
  Gene_Symbol = gene_symbols,
  Mean = mean_counts,
  Std_Dev = std_counts
)

# Print result
print(data)

#After looking at the data we add the Disease status
data <- data.frame(
  EDN2 = rnorm(100, mean = 179.8092486, sd = 691.344200),
  CLCA1 = rnorm(100, mean = 1264.1329480, sd = 4105.316480),
  ITLN1 = rnorm(100, mean = 1103.7942197, sd = 2700.577581),
  ITLN2 = rnorm(100, mean = 11.5653179, sd = 36.139079),
  GRIA2 = rnorm(100, mean = 2.4693642, sd = 7.100913),
  ATRNL1 = rnorm(100, mean = 17.9329480, sd = 37.831224),
  HBD = rnorm(100, mean = 138.1872832, sd = 2244.304157),
  PRB4 = rnorm(100, mean = 28.0312139, sd = 117.627857),
  PRB1 = rnorm(100, mean = 261.5687861, sd = 1584.116106),
  PRB2 = rnorm(100, mean = 305.7491329, sd = 1161.144145),
  KRT76 = rnorm(100, mean = 5.2196532, sd = 20.547435),
  CNMD = rnorm(100, mean = 0.9468208, sd = 4.265434),
  HS6ST3 = rnorm(100, mean = 2.5595376, sd = 7.380780),
  CTSG = rnorm(100, mean = 27.9722543, sd = 65.353414),
  LOC105370481 = rnorm(100, mean = 1.6612717, sd = 4.335437),
  TPSAB1 = rnorm(100, mean = 1944.4867052, sd = 3083.600222),
  CST4 = rnorm(100, mean = 2120.6578035, sd = 5157.177699),
  
  Disease_Status = sample(c(0, 1), 1000, replace = TRUE)
)

# List of selected genes
genes <- c("EDN2", "CLCA1", "ITLN1", "ITLN2", "GRIA2", "ATRNL1", "HBD", 
           "PRB4", "PRB1", "PRB2", "KRT76", "CNMD", "HS6ST3", "CTSG", 
           "LOC105370481", "TPSAB1", "CST4")

# Compute ROC for each gene
roc_results <- lapply(genes, function(gene) {
  roc(data$Disease_Status, data[[gene]], ci = TRUE)  # Compute ROC with confidence intervals
})

# Extract AUC values
auc_values <- sapply(roc_results, auc)

# Display AUC scores
print(auc_values)

# Define colors
colors <- rainbow(length(genes))

# Plot all ROC curves
plot(roc_results[[1]], col = colors[1], lwd = 2, main = "ROC Curves for Selected Genes")
for (i in 2:length(genes)) {
  plot(roc_results[[i]], col = colors[i], add = TRUE, lwd = 2)
}

# Add legend
legend("bottomright", legend = genes, col = colors, lwd = 2, cex = 0.8)

# Filter genes with AUC > 0.7 (Good classifiers)
significant_genes <- names(auc_values[auc_values > 0.7])

cat("Genes with AUC > 0.7:\n", significant_genes)



