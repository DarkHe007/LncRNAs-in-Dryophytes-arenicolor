# Load required packages, install if missing
if (!require("WGCNA")) {
  if (!require("BiocManager")) install.packages("BiocManager")
  BiocManager::install("WGCNA")  # WGCNA is available through Bioconductor
}
if (!require("tidyverse")) install.packages("tidyverse")  # Data manipulation and visualization
if (!require("doParallel")) install.packages("doParallel")  # Parallel processing

library(WGCNA)
library(tidyverse)
library(doParallel)

# Enable multi-threading for WGCNA functions to speed up computations
allowWGCNAThreads()

# DATA PREPROCESSING AND FILTERING

# Filter out lowly expressed genes (average expression > 1)
expr_filtered <- all[rowMeans(all) > 1, ]

# Filter genes by variance - remove bottom 20% least variable genes
gene_variance <- apply(expr_filtered, 1, var)
expr_filtered <- expr_filtered[gene_variance > quantile(gene_variance, 0.2), ]

# Further filter to keep only the top 75% most variable genes
gene_variance <- apply(expr_filtered, 1, var)  # Recompute variance after initial filtering
top_genes <- order(gene_variance, decreasing = TRUE)[1:round(nrow(expr_filtered) * 0.75)]
expr_filtered <- expr_filtered[top_genes, ]

# Save filtered expression matrix for future use
write.table(expr_filtered, "filtered_expression_matrix.tsv", sep="\t", quote=FALSE)

# WGCNA DATA PREPARATION

# Transpose expression matrix (WGCNA expects genes as columns, samples as rows)
datExpr <- t(expr_filtered)

# Check data for excessive missing values and genes with zero variance
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  # Remove problematic samples and genes if any were identified
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# SAMPLE QUALITY CONTROL

# Cluster samples to identify outliers
sampleTree <- hclust(dist(datExpr), method = "average")

# Visualize sample clustering to detect potential outliers
png("Sample_clustering.png", width = 2000, height = 1500, res = 200)
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", 
     cex.lab = 1.5, cex.axis = 1.5)
dev.off()

# NETWORK CONSTRUCTION PARAMETER SELECTION

# Test different soft-thresholding powers to find optimal value
powers <- 1:20  # Range of powers to test
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot scale independence and mean connectivity to help choose soft threshold
par(mfrow = c(1, 2))
# Scale independence plot - look for where curve flattens out
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2], type="b",
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit", 
     main="Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2], 
     labels=powers, col="red")

# Mean connectivity plot - check that connectivity doesn't drop too low
plot(sft$fitIndices[,1], sft$fitIndices[,5], type="b",
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", main="Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")

ggsave(filename="Soft_Power.png", height=6, width=12, units="in", dpi=300)

# Choose soft power - use estimated power or default to 6 if no clear estimate
softPower <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)

# NETWORK CONSTRUCTION

# Calculate adjacency matrix using chosen soft threshold
adjacency <- adjacency(datExpr, power = softPower)

# Transform adjacency into Topological Overlap Matrix (TOM)
TOM <- TOMsimilarity(adjacency)

# Convert TOM to dissimilarity (distance measure) for clustering
dissTOM <- 1 - TOM

# Save TOM for downstream analysis
write.csv(TOM, "TOM.csv", row.names = TRUE)

# GENE CLUSTERING AND MODULE IDENTIFICATION

# Cluster genes based on TOM dissimilarity
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Plot gene dendrogram
png("Gene_Dendrogram.png", width = 2000, height = 1500, res = 200)
plot(geneTree, main = "Gene Clustering on TOM-based Dissimilarity", 
     sub = "", xlab = "")
dev.off()

# Identify gene modules using dynamic tree cutting
dynamicMods <- cutreeDynamic(
  dendro = geneTree, distM = dissTOM,
  deepSplit = 2,                    # More sensitive splitting (0-4 scale)
  pamRespectsDendro = FALSE,        # Don't respect dendrogram for PAM stage
  minClusterSize = 30               # Minimum genes per module
)

# Convert module numbers to colors for visualization
moduleColors <- labels2colors(dynamicMods)

# Plot dendrogram with module colors
png("Dynamic_TreeCut.png", width = 2000, height = 1500, res = 200)
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# MODULE REFINEMENT

# Set height cutoff for merging similar modules
mergeCutHeight <- 0.25  # Modules with correlation > 0.75 will be merged

# Calculate module eigengenes (first principal components)
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes

# Cluster modules based on eigengene correlations
dissME <- 1 - cor(MEs)
hclustME <- hclust(as.dist(dissME), method = "average")

# Plot module clustering with merge cutoff line
png("Module_Eigengene_Clustering.png", width = 2000, height = 1500, res = 200)
plot(hclustME, main = "Module Eigengene Clustering", labels = FALSE)
abline(h = mergeCutHeight, col = "red")  # Line shows merge threshold
dev.off()

# Merge closely related modules
mergedModules <- mergeCloseModules(datExpr, moduleColors, 
                                   cutHeight = mergeCutHeight, verbose = 3)
moduleColors <- mergedModules$colors    # Updated module assignments
MEs <- mergedModules$newMEs             # Updated module eigengenes

# FINAL OUTPUT

# Create and save final module assignments
modules <- data.frame(Gene = colnames(datExpr), Module = moduleColors)
write.table(modules, "WGCNA_modules.tsv", sep="\t", quote=FALSE, row.names=FALSE)