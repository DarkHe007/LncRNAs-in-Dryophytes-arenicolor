# Load required packages for analysis
library(WGCNA)
library(tidyverse)
library(clusterProfiler)  # For enrichment analysis
library(org.Hs.eg.db)     # Human genome annotation database
library(DOSE)             # Disease Ontology Semantic and Enrichment analysis
library(biomaRt)          # Interface to BioMart databases
library(ggvenn)           # For creating Venn diagrams

# FILTER MODULES CONTAINING lncRNAs

# Read lncRNA gene names from file
lncrnas <- readLines("names.txt")

# Filter WGCNA modules to only include those containing lncRNAs
lncrnas_modules <- modules[modules$Gene %in% lncrnas, ]

print(lncrnas_modules)

# Count number of unique modules containing lncRNAs
length(table(lncrnas_modules$Module))

# Save lncRNA module assignments
write.table(file = "lncRNA_modules.csv", x = lncrnas_modules, col.names = F, quote = F, row.names = F)

# Extract modules that contain lncRNAs
selected_modules <- lncrnas_modules$Module

# Get all genes in modules that contain lncRNAs
genes_in_selected_modules <- modules[modules$Module %in% selected_modules, ]

# Save for downstream analysis
write.table(genes_in_selected_modules, "selected_modules_genes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

print(head(genes_in_selected_modules))

# MODULE-TRAIT RELATIONSHIPS ANALYSIS

# Define module eigengenes of interest (modules containing lncRNAs)
selected_ME <- c("MEblue", "MEcyan", "MEdarkgrey", "MEdarkorange", "MEdarkorange2",
                  "MEdarkred", "MEgreen", "MElightcoral", "MElightyellow", "MEmagenta",
                  "MEmagenta3", "MEorange", "MEpink", "MEplum1", "MEsalmon", "MEtan",
                  "MEturquoise", "MEviolet")

# Subset module eigengenes matrix to only include selected modules
subset_MEs <- MEs[, colnames(MEs) %in% selected_ME]

print(subset_MEs)

# Calculate correlation between module eigengenes and traits
selected_moduleTraitCor <- cor(subset_MEs, datTraits, use = "p")

# Create heatmap showing relationships between lncRNA-containing modules and traits
labeledHeatmap(
  Matrix = selected_moduleTraitCor,
  xLabels = colnames(datTraits),
  yLabels = selected_ME,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),  # Color scale from blue (negative) to red (positive)
  textMatrix = signif(selected_moduleTraitCor, 2),  # Display correlations with 2 significant digits
  main = "lncRNA Module-Trait Relationships"
)

# LINEAR MODEL ANALYSIS: MODULE EIGENGENES VS DEVELOPMENTAL STAGES

# Fit linear models for each module eigengene against stages
results <- lapply(colnames(MEs), function(module) {
  # Fit linear model: stages ~ module eigengene
  model <- lm(stages ~ MEs[[module]])
  
  # Extract coefficients and p-values
  summary_model <- summary(model)
  coef <- summary_model$coefficients[2, 1]  # Slope (effect of module on stage)
  p_value <- summary_model$coefficients[2, 4]  # P-value for the slope
  
  # Return results as a data frame
  data.frame(Module = module, Coefficient = coef, P_value = p_value)
})

# Combine results from all modules into single data frame
results_df <- do.call(rbind, results)

print(results_df)

# Create bar plot of module-stage relationships
ggplot(results_df, aes(x = Module, y = Coefficient, fill = P_value < 0.05)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Module-Stage Relationships",
       x = "Module",
       y = "Coefficient",
       fill = "Significant (p < 0.05)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Create dot plot alternative visualization
ggplot(results_df, aes(x = Module, y = Coefficient, color = P_value < 0.05)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Module-Stage Relationships (Linear Model)",
       x = "Module",
       y = "Coefficient",
       color = "Significant (p < 0.05)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

############################################################################
# PROTEIN MAPPING AND FUNCTIONAL ENRICHMENT ANALYSIS
############################################################################

# Read BLAST results against human reference database
blast_results <- read.csv("Referencia_uniprot.tsv", 
                          header = F, 
                          sep = "\t")

# Set column names for BLAST results
colnames(blast_results) <- c("query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits")

# Extract unique UniProt IDs from BLAST results
uniprot_ids <- unique(blast_results$target)

# Filter BLAST results for high-confidence hits
blast_filtered <- blast_results %>% 
  filter(fident > 0.9) %>%        # Minimum 90% identity
  filter(evalue < 0.001) %>%      # Maximum E-value threshold
  dplyr::select(1,2) %>%          # Keep only query and target columns
  distinct()                      # Remove duplicates

# Get target modules (modules containing lncRNAs)
targetModules <- unique(lncrnas_modules$Module)

# Extract all genes in the target modules
genesInModules <- modules$Gene[modules$Module %in% targetModules]

head(genesInModules)

# Save gene list for external analysis
write.table(genesInModules, file="selected_modules_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Create list structure storing genes for each module
module_genes_list <- list()

# Populate list with genes for each module
for (module in targetModules) {
  genes <- modules$Gene[modules$Module == module]
  module_genes_list[[module]] <- genes  # Store genes in the list
}

# BIO-MART ANNOTATION: MAP UNIPROT IDs TO GENE SYMBOLS

# Connect to Ensembl BioMart database
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Retrieve gene annotations for UniProt IDs
gene_annotations <- getBM(
  attributes = c("uniprotswissprot", "external_gene_name", "entrezgene_id"),
  filters = "uniprotswissprot",
  values = uniprot_ids,
  mart = ensembl
)

# MAP BLAST TARGETS TO MODULES

# Initialize list to store protein targets for each module
module_targets_list <- list()

# For each module, find protein targets of its genes
for (module in names(module_genes_list)) {
  genes <- module_genes_list[[module]]
  targets <- blast_filtered$target[blast_filtered$query %in% genes]
  module_targets_list[[module]] <- targets
}

# MAP ANNOTATIONS TO MODULES

# Initialize list to store gene annotations for each module
module_mappings_list <- list()

# For each module, get annotations for its protein targets
for (module in names(module_targets_list)) {
  # Get UniProt IDs for current module
  uniprot_ids <- module_targets_list[[module]]
  
  # Subset annotations for module's UniProt IDs
  module_mappings <- subset(gene_annotations, uniprotswissprot %in% uniprot_ids)
  
  # Store in list
  module_mappings_list[[module]] <- module_mappings
}

print(module_mappings_list)

# ENHANCE ANNOTATIONS WITH BLAST INFORMATION

# Merge annotation data with BLAST results to include query information
for (module in names(module_mappings_list)) {
  module_mappings <- module_mappings_list[[module]]
  
  # Merge annotations with BLAST results
  merged_data <- merge(module_mappings, blast_filtered, by.x = "uniprotswissprot", by.y = "target", all.x = TRUE)
  
  module_mappings_list[[module]] <- merged_data
}

print(module_mappings_list)

# Combine all module mappings into single data frame
all_mappings <- do.call(rbind, module_mappings_list)
print(all_mappings)

###########################################################################
# FUNCTIONAL ENRICHMENT ANALYSIS (GO AND KEGG)
###########################################################################

# Function to perform GO and KEGG enrichment analysis for a module
run_enrichment <- function(entrez_ids, module_name) {
  # Remove NA values from Entrez IDs
  entrez_ids <- na.omit(entrez_ids)
  
  # GO Enrichment Analysis - Biological Process, Cellular Component, Molecular Function
  go_results <- enrichGO(gene         = entrez_ids,
                         OrgDb        = org.Hs.eg.db,
                         keyType      = "ENTREZID",
                         ont          = "ALL",   # All three GO categories
                         pAdjustMethod = "bonferroni",  # Multiple testing correction
                         pvalueCutoff  = 0.01,   # Significance threshold
                         readable      = TRUE)   # Convert Entrez IDs to gene symbols
  
  # KEGG Pathway Enrichment Analysis
  kegg_results <- enrichKEGG(gene          = entrez_ids,
                             organism      = "hsa",  # Homo sapiens
                             pAdjustMethod = "bonferroni",
                             pvalueCutoff  = 0.01)
  
  # Save results to files
  write.table(as.data.frame(go_results), file=paste0(module_name, "_GO_enrichment.txt"), sep="\t", row.names=FALSE)
  write.table(as.data.frame(kegg_results), file=paste0(module_name, "_KEGG_enrichment.txt"), sep="\t", row.names=FALSE)
  
  return(list(GO=go_results, KEGG=kegg_results))
}

# Run enrichment analysis for all modules
go_kegg_results <- list()

for (mod in names(module_mappings_list)) {
  # Extract Entrez IDs for current module
  entrez_ids <- unique(module_mappings_list[[mod]]$entrezgene_id)
  
  # Run enrichment analysis if there are valid Entrez IDs
  if (length(na.omit(entrez_ids)) > 0) {
    go_kegg_results[[mod]] <- run_enrichment(entrez_ids, mod)
  }
}

# FILTER FOR SIGNIFICANT MODULES

# Define significance threshold
significance_threshold <- 0.01

# Create list to store only modules with significant enrichment
significant_modules <- list()

for (mod in names(go_kegg_results)) {
  go_results <- go_kegg_results[[mod]]$GO
  kegg_results <- go_kegg_results[[mod]]$KEGG
  
  # Check if module has significant terms in both GO and KEGG
  if (any(go_results$p.adjust < significance_threshold) & any(kegg_results$p.adjust < significance_threshold)) {
    significant_modules[[mod]] <- go_kegg_results[[mod]]
  }
}

print(names(significant_modules))

# Create dot plot for significant GO terms in magenta module
dotplot(significant_modules$magenta$GO)

# CREATE SUMMARY TABLES OF SIGNIFICANT RESULTS

# Initialize data frames for significant results
significant_GO_df <- data.frame()
significant_KEGG_df <- data.frame()

for (mod in names(go_kegg_results)) {
  go_results <- go_kegg_results[[mod]]$GO
  kegg_results <- go_kegg_results[[mod]]$KEGG
  
  # Filter for significant results
  go_significant <- go_results[go_results$p.adjust < significance_threshold, ]
  kegg_significant <- kegg_results[kegg_results$p.adjust < significance_threshold, ]
  
  # Add module information and combine
  if (nrow(go_significant) > 0) {
    go_significant$Module <- mod
    significant_GO_df <- rbind(significant_GO_df, go_significant)
  }
  
  if (nrow(kegg_significant) > 0) {
    kegg_significant$Module <- mod
    significant_KEGG_df <- rbind(significant_KEGG_df, kegg_significant)
  }
}

head(significant_GO_df)
head(significant_KEGG_df)

# Save significant results
write.csv(significant_GO_df, "significant_GO_results.csv", row.names = FALSE)
write.csv(significant_KEGG_df, "significant_KEGG_results.csv", row.names = FALSE)

############################################################################
# VISUALIZATION OF ENRICHMENT RESULTS
############################################################################

# Define color scheme for modules
module_color_map <- c(
  "blue" = "#0104ba",
  "cyan" = "#19e6ca",
  "darkorange2" = "#dc6e00",
  "green" = "#188003",
  "lightcoral" = "#ff4848",
  "yellow" = "#d8e41c",
  "magenta" = "#c600b4",
  "pink" = "#f5afee",
  "plum" = "#8b3d83",
  "salmon" = "#e58c6b",
  "tan" = "#fdd9c2",
  "turquoise" = "#3e8585"
)

# Convert GeneRatio from string to numeric for plotting
significant_GO_df <- significant_GO_df %>%
  mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))

significant_KEGG_df <- significant_KEGG_df %>%
  mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))

# Select top 30 GO terms for visualization
top_GO <- significant_GO_df %>%
  arrange(p.adjust) %>%
  head(30)

# Create bubble plot of GO term enrichment
ggplot(top_GO, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = -log10(p.adjust))) +
  geom_point(aes(size = Count)) +  # Point size represents number of genes
  scale_color_gradient(low = "blue", high = "red") +  # Color represents significance
  facet_grid(rows = vars(Module), scales = "free_y", space = "free_y") +  # Separate by module
  labs(title = "Top 30 GO Term Enrichments",
       x = "Gene Ratio",
       y = "GO Term",
       color = "-log10(adj.PValue)",
       size = "Gene Count") + 
  theme(axis.text = element_text(size = 12, colour = "black", family = "serif"),
        strip.text.y = element_text(angle = 0, face = "bold", hjust = 0, family = "serif"),
        strip.background = element_rect(fill = "darkgrey"),
        plot.title = element_text(size = 20, hjust = 0.5, family = "serif", face = "bold"),
        axis.title = element_text(size = 15, face = "bold", family = "serif", colour = "black"),
        legend.title = element_text(size = 15, face = "bold", family = "serif", colour = "black"),
        legend.text = element_text(size = 12, colour = "black", family = "serif"))

# Alternative visualization colored by module
ggplot(top_GO, aes(x = sort(GeneRatio, decreasing = T), y = reorder(Description, Count), color = Module)) +
  geom_point(aes(size = Count)) +
  scale_color_manual(values = module_color_map) +  # Use predefined module colors
  labs(title = "Top 30 GO Term Enrichments",
       x = "Gene Ratio",
       y = "GO Term",
       color = "Module",
       size = "Gene Count") + 
  theme(axis.text.y = element_text(size = 10))